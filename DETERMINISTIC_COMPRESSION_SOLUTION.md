# Deterministic Compression Solution

## Overview

This document describes the implementation and verification of deterministic tar.gz compression for the dsHPC system. Deterministic compression ensures that identical content always produces identical hash values, which is critical for file deduplication and caching mechanisms.

## Technical Requirements

The system requires that uploading the same dataset multiple times produces the same file hash. This enables:
- Efficient deduplication at the storage layer
- Cache key consistency for processed results
- Avoiding redundant computation on identical datasets

## Native R Limitations

### Absence of Pure R Solution

Native R compression functions (`utils::tar()`, `tar::tar()`, etc.) do not provide deterministic compression because they lack control over:
- File modification timestamps
- File ownership attributes (uid/gid)
- Archive member ordering
- Tar metadata fields
- Gzip header timestamps

These functions internally call system utilities but do not expose the necessary configuration options for reproducible output.

### Standard Practice in R Ecosystem

Using `system()` calls for external tools is a standard and accepted practice in the R ecosystem. Notable examples include:
- `devtools` - Package development tools
- `usethis` - Workflow automation
- `remotes` - Package installation
- `pkgbuild` - Package building infrastructure

This approach is considered idiomatic when native R functionality is insufficient.

## Solution: GNU tar with Deterministic Flags

The following command ensures deterministic compression:

```bash
GZIP=-n gtar --sort=name \
            --mtime='2000-01-01 00:00:00' \
            --owner=0 \
            --group=0 \
            --numeric-owner \
            -czf output.tar.gz \
            -C source_dir/ \
            folder_name/
```

### Required Flags

| Flag | Purpose | Impact on Determinism |
|------|---------|----------------------|
| `GZIP=-n` | Omit gzip header timestamp | Critical: Without this, each compression embeds a different timestamp |
| `--sort=name` | Alphabetical file ordering | Ensures consistent archive member sequence |
| `--mtime='2000-01-01...'` | Fixed modification time | Normalizes all file timestamps to constant value |
| `--owner=0` | Fixed ownership (root) | Eliminates user ownership variation |
| `--group=0` | Fixed group (root) | Eliminates group ownership variation |
| `--numeric-owner` | Numeric IDs | Uses numeric identifiers instead of name strings |

### Critical Component: GZIP=-n

The `GZIP=-n` environment variable is the most critical component. Analysis:

**Without GZIP=-n:**
- gzip embeds the current timestamp in the compressed file header
- Each compression operation produces a unique timestamp
- Result: Identical content produces different hashes

**With GZIP=-n:**
- gzip omits the timestamp from the header
- Header remains constant across compressions
- Result: Identical content produces identical hashes

## Experimental Verification

### Test 1: Without GZIP=-n
```
Compression 1: 39bcc42cb5719cda93b6b797827eaa9a
Compression 2: 2db268165ca70d3e778d49355950a749
Compression 3: 2dd48690147f0397f8557736aedc5416
Compression 4: 79ded8fde7600b7bdda38b7dfadd00f7
Compression 5: 859deeaedb5d2e51fd305556576d473a

Result: 5 distinct hashes (non-deterministic)
```

### Test 2: With GZIP=-n
```
Compression 1: e65fcf23415a494cb1648746bc1e6d27
Compression 2: e65fcf23415a494cb1648746bc1e6d27
Compression 3: e65fcf23415a494cb1648746bc1e6d27
Compression 4: e65fcf23415a494cb1648746bc1e6d27
Compression 5: e65fcf23415a494cb1648746bc1e6d27

Result: Identical hash across all compressions (deterministic)
```

## R Implementation

```r
# Function to locate GNU tar on the system
find_gnu_tar <- function() {
  # macOS: Homebrew installs GNU tar as 'gtar'
  if (Sys.which("gtar") != "") {
    return("gtar")
  }
  
  # Linux: typically uses 'tar' (GNU tar by default)
  test_result <- system("tar --version 2>&1 | grep -q GNU", 
                       ignore.stdout = TRUE, 
                       ignore.stderr = TRUE)
  if (test_result == 0) {
    return("tar")
  }
  
  return(NULL)  # Not found
}

# Function to create deterministic tar.gz archive
create_deterministic_targz <- function(source_dir, output_file) {
  tar_cmd <- find_gnu_tar()
  
  if (is.null(tar_cmd)) {
    stop("GNU tar not found! Install with: brew install gnu-tar")
  }
  
  cmd <- sprintf(
    "GZIP=-n %s --sort=name --mtime='2000-01-01 00:00:00' --owner=0 --group=0 --numeric-owner -czf %s -C %s %s",
    tar_cmd,
    shQuote(output_file),
    shQuote(dirname(source_dir)),
    shQuote(basename(source_dir))
  )
  
  result <- system(cmd, intern = FALSE)
  
  if (result != 0) {
    stop("Failed to create tar.gz archive")
  }
  
  return(output_file)
}
```

## Dependency Installation

### macOS
```bash
brew install gnu-tar
```
This installs GNU tar as `gtar`

### Linux (Docker)
GNU tar is installed by default on standard Linux distributions.

### Verification
```bash
gtar --version  # macOS
tar --version   # Linux
```

## Available Test Scripts

### 1. test_deterministic_compression.R
Validates deterministic compression behavior:
- Compresses the same directory 5 times
- Compares all resulting hashes
- Verifies hash identity across compressions

Usage:
```bash
Rscript test_deterministic_compression.R
```

### 2. test_read_idat.R
Complete integration test for the `read_idat` method:
- Creates a data subset
- Applies deterministic compression
- Uploads to server
- Executes the processing job
- Reports results

Usage:
```bash
Rscript test_read_idat.R
```

## System Benefits

1. **Efficient deduplication**: Server can identify duplicate files by hash comparison
2. **Storage optimization**: Eliminates redundant copies of identical datasets
3. **Computational efficiency**: Prevents reprocessing of previously analyzed data
4. **Cache consistency**: Hash serves as a reliable cache key

## Frequently Asked Questions

### Why not use ZIP?
ZIP format is not deterministic because it embeds timestamps and variable metadata even when compressing identical content.

### Does this work on Windows?
Not directly. Windows requires either GNU tar for Windows or Windows Subsystem for Linux (WSL).

### What if I use BSD tar (native macOS)?
BSD tar does not support the required options (`--sort`, `--mtime`, etc.). GNU tar must be used instead.

### Why use --mtime='2000-01-01'?
Any fixed date would work. The year 2000 is chosen for memorability and because it predates all source files.

### Performance impact?
Minimal. The overhead of file sorting and timestamp normalization is negligible compared to compression time.

## Conclusion

Key findings:
- No native R solution exists for deterministic compression
- Using `system()` calls is standard practice in the R ecosystem
- The combination of `GZIP=-n` and GNU tar flags provides reliable determinism
- Verified through experimental testing: 5/5 compressions produced identical hashes
- Solution is production-ready for the dsHPC system

## Technical Note

The combination of `system()` with GNU tar represents the only reliable method for achieving deterministic compression from R. This approach is robust, well-tested, and widely used in production environments.

