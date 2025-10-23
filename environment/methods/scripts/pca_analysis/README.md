# PCA Analysis Method

## Description
Performs Principal Component Analysis (PCA) on methylation beta values from the `read_idat` method output. This method is designed to be chained after `read_idat` using meta-jobs.

## Features
- Automatically decodes Parquet or JSON format from previous job
- Converts beta values to M-values for PCA
- Filters out sex chromosome CpGs (chrX, chrY)
- Selects most variable CpGs for analysis
- Generates publication-quality PCA plot
- Returns variance explained by each PC
- Returns PCA coordinates for downstream analysis

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `top_CpG` | integer | 10000 | Number of most variable CpGs to use for PCA |
| `array` | string | "EPIC" | Array type for annotations: "450k", "EPIC", or "EPICv2" |

## Input
This method expects the output from `read_idat` as input, containing:
- `data.betas`: Beta values matrix (CpGs × Samples)
- `data.pheno`: Sample metadata with `Case_Cont` and `Sex` columns

The input can be in either JSON or Parquet format (automatically detected).

## Output

Returns a JSON object with:

```json
{
  "status": "success",
  "message": "PCA analysis completed successfully",
  "data": {
    "plot": "<base64-encoded PNG>",
    "variance_explained": {
      "PC1": 25.5,
      "PC2": 15.3,
      ...
    },
    "pca_coordinates": {
      "PC1": [...],
      "PC2": [...],
      ...
    },
    "n_samples": 48,
    "n_cpgs_used": 10000,
    "n_cpgs_filtered": 850000
  },
  "parameters": {
    "top_CpG": 10000,
    "array": "EPIC"
  }
}
```

## Usage with Meta-Jobs

### Example: Chain read_idat → pca_analysis

```r
library(dsHPC)

# Create resource client
resource <- list(
  name = "test_dshpc",
  url = "http://localhost:8001",
  format = "dshpc.api",
  secret = "your_api_key"
)

client <- HPCResourceClient$new(resource)

# Create API config
api_config <- dsHPC:::create_api_config(
  base_url = "http://localhost",
  port = 8001,
  api_key = "your_api_key",
  auth_header = "X-API-Key",
  auth_prefix = ""
)

# Upload IDAT archive
idat_content <- readBin("your_data.tar.gz", "raw", file.info("your_data.tar.gz")$size)
client$uploadFile(idat_content, "your_data.tar.gz")

# Define processing chain
chain <- list(
  list(
    method_name = "read_idat",
    parameters = list(
      output_format = "hybrid"  # Use Parquet for efficiency
    )
  ),
  list(
    method_name = "pca_analysis",
    parameters = list(
      top_CpG = 10000,
      array = "EPIC"
    )
  )
)

# Execute chain
result <- execute_processing_chain(
  config = api_config,
  content = idat_content,
  method_chain = chain,
  upload_filename = "your_data.tar.gz",
  timeout = 600  # 10 minutes
)

# Extract and display plot
plot_png <- base64enc::base64decode(result$data$plot, output = "plot.png")
system("open plot.png")  # macOS
# Or use: browseURL("plot.png")

# View variance explained
print(result$data$variance_explained)

# Access PCA coordinates
pca_coords <- as.data.frame(result$data$pca_coordinates)
head(pca_coords)
```

## Dependencies

### R Packages
- SummarizedExperiment
- ggplot2
- dplyr
- tibble
- arrow (for Parquet support)
- base64enc
- jsonlite
- Annotation packages (loaded dynamically):
  - IlluminaHumanMethylation450kmanifest
  - IlluminaHumanMethylation450kanno.ilmn12.hg19
  - IlluminaHumanMethylationEPICmanifest
  - IlluminaHumanMethylationEPICanno.ilm10b4.hg19
  - IlluminaHumanMethylationEPICv2manifest
  - IlluminaHumanMethylationEPICv2anno.20a1.hg38

## Plot Details

The generated PCA plot includes:
- PC1 vs PC2 scatter plot
- Points colored by `Case_Cont` (case/control status)
- Points shaped by `Sex`
- 99% confidence ellipses for each group
- Variance explained percentages in axis labels
- Clean, minimal theme

## Performance Notes

- Filtering sex chromosomes reduces dimensionality
- Using top variable CpGs (default 10000) improves PCA quality
- M-values (log-odds) are more suitable for PCA than beta values
- Parquet input format provides faster data loading

## Error Handling

The method handles:
- Missing or malformed input data
- Incompatible array types
- Sample/CpG mismatches between betas and pheno
- Both JSON and Parquet input formats

## Caching

When used in meta-jobs, results are cached based on:
- Previous job output (file hash)
- Method parameters (`top_CpG`, `array`)
- Method version and runtime environment

Identical chains will reuse cached results instantly.

