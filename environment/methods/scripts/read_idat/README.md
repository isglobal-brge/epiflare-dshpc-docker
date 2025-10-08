# read_idat Method

## Description

Processes Illumina methylation array IDAT files using the minfi Bioconductor package. This method performs functional normalization, quality control filtering, and returns Beta values along with phenotype data.

## Input Requirements

The input must be a **tar.gz archive file** (`.tar.gz` or `.tgz`) containing:

1. **pheno.csv**: A CSV file with sample metadata containing at least:
   - `PID`: Patient/Sample ID
   - `array_id`: Array identifier
   - `Sex`: Sample sex (required for functional normalization)

2. **IDATs/**: A subfolder containing the raw IDAT files
   - Files should be named: `{PID}_{array_id}_Red.idat` and `{PID}_{array_id}_Grn.idat`

### Example Archive Structure
```
data_folder/
├── pheno.csv
└── IDATs/
    ├── Sample1_12345678_Red.idat
    ├── Sample1_12345678_Grn.idat
    ├── Sample2_87654321_Red.idat
    └── Sample2_87654321_Grn.idat
```

### Creating Deterministic Archives

For consistent file hashing (same content = same hash), you **must** use deterministic tar.gz compression:

```bash
tar --sort=name \
    --mtime='2000-01-01 00:00:00' \
    --owner=0 --group=0 \
    --numeric-owner \
    -czf data.tar.gz data_folder/
```

**Important**: This ensures that compressing the same content will **always produce the same hash**, which is critical for the file deduplication system.

#### Supported Formats
- `.tar.gz` or `.tgz` **only** (ZIP is not supported because it's non-deterministic)

## Processing Steps

1. **Load Data**: Reads IDAT files using minfi
2. **Detection P-values**: Calculates detection p-values for quality assessment
3. **Functional Normalization**: Applies `preprocessFunnorm` from minfi
4. **Sample Filtering**: Removes samples with mean detection p-value ≥ 0.01
5. **CpG Filtering**: Removes CpGs with detection p-value ≥ 0.01 in any sample
6. **Extract Results**: Returns Beta values matrix and phenotype data

## Output Files

The method generates two CSV files:

1. **Betas.csv**: Matrix of Beta values (methylation levels)
   - Rows: CpG probes
   - Columns: Samples
   - Values: Beta values (0-1, representing methylation proportion)

2. **pheno.csv**: Phenotype/sample metadata
   - Contains the processed sample information from the GRangesSet

## Quality Control

- **Sample QC**: Samples with poor detection (mean detP ≥ 0.01) are removed
- **Probe QC**: CpG probes failing detection threshold (detP ≥ 0.01) in any sample are removed
- **Normalization**: Functional normalization accounts for technical variation

## Parameters

This method has no configurable parameters. Quality control thresholds are fixed at:
- Detection p-value threshold: 0.01

## Output JSON

The method returns a JSON object with:
```json
{
  "status": "success",
  "message": "IDAT files processed successfully",
  "output_files": {
    "betas": "Betas.csv",
    "pheno": "pheno.csv"
  },
  "summary": {
    "samples_processed": 10,
    "cpgs_retained": 850000,
    "samples_filtered": 2
  }
}
```

## Dependencies

- minfi
- jsonlite
- readr
- dplyr

## Notes

- This method uses functional normalization which requires sex information
- Processing time depends on the number of samples and array type
- Large datasets will produce large output files
- All quality control metrics are fixed (not user-configurable)

## References

- Aryee MJ, et al. (2014). Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics.
- Fortin JP, et al. (2014). Functional normalization of 450k methylation array data improves replication in large cancer studies. Genome Biology.

