#!/usr/bin/env Rscript

# Load necessary libraries
library(minfi)
library(jsonlite)
library(readr)
library(dplyr)
library(arrow)
library(base64enc)

# Function to extract compressed file
extract_archive <- function(archive_path, extract_dir) {
  # Only accept tar.gz or tgz
  if (!grepl("\\.tar\\.gz$|\\.tgz$", archive_path, ignore.case = TRUE)) {
    stop("Unsupported archive format. Only .tar.gz or .tgz files are accepted (for deterministic compression)")
  }
  
  # Extract tar.gz
  message("************* Extracting tar.gz archive")
  untar(archive_path, exdir = extract_dir)
  
  # Find the extracted folder (usually there's one top-level folder)
  extracted_items <- list.files(extract_dir, full.names = TRUE)
  
  if (length(extracted_items) == 0) {
    stop("Archive is empty or extraction failed")
  }
  
  # If there's only one item and it's a directory, use it
  if (length(extracted_items) == 1 && dir.exists(extracted_items[1])) {
    return(extracted_items[1])
  }
  
  # Otherwise, return the extraction directory itself
  return(extract_dir)
}

# Function to convert dataframe to Parquet and encode as base64
dataframe_to_parquet_base64 <- function(df) {
  tryCatch({
    # Create a temporary file for Parquet
    temp_file <- tempfile(fileext = ".parquet")
    
    # Write dataframe to Parquet with compression
    arrow::write_parquet(
      df,
      temp_file,
      compression = "snappy",  # Fast compression, good balance
      use_dictionary = TRUE    # Better compression for repeated values
    )
    
    # Encode the Parquet file to base64
    parquet_base64 <- base64enc::base64encode(temp_file)
    
    # Get file size for metadata
    file_size <- file.info(temp_file)$size
    
    # Clean up
    unlink(temp_file)
    
    return(list(
      data = parquet_base64,
      format = "parquet",
      compression = "snappy",
      size_bytes = file_size,
      rows = nrow(df),
      cols = ncol(df)
    ))
  }, error = function(e) {
    # Fallback to JSON if Parquet fails
    warning(paste("Failed to convert to Parquet:", e$message))
    return(list(
      data = df,
      format = "json",
      error = e$message
    ))
  })
}

# Main function to read and process IDAT files
read_idat <- function(archive_path, output_format = "hybrid") {
  tryCatch({
    # Create temporary directory for extraction
    temp_dir <- tempdir()
    extract_dir <- file.path(temp_dir, paste0("idat_extract_", as.integer(Sys.time())))
    dir.create(extract_dir, recursive = TRUE)
    
    message(paste("************* Extracting archive:", basename(archive_path)))
    
    # Extract the archive
    data_folder <- extract_archive(archive_path, extract_dir)
    
    message(paste("************* Data extracted to:", data_folder))
    
    # Check if pheno.csv exists
    pheno_path <- file.path(data_folder, "pheno.csv")
    if (!file.exists(pheno_path)) {
      stop("pheno.csv not found in the extracted folder")
    }
    
    # Check if IDATs folder exists
    idats_folder <- file.path(data_folder, "IDATs")
    if (!dir.exists(idats_folder)) {
      stop("IDATs/ subfolder not found in the extracted folder")
    }
    
    # Read targets
    targets <- readr::read_csv(pheno_path, show_col_types = FALSE) %>%
      mutate(
        "Basename" = paste0(idats_folder, "/", PID, "_", array_id)
      )
    
    message("************* Load RGset")
    RGset <- minfi::read.metharray.exp(targets = targets)
    sampleNames(RGset) <- targets$PID
    
    detP <- minfi::detectionP(RGset)
    
    message("************* Load GRset (preprocessFunnorm)")
    GRset <- minfi::preprocessFunnorm(RGset, sex = targets$Sex)
    
    # Partial processing
    # First, Sample filtering according to detP.
    keep <- colMeans(detP) < 0.01
    GRset <- GRset[, keep]
    detP <- detP[, keep]
    
    # Now, CpG filtering.
    if (!(identical(rownames(detP), rownames(GRset)))) {
      message("Not identical CpGs")
      common_cpg <- intersect(rownames(detP), rownames(GRset))
      detP <- detP[common_cpg, ]
      GRset <- GRset[common_cpg, ]
    }
    keep <- rowSums(detP < 0.01) == ncol(GRset)
    GRset <- GRset[keep, ]
    
    # Store in separate tables.
    pheno <- as.data.frame(colData(GRset))
    Betas <- getBeta(GRset)
    
    # Convert Betas matrix to data frame with CpG names as a column
    message("************* Preparing Betas data for output")
    betas_df <- as.data.frame(Betas)
    betas_df$CpG <- rownames(Betas)
    # Move CpG column to first position
    betas_df <- betas_df[, c("CpG", setdiff(names(betas_df), "CpG"))]
    
    # Add sample IDs as column for pheno data
    message("************* Preparing pheno data for output")
    pheno$Sample_ID <- rownames(pheno)
    # Move Sample_ID to first position
    pheno <- pheno[, c("Sample_ID", setdiff(names(pheno), "Sample_ID"))]
    
    # Clean up temporary extraction directory
    message("************* Cleaning up temporary files")
    unlink(extract_dir, recursive = TRUE)
    
    # Decide output format based on parameter
    if (output_format == "parquet" || output_format == "hybrid") {
      message("************* Converting data to Parquet format")
      
      # Convert both dataframes to Parquet
      betas_parquet <- dataframe_to_parquet_base64(betas_df)
      pheno_parquet <- dataframe_to_parquet_base64(pheno)
      
      # Create result with Parquet data
      result <- list(
        status = "success",
        message = "IDAT files processed successfully",
        format = ifelse(
          betas_parquet$format == "parquet" && pheno_parquet$format == "parquet",
          "parquet",
          "mixed"
        ),
        data = list(
          betas = betas_parquet,
          pheno = pheno_parquet
        ),
        summary = list(
          samples_processed = ncol(Betas),
          cpgs_retained = nrow(Betas),
          # samples_filtered refers to CpGs filtered out based on detection p-values
          # keep <- rowSums(detP < 0.01) == ncol(GRset) filters CpGs where all samples have detP < 0.01
          # sum(!keep) counts how many CpGs were removed due to poor detection p-values
          samples_filtered = sum(!keep),
          data_size = list(
            betas = paste0(round(betas_parquet$size_bytes / 1024 / 1024, 2), " MB"),
            pheno = paste0(round(pheno_parquet$size_bytes / 1024, 2), " KB")
          )
        )
      )
      
      message(sprintf("************* Betas Parquet size: %.2f MB", betas_parquet$size_bytes / 1024 / 1024))
      message(sprintf("************* Pheno Parquet size: %.2f KB", pheno_parquet$size_bytes / 1024))
      
    } else {
      # Legacy JSON format
      result <- list(
        status = "success",
        message = "IDAT files processed successfully",
        format = "json",
        data = list(
          betas = betas_df,
          pheno = pheno
        ),
        summary = list(
          samples_processed = ncol(Betas),
          cpgs_retained = nrow(Betas),
          samples_filtered = sum(!keep)
        )
      )
    }
    
    return(result)
    
  }, error = function(e) {
    result <- list(
      status = "error",
      error = paste("Failed to process IDAT files:", e$message)
    )
    return(result)
  })
}

# Main entry point for the script
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 3) {
    result <- list(
      status = "error",
      message = "Usage: Rscript main.R <input_file> <metadata_file> <params_file>"
    )
    cat(toJSON(result, auto_unbox = TRUE))
    quit(status = 1)
  }
  
  input_file <- args[1]
  metadata_file <- args[2]
  params_file <- args[3]
  
  # Load metadata (contains workspace_dir and file paths)
  metadata <- tryCatch({
    fromJSON(readLines(metadata_file, warn = FALSE))
  }, error = function(e) {
    result <- list(
      status = "error",
      message = paste("Error loading metadata:", e$message)
    )
    cat(toJSON(result, auto_unbox = TRUE))
    quit(status = 1)
  })
  
  # Check if the input file exists
  if (!file.exists(input_file)) {
    result <- list(
      status = "error",
      message = paste("Error: Input file", input_file, "does not exist")
    )
    cat(toJSON(result, auto_unbox = TRUE))
    return(FALSE)
  }
  
  # Check if the params file exists
  if (!file.exists(params_file)) {
    result <- list(
      status = "error",
      message = paste("Error: Parameters file", params_file, "does not exist")
    )
    cat(toJSON(result, auto_unbox = TRUE))
    return(FALSE)
  }
  
  # Load the parameters (even if empty, we still read it for consistency)
  tryCatch({
    params <- fromJSON(readLines(params_file))
  }, error = function(e) {
    result <- list(
      status = "error",
      message = paste("Error loading parameters:", e$message)
    )
    cat(toJSON(result, auto_unbox = TRUE))
    return(FALSE)
  })
  
  # Check for output format parameter (default to hybrid for optimal size)
  output_format <- ifelse(
    !is.null(params$output_format),
    params$output_format,
    "hybrid"  # Default to hybrid format (Parquet for data, JSON for metadata)
  )
  
  # Process the IDAT files
  result <- read_idat(input_file, output_format)
  
  # Add parameters applied to the result
  result$parameters_applied <- params
  
  # Output the result as JSON
  # For large outputs, write to stdout efficiently
  tryCatch({
    # Convert to JSON without pretty printing for large outputs
    # Pretty printing can cause issues with very large strings
    json_output <- toJSON(result, auto_unbox = TRUE, pretty = FALSE)
    
    # Write the JSON output
    # Use writeLines for better handling of large strings
    writeLines(json_output, con = stdout())
    
    return(result$status == "success")
  }, error = function(e) {
    # If JSON serialization fails, output an error
    error_result <- list(
      status = "error",
      message = paste("Failed to serialize output:", e$message)
    )
    cat(toJSON(error_result, auto_unbox = TRUE))
    return(FALSE)
  })
}

# Run main function
if (sys.nframe() == 0) {
  success <- main()
  # Exit with appropriate code: 0 for success, 1 for error
  quit(status = if(success) 0 else 1)
}

