#!/usr/bin/env Rscript

# Load necessary libraries
library(minfi)
library(jsonlite)
library(readr)
library(dplyr)

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

# Main function to read and process IDAT files
read_idat <- function(archive_path) {
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
    
    # Save output files
    output_dir <- getwd()
    betas_file <- file.path(output_dir, "Betas.csv")
    pheno_file <- file.path(output_dir, "pheno.csv")
    
    message("************* Saving Betas matrix")
    write.csv(Betas, betas_file, row.names = TRUE)
    
    message("************* Saving pheno data")
    write.csv(pheno, pheno_file, row.names = TRUE)
    
    # Clean up temporary extraction directory
    message("************* Cleaning up temporary files")
    unlink(extract_dir, recursive = TRUE)
    
    result <- list(
      status = "success",
      message = "IDAT files processed successfully",
      output_files = list(
        betas = "Betas.csv",
        pheno = "pheno.csv"
      ),
      summary = list(
        samples_processed = ncol(Betas),
        cpgs_retained = nrow(Betas),
        samples_filtered = sum(!keep)
      )
    )
    
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
  
  if (length(args) < 2) {
    result <- list(
      status = "error",
      message = "Usage: Rscript main.R <input_file> <params_file>"
    )
    cat(toJSON(result, auto_unbox = TRUE))
    return(FALSE)
  }
  
  input_file <- args[1]
  params_file <- args[2]
  
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
  
  # Process the IDAT files
  result <- read_idat(input_file)
  
  # Add parameters applied to the result
  result$parameters_applied <- params
  
  # Output the result as JSON
  cat(toJSON(result, auto_unbox = TRUE, pretty = TRUE))
  return(result$status == "success")
}

# Run main function
if (sys.nframe() == 0) {
  success <- main()
  # Always exit with code 0, even if there are errors
  quit(status = 0)
}

