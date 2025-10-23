#!/usr/bin/env Rscript

# Load necessary libraries
library(SummarizedExperiment)
library(dplyr)
library(arrow)
library(base64enc)
library(jsonlite)

# Function to convert beta values to M values
beta_to_m <- function(beta_matrix, eps = 1e-6) {
  beta <- pmin(pmax(beta_matrix, eps), 1 - eps)
  log2(beta / (1 - beta))
}

# Function to decode Parquet data from base64
decode_parquet_base64 <- function(base64_string) {
  tryCatch({
    # Create temporary file
    temp_file <- tempfile(fileext = ".parquet")
    
    # Decode base64 string to raw bytes
    raw_data <- base64enc::base64decode(what = base64_string)
    
    # Write raw bytes to file
    writeBin(raw_data, temp_file)
    
    # Read Parquet file
    df <- arrow::read_parquet(temp_file)
    
    # Clean up
    unlink(temp_file)
    
    return(df)
  }, error = function(e) {
    stop(paste("Failed to decode Parquet data:", e$message))
  })
}

# Main PCA function
pca_topVariable <- function(
    betas,
    pheno,
    top_CpG = 10000,
    array = "EPIC"
) {
  
  message("************* Creating SummarizedExperiment")
  SE <- SummarizedExperiment(
    assays = SimpleList(
      "Beta" = as.matrix(betas),
      "M" = as.matrix(beta_to_m(betas))
    ),
    colData = pheno
  )

  message(paste("************* Loading", array, "annotations"))
  if (array == "450k") {
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylation450kmanifest)
    annotations <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else if (array == "EPIC") {
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    library(IlluminaHumanMethylationEPICmanifest)
    annotations <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else if (array == "EPICv2") {
    library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    library(IlluminaHumanMethylationEPICv2manifest)
    annotations <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  } else {
    stop(paste("Unsupported array type:", array))
  }
  
  message("************* Filtering to common CpGs")
  common_cpg <- intersect(rownames(SE), rownames(annotations))
  annotations <- annotations[common_cpg, ]
  SE <- SE[common_cpg, ]
  rowData(SE) <- annotations

  message("************* Filtering sex chromosomes")
  SE <- SE[!rowData(SE)$chr %in% c("chrX", "chrY"), ]

  message(paste("************* Selecting top", top_CpG, "most variable CpGs"))
  variance_betas <- apply(X = assay(SE, "M"), MARGIN = 1, var)
  top_var <- sort(variance_betas, decreasing = TRUE)[1:min(top_CpG, length(variance_betas))]

  message("************* Running PCA")
  pca_res <- prcomp(t(assay(SE, "M")[names(top_var), ]))

  message("************* Calculating variance explained")
  var_explained <- round(((pca_res$sdev**2) * 100) / (sum(pca_res$sdev ** 2)), 2)
  names(var_explained) <- colnames(pca_res$x)
  # Convert to named list for better JSON serialization
  var_explained_list <- as.list(var_explained)
  
  message("************* Preparing PCA results with metadata")
  # Get PCA coordinates
  pca_coords <- as.data.frame(pca_res$x)
  
  # Get sample metadata (colData from SummarizedExperiment)
  sample_metadata <- as.data.frame(colData(SE))
  
  # Ensure row names match
  sample_metadata <- sample_metadata[rownames(pca_coords), , drop = FALSE]
  
  return(list(
    pca_coordinates = pca_coords,
    sample_metadata = sample_metadata,
    variance_explained = var_explained_list,
    n_samples = ncol(SE),
    n_cpgs_used = length(top_var),
    n_cpgs_filtered = nrow(SE)
  ))
}

# Main function to handle input from previous job
perform_pca_analysis <- function(previous_result, top_CpG = 10000, array = "EPIC") {
  tryCatch({
    # Extract betas and pheno from previous result
    if (is.null(previous_result$data$betas) || is.null(previous_result$data$pheno)) {
      stop("Previous result must contain 'data$betas' and 'data$pheno'")
    }
    
    message("************* Extracting data from previous result")
    
    # Check if data is in Parquet format (base64 encoded)
    if (!is.null(previous_result$data$betas$format) && 
        previous_result$data$betas$format == "parquet") {
      message("************* Decoding Parquet betas data")
      betas_df <- decode_parquet_base64(previous_result$data$betas$data)
    } else {
      message("************* Using JSON betas data")
      betas_df <- as.data.frame(previous_result$data$betas)
    }
    
    if (!is.null(previous_result$data$pheno$format) && 
        previous_result$data$pheno$format == "parquet") {
      message("************* Decoding Parquet pheno data")
      pheno_df <- decode_parquet_base64(previous_result$data$pheno$data)
    } else {
      message("************* Using JSON pheno data")
      pheno_df <- as.data.frame(previous_result$data$pheno)
    }
    
    message(paste("************* Betas dimensions:", nrow(betas_df), "x", ncol(betas_df)))
    message(paste("************* Pheno dimensions:", nrow(pheno_df), "x", ncol(pheno_df)))
    
    # Convert betas_df to matrix (first column is CpG names)
    cpg_names <- betas_df[[1]]  # First column (CpG or similar)
    betas_matrix <- as.matrix(betas_df[, -1])
    rownames(betas_matrix) <- cpg_names
    
    # Set pheno rownames from first column (Sample_ID or similar)
    sample_ids <- pheno_df[[1]]
    pheno_df <- pheno_df[, -1]
    rownames(pheno_df) <- sample_ids
    
    # Ensure column names match
    common_samples <- intersect(colnames(betas_matrix), rownames(pheno_df))
    if (length(common_samples) == 0) {
      stop("No common samples between betas and pheno")
    }
    
    message(paste("************* Using", length(common_samples), "common samples"))
    betas_matrix <- betas_matrix[, common_samples]
    pheno_df <- pheno_df[common_samples, , drop = FALSE]
    
    # Perform PCA analysis
    pca_result <- pca_topVariable(
      betas = betas_matrix,
      pheno = pheno_df,
      top_CpG = top_CpG,
      array = array
    )
    
    result <- list(
      status = "success",
      message = "PCA analysis completed successfully",
      data = pca_result,
      parameters = list(
        top_CpG = top_CpG,
        array = array
      )
    )
    
    return(result)
    
  }, error = function(e) {
    result <- list(
      status = "error",
      error = paste("Failed to perform PCA analysis:", e$message)
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
  
  # Check if the input file exists (this is the result from previous job)
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
  
  # Load the previous result
  tryCatch({
    previous_result <- fromJSON(readLines(input_file, warn = FALSE))
  }, error = function(e) {
    result <- list(
      status = "error",
      message = paste("Error loading previous result:", e$message)
    )
    cat(toJSON(result, auto_unbox = TRUE))
    return(FALSE)
  })
  
  # Load the parameters
  tryCatch({
    params <- fromJSON(readLines(params_file, warn = FALSE))
  }, error = function(e) {
    result <- list(
      status = "error",
      message = paste("Error loading parameters:", e$message)
    )
    cat(toJSON(result, auto_unbox = TRUE))
    return(FALSE)
  })
  
  # Get parameters with defaults
  top_CpG <- ifelse(!is.null(params$top_CpG), as.integer(params$top_CpG), 10000)
  array <- ifelse(!is.null(params$array), params$array, "EPIC")
  
  # Perform PCA analysis
  result <- perform_pca_analysis(previous_result, top_CpG, array)
  
  # Add parameters applied to the result
  result$parameters_applied <- params
  
  # Output the result as JSON
  tryCatch({
    json_output <- toJSON(result, auto_unbox = TRUE, pretty = FALSE)
    writeLines(json_output, con = stdout())
    return(result$status == "success")
  }, error = function(e) {
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
  quit(status = 0)
}

