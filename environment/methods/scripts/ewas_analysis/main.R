#!/usr/bin/env Rscript

# Load necessary libraries
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(EpiDISH)
library(limma)
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
    temp_file <- tempfile(fileext = ".parquet")
    raw_data <- base64enc::base64decode(what = base64_string)
    writeBin(raw_data, temp_file)
    df <- arrow::read_parquet(temp_file)
    unlink(temp_file)
    return(df)
  }, error = function(e) {
    stop(paste("Failed to decode Parquet data:", e$message))
  })
}

# Main eWAS function
episignature_call <- function(
  betas,
  pheno,
  covariates = c("Case_Cont", "Sex"),
  epi_type = "BC",
  padj_threshold = 0.01,
  deltaB_threshold = 0.3
) {
  
  message("************* Creating SummarizedExperiment for eWAS")
  SE <- SummarizedExperiment(
    assays = list(Beta = betas), 
    colData = pheno
  )
  
  # Filter to Case/Control samples
  message("************* Filtering to Case/Control samples")
  SE <- SE[, pheno$Case_Cont %in% c("Control", "Case")]
  pheno <- as.data.frame(colData(SE))
  pheno$Case_Cont <- factor(pheno$Case_Cont, levels = c("Control", "Case"))
  
  # Apply epi_type-specific processing
  if (epi_type == "BC" || epi_type == "bc") {
    message("************* Blood Cell deconvolution")
    covariates <- c(covariates, c("Bcell", "NK", "CD4T", "CD8T", "Mono", "Neu"))
    
    # Load reference
    data("centDHSbloodDMC.m", package = "EpiDISH")
    
    # Estimate cell proportions
    message("  Running EpiDISH...")
    ans <- EpiDISH::epidish(
      beta.m = assay(SE, "Beta"), 
      ref.m = centDHSbloodDMC.m, 
      method = "RPC"
    )$estF
    
    # Add cell proportions to pheno
    # Use rownames directly - they're already Sample IDs
    cell_props <- as.data.frame(ans) %>%
      dplyr::select(-Eosino) %>%
      `colnames<-`(c("Bcell", "NK", "CD4T", "CD8T", "Mono", "Neu"))
    
    # Combine by rownames (both have same sample IDs as rownames)
    pheno <- cbind(pheno, cell_props[rownames(pheno), ])
  }
  
  message("************* Converting to M values")
  methylation_data <- beta_to_m(assay(SE, "Beta"))
  
  message("************* Building model matrix")
  covariate_data <- pheno[, covariates, drop = FALSE]
  model_matrix <- model.matrix(~ ., data = covariate_data)
  
  message("************* Fitting linear model")
  fit <- limma::lmFit(object = methylation_data, design = model_matrix)
  fit <- limma::eBayes(fit = fit)
  
  message("************* Extracting top table")
  fit_top <- limma::topTable(fit, coef = "Case_ContCase", number = Inf) %>%
    dplyr::mutate(id = rownames(.))
  
  message("************* Calculating delta beta")
  deltaB <- tibble::tibble(
    id = rownames(SE),
    deltaBeta = rowMeans(assay(SE, "Beta")[, SE$Case_Cont == "Case"], na.rm = TRUE) - 
                rowMeans(assay(SE, "Beta")[, SE$Case_Cont == "Control"], na.rm = TRUE)
  )
  
  fit_top <- fit_top %>%
    dplyr::left_join(deltaB, by = "id")
  
  message(paste("************* Applying filters: padj <", padj_threshold, ", |deltaB| >=", deltaB_threshold))
  fit_top_flt <- fit_top %>%
    dplyr::filter(
      abs(deltaBeta) >= deltaB_threshold,
      adj.P.Val < padj_threshold
    )
  
  message(paste("************* Found", nrow(fit_top_flt), "significant CpGs out of", nrow(fit_top)))
  
  return(list(
    results = fit_top_flt,
    summary = list(
      n_cpgs_tested = nrow(fit_top),
      n_cpgs_significant = nrow(fit_top_flt),
      n_samples_case = sum(pheno$Case_Cont == "Case"),
      n_samples_control = sum(pheno$Case_Cont == "Control")
    )
  ))
}

# Main function to handle input from previous job
perform_ewas_analysis <- function(previous_result, covariates, epi_type, padj_threshold, deltaB_threshold) {
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
    cpg_names <- betas_df[[1]]
    betas_matrix <- as.matrix(betas_df[, -1])
    rownames(betas_matrix) <- cpg_names
    
    # Set pheno rownames from first column (Sample_ID or similar)
    sample_ids <- pheno_df[[1]]
    pheno_df <- pheno_df[, -1, drop = FALSE]
    rownames(pheno_df) <- sample_ids
    
    # Ensure column names match
    common_samples <- intersect(colnames(betas_matrix), rownames(pheno_df))
    if (length(common_samples) == 0) {
      stop("No common samples between betas and pheno")
    }
    
    message(paste("************* Using", length(common_samples), "common samples"))
    betas_matrix <- betas_matrix[, common_samples]
    pheno_df <- pheno_df[common_samples, , drop = FALSE]
    
    # Perform eWAS analysis
    ewas_result <- episignature_call(
      betas = betas_matrix,
      pheno = pheno_df,
      covariates = covariates,
      epi_type = epi_type,
      padj_threshold = padj_threshold,
      deltaB_threshold = deltaB_threshold
    )
    
    result <- list(
      status = "success",
      message = "eWAS analysis completed successfully",
      data = ewas_result,
      parameters = list(
        covariates = covariates,
        epi_type = epi_type,
        padj_threshold = padj_threshold,
        deltaB_threshold = deltaB_threshold
      )
    )
    
    return(result)
    
  }, error = function(e) {
    result <- list(
      status = "error",
      error = paste("Failed to perform eWAS analysis:", e$message)
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
  
  # Check if files exist
  if (!file.exists(input_file)) {
    result <- list(
      status = "error",
      message = paste("Error: Input file", input_file, "does not exist")
    )
    cat(toJSON(result, auto_unbox = TRUE))
    quit(status = 1)
  }
  
  if (!file.exists(params_file)) {
    result <- list(
      status = "error",
      message = paste("Error: Parameters file", params_file, "does not exist")
    )
    cat(toJSON(result, auto_unbox = TRUE))
    quit(status = 1)
  }
  
  # Load the previous result
  previous_result <- tryCatch({
    fromJSON(readLines(input_file, warn = FALSE))
  }, error = function(e) {
    result <- list(
      status = "error",
      message = paste("Error loading previous result:", e$message)
    )
    cat(toJSON(result, auto_unbox = TRUE))
    quit(status = 1)
  })
  
  # Load the parameters
  params <- tryCatch({
    fromJSON(readLines(params_file, warn = FALSE))
  }, error = function(e) {
    result <- list(
      status = "error",
      message = paste("Error loading parameters:", e$message)
    )
    cat(toJSON(result, auto_unbox = TRUE))
    quit(status = 1)
  })
  
  # Get parameters with defaults
  covariates <- if (!is.null(params$covariates)) params$covariates else c("Case_Cont", "Sex")
  epi_type <- if (!is.null(params$epi_type)) params$epi_type else "BC"
  padj_threshold <- if (!is.null(params$padj_threshold)) as.numeric(params$padj_threshold) else 0.01
  deltaB_threshold <- if (!is.null(params$deltaB_threshold)) as.numeric(params$deltaB_threshold) else 0.3
  
  # Perform eWAS analysis
  result <- perform_ewas_analysis(previous_result, covariates, epi_type, padj_threshold, deltaB_threshold)
  
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
    quit(status = 1)
  })
}

# Run main function
if (sys.nframe() == 0) {
  success <- main()
  quit(status = if(success) 0 else 1)
}

