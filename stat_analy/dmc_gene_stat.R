###############################################################################
# DMC-Gene Statistical Analysis
# 
# This script analyzes the relationship between differentially expressed genes 
# (DEGs) and differentially methylated cytosines (DMCs) in gene, promoter, 
# and upstream regions using Bayesian statistical modeling.
# 
# Main functions:
# - meth_in_region: Count methylation sites within genomic regions
# - get_DMC_in_regions: Extract DMC data for specific regions
# - test_Expr_DMC/Mval/Bval: Bayesian models for expression-methylation relationships
# 
# Author: [Xiaoyu Zhang]
# Date: [2025-08-08]
###############################################################################

# Load required libraries
library(rethinking)
library(tidyverse)
library(DESeq2)
library(methylKit)

# Configure parallel processing
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

# Configuration constants
NUM_SAMPLES <- 6
MIN_EXPRESSED_SAMPLES <- 3
NUM_CHAINS <- 4
NUM_CORES <- 4
SAMPLE_SIZE_SUBSET <- 1000
# ============================================================================
# DATA LOADING AND PREPROCESSING
# ============================================================================

# Load gene regions data (assumes gene_regions is available from methyl_by_regions.R)
if (!exists("gene_regions")) {
  stop("gene_regions object not found. Please run ../methyl_by_regions.R first")
}

# Validate gene_regions structure
required_cols <- c("ID", "chr", "start", "end")
if (!all(required_cols %in% colnames(gene_regions))) {
  stop(paste("gene_regions missing required columns:", 
             paste(setdiff(required_cols, colnames(gene_regions)), collapse = ", ")))
}

# Remove duplicate gene regions
gene_regions_uq <- gene_regions %>% 
  distinct(ID, .keep_all = TRUE)

# Load differential expression data
deg_data_file <- "./RData/deg.RData"
if (!file.exists(deg_data_file)) {
  stop(paste("Expression data file not found:", deg_data_file))
}
load(deg_data_file)

# Calculate DESeq2 scaling factors
samples <- data.frame(
  samples = colnames(deg_data[, 2:7]),
  strain = rep(c("EH1516", "EH217"), each = 3)
)

dds <- DESeqDataSetFromMatrix(
  countData = deg_data[, 2:7], 
  colData = samples, 
  design = ~1
)
dds <- estimateSizeFactors(dds)
scaling_factors <- log(sizeFactors(dds))

# Standardize column names
colnames(deg_data)[c(8, 9, 14, 15)] <- c("expr_base", "std", "logfc", "std.1")

# Merge gene regions with expression data
gene_data <- merge(gene_regions_uq, deg_data, by = "ID", no.dups = TRUE)

# Load methylation data
dmc_files <- c("./RData/dmc.RData", "./RData/dmc_list.RData")
missing_files <- dmc_files[!file.exists(dmc_files)]
if (length(missing_files) > 0) {
  stop(paste("Methylation data files not found:", paste(missing_files, collapse = ", ")))
}

load("./RData/dmc.RData")
load("./RData/dmc_list.RData")

# Validate methylation data structure
if (!exists("dmc_data_list") || !is.list(dmc_data_list)) {
  stop("dmc_data_list not found or not a list")
}

# Define strain coding
strain_vector <- rep(0:1, each = 3)  # 0: EH1516, 1: EH217
# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Calculate methylation statistics within genomic regions
#' 
#' @param gr Data frame with genomic regions (chr, start, end)
#' @return Data frame with methylation counts and levels
#' @details Optimized version with vectorized operations and error handling
meth_in_region <- function(gr) {
  if (nrow(gr) == 0) {
    warning("Empty genomic regions provided")
    return(data.frame())
  }
  
  nr <- nrow(gr)
  
  # Pre-allocate result matrix for better performance
  result_cols <- c("ncpg", "dmc_neg", "dmc_pos", "b1", "b2", "dmc_b1", "dmc_b2")
  m_count <- data.frame(matrix(0, nrow = nr, ncol = length(result_cols)))
  colnames(m_count) <- result_cols
  
  # Helper function to calculate methylation levels
  calc_meth_level <- function(data, num_cols, cov_cols) {
    if (nrow(data) == 0) return(0)
    total_cs <- rowSums(data[, num_cols, drop = FALSE])
    total_cov <- rowSums(data[, cov_cols, drop = FALSE])
    if (sum(total_cov) == 0) return(0)
    return(sum(total_cs) / sum(total_cov))
  }
  
  # Process each region
  for (i in seq_len(nr)) {
    chr <- as.character(gr[i, 1])
    start <- gr[i, 2]
    end <- gr[i, 3]
    
    meth_d <- dmc_data_list[[chr]]
    
    if (!is.null(meth_d) && nrow(meth_d) > 0) {
      # Vectorized region filtering
      region_mask <- meth_d$start >= start & meth_d$end <= end
      m_in_g <- meth_d[region_mask, ]
      
      if (nrow(m_in_g) > 0) {
        nc <- nrow(m_in_g)
        
        # Vectorized DMC counting
        dmc_mask <- m_in_g$dmc
        meth_diff <- m_in_g$m2 - m_in_g$m1
        n_neg <- sum(dmc_mask & meth_diff < 0, na.rm = TRUE)
        n_pos <- sum(dmc_mask & meth_diff > 0, na.rm = TRUE)
        
        # Calculate methylation levels for all CpGs
        eh1516_cols <- c("numCs1", "numCs2")
        eh1516_cov_cols <- c("coverage1", "coverage2")
        eh217_cols <- c("numCs3", "numCs4", "numCs5")
        eh217_cov_cols <- c("coverage3", "coverage4", "coverage5")
        
        b1 <- calc_meth_level(m_in_g, eh1516_cols, eh1516_cov_cols)
        b2 <- calc_meth_level(m_in_g, eh217_cols, eh217_cov_cols)
        
        # Calculate methylation levels for DMCs only
        dmc_b1 <- 0
        dmc_b2 <- 0
        if (sum(dmc_mask) > 0) {
          dmc_data <- m_in_g[dmc_mask, ]
          dmc_b1 <- calc_meth_level(dmc_data, eh1516_cols, eh1516_cov_cols)
          dmc_b2 <- calc_meth_level(dmc_data, eh217_cols, eh217_cov_cols)
        }
        
        m_count[i, ] <- c(nc, n_neg, n_pos, b1, b2, dmc_b1, dmc_b2)
      }
    }
  }
  
  return(m_count)
}

#' Extract DMC data for genomic regions
#' 
#' @param regions Data frame with genomic regions
#' @return Data frame with region data merged with methylation statistics
get_DMC_in_regions <- function(regions) {
  if (!requireNamespace("methylKit", quietly = TRUE)) {
    stop("methylKit package is required but not installed")
  }
  
  # Merge regions with expression data
  merged_data <- merge(regions, deg_data, by = "ID", no.dups = TRUE)
  
  # Calculate methylation statistics for each region
  meth_stats <- meth_in_region(merged_data[, 2:4])
  combined_data <- cbind(merged_data, meth_stats)
  
  # Filter for regions with CpG sites and sufficient expression
  filtered_data <- combined_data %>%
    filter(
      ncpg > 0,
      rowSums(.[, 8:13] > 0) >= MIN_EXPRESSED_SAMPLES
    )
  
  return(filtered_data)
}

# ============================================================================
# DATA PROCESSING WORKFLOWS
# ============================================================================

# Get methylation data in gene regions
gene_data_dmc <- get_DMC_in_regions(gene_regions_uq)

#' Reshape expression data for Bayesian modeling
#' 
#' @param expr_data Data frame with expression data
#' @param nsample Number of samples
#' @param strain_coding Vector indicating strain for each sample
#' @return Reshaped data frame suitable for modeling
reshape_expression_data <- function(expr_data, nsample = NUM_SAMPLES, strain_coding = strain_vector) {
  if (ncol(expr_data) < nsample + 1) {
    stop("Expression data has insufficient columns for specified number of samples")
  }
  
  # More efficient reshape using tidyr approach
  gene_ids <- expr_data[, 1]
  expr_matrix <- expr_data[, 2:(nsample + 1)]
  
  # Create long format data
  reshaped_data <- expand.grid(
    gene_idx = seq_len(nrow(expr_data)),
    sample = seq_len(nsample)
  )
  
  # Add expression values, strain info, and scaling factors
  reshaped_data$gid <- rep(gene_ids, each = nsample)
  reshaped_data$expr <- as.vector(t(expr_matrix))
  reshaped_data$strain <- rep(strain_coding, nrow(expr_data))
  reshaped_data$sf <- rep(scaling_factors, nrow(expr_data))
  
  # Reorder columns for consistency
  reshaped_data <- reshaped_data[, c("gid", "expr", "sample", "strain", "sf", "gene_idx")]
  colnames(reshaped_data)[ncol(reshaped_data)] <- "id"
  
  return(reshaped_data)
}

#' Create data list for Bayesian modeling
#' 
#' @param reshaped_data Reshaped expression data
#' @param region_data Original region data with methylation info
#' @param model_type Type of model ("dmc", "mval", "bval")
#' @return List suitable for rethinking::ulam
create_model_data <- function(reshaped_data, region_data, model_type = "dmc") {
  e_bar <- log(median(reshaped_data$expr))
  n_genes <- max(reshaped_data$id)
  n_samples <- nrow(reshaped_data)
  
  base_data <- list(
    G = reshaped_data$id,
    E = reshaped_data$expr,
    S = reshaped_data$sample,
    T = reshaped_data$strain,
    sf = reshaped_data$sf,
    e_bar = rep(e_bar, n_samples),
    ng = n_genes
  )
  
  # Add model-specific variables
  if (model_type == "dmc") {
    base_data$XC <- normalize(log(rep(region_data$ncpg, NUM_SAMPLES)))
    base_data$XDP <- normalize(log(rep(region_data$dmc_pos + 1, NUM_SAMPLES)))
    base_data$XDN <- normalize(log(rep(region_data$dmc_neg + 1, NUM_SAMPLES)))
  } else if (model_type == "mval") {
    base_data$dm <- c(
      rep(region_data$m1 - region_data$m1, 3),  # EH1516 samples
      rep(region_data$m2 - region_data$m1, 3)   # EH217 samples
    )
  } else if (model_type == "bval") {
    base_data$db <- c(
      rep(region_data$b1 - region_data$b1, 3),  # EH1516 samples
      rep(region_data$b2 - region_data$b1, 3)   # EH217 samples
    )
  }
  
  return(base_data)
}

#' Test expression-DMC relationship using Bayesian modeling
#' 
#' @param data Input data frame
#' @param model_type Type of model ("dmc", "mval", "bval")
#' @return Fitted Bayesian model
test_expression_methylation <- function(data, model_type = "dmc") {
  # Prepare data
  data$id <- seq_len(nrow(data))
  expr_cols <- if (ncol(data) >= 30) c(1, 8:13, 2:7, 26:30) else c(1, 8:13, 2:7)
  expr_data <- data[, expr_cols]
  
  # Reshape data
  reshaped_data <- reshape_expression_data(expr_data)
  
  # Add methylation variables
  reshaped_data$ncpg <- rep(data$ncpg, NUM_SAMPLES)
  reshaped_data$dmc_pos <- rep(data$dmc_pos, NUM_SAMPLES)
  reshaped_data$dmc_neg <- rep(data$dmc_neg, NUM_SAMPLES)
  if ("DE" %in% colnames(data)) {
    reshaped_data$DE <- rep(data$DE, NUM_SAMPLES)
  }
  
  # Create model data
  model_data <- create_model_data(reshaped_data, data, model_type)
  
  # Define and fit model based on type
  if (model_type == "dmc") {
    model <- ulam(
      alist(
        E ~ dgampois(lambda, phi),
        log(lambda) <- e_bar + f[S] + e[G] + bDP * XDP * T + bDN * XDN * T,
        vector[6]: f ~ normal(0, 0.2),
        vector[ng]: e ~ normal(0, 3),
        bDP ~ normal(0, 1.5),
        bDN ~ normal(0, 1.5),
        phi ~ dexp(1)
      ), 
      data = model_data, 
      chains = NUM_CHAINS, 
      cores = NUM_CORES
    )
  } else if (model_type == "mval") {
    model <- ulam(
      alist(
        E ~ dgampois(lambda, phi),
        log(lambda) <- e_bar + f[S] + e[G] + bM * dm,
        vector[6]: f ~ normal(0, 0.1),
        vector[ng]: e ~ normal(0, 3),
        bM ~ normal(0, 1.5),
        phi ~ dexp(1)
      ), 
      data = model_data, 
      chains = NUM_CHAINS, 
      cores = NUM_CORES
    )
  } else if (model_type == "bval") {
    model <- ulam(
      alist(
        E ~ dgampois(lambda, phi),
        log(lambda) <- e_bar + f[S] + e[G] + bB * db,
        vector[6]: f ~ normal(0, 0.1),
        vector[ng]: e ~ normal(0, 3),
        bB ~ normal(0, 1.5),
        phi ~ dexp(1)
      ), 
      data = model_data, 
      chains = NUM_CHAINS, 
      cores = NUM_CORES
    )
  }
  
  return(model)
}

# Wrapper functions for backward compatibility
test_Expr_DMC <- function(d) {
  return(test_expression_methylation(d, "dmc"))
}

test_Expr_Mval <- function(d) {
  return(test_expression_methylation(d, "mval"))
}

test_Expr_Bval <- function(d) {
  return(test_expression_methylation(d, "bval"))
}

# ============================================================================
# MAIN ANALYSIS WORKFLOW
# ============================================================================

# Fit basic DMC models
message("Fitting DMC models for gene regions...")
gene_models <- list(
  all_genes = test_Expr_DMC(gene_data_dmc[, 1:29]),
  dmc_genes_only = test_Expr_DMC(
    gene_data_dmc[, 1:29] %>% filter(dmc_pos > 0 | dmc_neg > 0)
  )
)

# Fit methylation value models  
message("Fitting methylation value models...")
mval_models <- list(
  dmc_genes_mval = test_Expr_Mval(
    gene_data_dmc %>% filter(dmc_neg > 0 | dmc_pos > 0)
  ),
  dmc_genes_bval = test_Expr_Bval(
    gene_data_dmc %>% filter((dmc_neg + dmc_pos) > 0)
  )
)

# Display model summaries
message("Model summaries:")
cat("\n=== Basic DMC Model (All Genes) ===\n")
precis(gene_models$all_genes, depth = 2, prob = 0.95,
       pars = c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]", "f[6]"))

cat("\n=== DMC Model (DMC Genes Only) ===\n")
precis(gene_models$dmc_genes_only)


# ============================================================================
# REGIONAL ANALYSIS (PROMOTERS AND UPSTREAM REGIONS)
# ============================================================================

#' Analyze DMCs in specific genomic regions
#' 
#' @param region_name Name of the region for output
#' @param region_data Region data from getPromoterRegions
#' @return List of fitted models
analyze_region <- function(region_name, region_data) {
  message(paste("Analyzing", region_name, "regions..."))
  
  # Clean region data
  if (ncol(region_data) > 29 && !is.null(region_data[, 7])) {
    region_data[, 7] <- NULL  # Remove problematic column if present
  }
  
  # Filter for regions with DMCs
  dmc_regions <- region_data %>% filter(dmc_neg > 0 | dmc_pos > 0)
  
  if (nrow(dmc_regions) == 0) {
    warning(paste("No DMCs found in", region_name, "regions"))
    return(NULL)
  }
  
  # Fit models
  models <- list(
    all_regions = test_Expr_DMC(region_data[, 1:29]),
    dmc_regions_only = test_Expr_DMC(dmc_regions[, 1:29]),
    dmc_mval_model = test_Expr_Mval(dmc_regions)
  )
  
  return(models)
}

# Get scaffold information
if (!exists("scaffold_ranges")) {
  scaffold_ranges <- read_scaffold_gr("data/Ehux_genome.fasta.len")
}

# Define regions to analyze
region_configs <- list(
  promoters = list(up = 1000, down = 1000, name = "promoter"),
  upstream_2kb = list(up = 2000, down = 0, name = "upstream 2kb"),
  upstream_1kb = list(up = 1000, down = 0, name = "upstream 1kb")
)

# Analyze each region type
regional_models <- list()

for (config_name in names(region_configs)) {
  config <- region_configs[[config_name]]
  
  # Get region coordinates
  regions <- getPromoterRegions(
    up = config$up, 
    down = config$down, 
    bed_file = "data/Ehux_genbank.bed"
  ) %>% 
    trimRegionByScaffolds()
  
  # Get DMC data for regions
  region_dmc_data <- get_DMC_in_regions(regions)
  
  # Analyze regions
  regional_models[[config_name]] <- analyze_region(config$name, region_dmc_data)
}

# Display regional analysis summaries
message("\nRegional Analysis Complete!")
for (region_name in names(regional_models)) {
  if (!is.null(regional_models[[region_name]])) {
    cat(paste("\n=== Models for", region_name, "===\n"))
    cat("All regions model summary:\n")
    print(summary(regional_models[[region_name]]$all_regions))
  }
}
