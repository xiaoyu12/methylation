#library(Gviz)
#library(ggbio)
library(ggrepel)

# Add a column of gene IDs to the GRanges or data.frame using a transcript ->
# gene mapping table.
addMapID <- function(ranges, mapping) {
  ranges$ID <- mapping[as.character(ranges$name), "ID"]
  ranges <- as.data.frame(ranges)
  ranges <- ranges[!is.na(ranges$ID), ]           # drop entries without mapping
  ranges <- ranges[!duplicated(ranges$ID), ]       # collapse alternative isoforms
  ranges <- as(ranges, "GRanges")
  return(ranges)
}

# -----------------------------------------------------------------------------
# Correlation helpers ---------------------------------------------------------
# -----------------------------------------------------------------------------

# Compute a sample-by-sample methylation correlation matrix from a methylBase.
compute_methyl_correlation <- function(meth) {
  if (!inherits(meth, "methylBase")) {
    stop("'meth' must be a methylBase object")
  }
  dat <- getData(meth)
  coverage_idx <- slot(meth, "coverage.index")
  numCs_idx <- slot(meth, "numCs.index")
  coverage <- dat[, coverage_idx]
  numCs <- dat[, numCs_idx]
  # guard against division by zero
  beta <- numCs / pmax(coverage, 1)
  colnames(beta) <- slot(meth, "sample.ids")
  stats::cor(beta, use = "pairwise.complete.obs")
}

# Backwards-compatible alias using historical naming style
getMethCorrelation <- compute_methyl_correlation


# Analyse DMRs overlapped with gene features and expression profiles.
# DMR can be either a GRanges or a data.frame coercible to GRanges.
overlap_dmr_gene_features <- function(DMR, features, gene_exps, mapping,
                                      selected_ids = NULL, ignore.strand = TRUE,
                                      plot_annotation = TRUE) {
  DMR.gr <- as(DMR, "GRanges")
  DMR.Ann <- annotateWithGeneParts(DMR.gr, features, intersect.chr = TRUE,
                                   ignore.strand = ignore.strand)
  if (isTRUE(plot_annotation)) {
    genomation::plotTargetAnnotation(DMR.Ann)
  }

  DMR.TSS <- getAssociationWithTSS(DMR.Ann)
  DMR.df <- if (is.data.frame(DMR)) DMR else as.data.frame(DMR)
  rows <- DMR.TSS$target.row
  DMR.TSS <- cbind(DMR.df[rows, , drop = FALSE], DMR.TSS)

  transcript_ids <- as.character(DMR.TSS$feature.name)
  if (!"ID" %in% colnames(mapping)) {
    stop("'mapping' must contain a column named 'ID'")
  }
  DMR.TSS$ID <- mapping[transcript_ids, "ID"]

  # Support both limma (logFC) and DESeq2 (LFC) outputs
  logfc_col <- intersect(c("logFC", "LFC"), colnames(gene_exps))
  if (length(logfc_col) == 0) {
    stop("gene_exps must contain either a 'logFC' or 'LFC' column")
  }
  DMR.TSS$logFC <- gene_exps[as.character(DMR.TSS$ID), logfc_col[1]]

  DMR.TSS.sel <- NULL
  if (!is.null(selected_ids)) {
    ids <- intersect(DMR.TSS$ID, selected_ids)
    DMR.TSS.sel <- DMR.TSS[DMR.TSS$ID %in% ids, , drop = FALSE]
  }

  list(ann = DMR.Ann, tss = DMR.TSS, sel = DMR.TSS.sel)
}

# Backwards-compatible wrappers ------------------------------------------------

overlap_DMR_GeneFeatures <- function(DMR, features, gene_exps,
                                     selected_ids = NULL, ignore.strand = TRUE,
                                     mapping = get0("mapping", ifnotfound = NULL),
                                     plot_annotation = TRUE) {
  if (is.null(mapping)) {
    stop("Provide a 'mapping' data.frame or ensure a global object named 'mapping' exists.")
  }
  overlap_dmr_gene_features(DMR, features, gene_exps, mapping,
                            selected_ids = selected_ids,
                            ignore.strand = ignore.strand,
                            plot_annotation = plot_annotation)
}

overlap_DMR_GeneFeatures2 <- function(DMR, features, gene_exps, mapping,
                                      selected_ids = NULL, ignore.strand = TRUE,
                                      plot_annotation = TRUE) {
  overlap_dmr_gene_features(DMR, features, gene_exps, mapping,
                            selected_ids = selected_ids,
                            ignore.strand = ignore.strand,
                            plot_annotation = plot_annotation)
}

# Calculate methylation m values for methylBase object
calc_M_values <- function (meth, nsamples = length(slot(meth, "sample.ids"))) {
  dat <- getData(meth)
  result <- dat[, seq_len(4)]
  numCs_idx <- slot(meth, "numCs.index")
  numTs_idx <- slot(meth, "numTs.index")
  ns <- min(nsamples, length(numCs_idx))
  for (i in seq_len(ns)) {
    colname <- paste0("m", i)
    result[[colname]] <- log2((dat[, numCs_idx[i]] + 1) / (dat[, numTs_idx[i]] + 1))
  }
  return(result)
}

# Calculate methylation beta values for methylBase object
calc_Beta_values <- function (meth, nsamples = length(slot(meth, "sample.ids"))) {
  dat <- getData(meth)
  result <- dat[, seq_len(4)]
  coverage_idx <- slot(meth, "coverage.index")
  numCs_idx <- slot(meth, "numCs.index")
  ns <- min(nsamples, length(numCs_idx))
  for (i in seq_len(ns)) {
    colname <- paste0("beta", i)
    coverage <- dat[, coverage_idx[i]]
    result[[colname]] <- ifelse(coverage == 0, NA_real_, dat[, numCs_idx[i]] / coverage)
  }
  return(result)
}
# Convert a region methylRaw list into a data.frame
# methylraw:  Large methylRawList 
# regions: data frame of genomic regions
methyl_to_data_frame <- function(methylraw, regions) {
  n <- length(methylraw)  # length of the list
  methyl_list <- as.list(1:n)
  for (i in 1:n) {
    t <- merge(getData(methylraw[[i]]), regions) # merge data frames according to "chr", "start", "end" and "strand"
    t <- t[!duplicated(t$ID), ]     # Remove duplicated ID's in the data frame
    t <- t[t$coverage >= 100, ]
    t$beta <- t$numCs / t$coverage
    t$m <- log2((t$numCs + 1) / (t$numTs + 1))  # Add one to avoid dividing by 0
    rownames(t) <- t$ID
    methyl_list[[i]] <- t
  }
  
  # Find shared IDs of data frames in methyl_list
  ids <- methyl_list[[1]]$ID
  for (i in 2:n) {
    ids <- intersect(ids, methyl_list[[i]]$ID)
  }
  # Get m values of the first sample
  data <- as.data.frame(methyl_list[[1]][as.character(ids), "m"])
  for (i in 2:n) {
    data <- cbind(data, methyl_list[[i]][as.character(ids), "m"])
  }
  colnames(data) <- paste0("m", 1:n)
  return(data)
}
  
# Get methylation values for gene annotated features: e.g. promoters, exons ...
getFeatureMethyl <- function(m.obj, g.features, s.names, mapping = NULL,
                             lo.count = 100, id_column = "ID") {
  nsamples <- length(m.obj)
  if (length(s.names) != nsamples) {
    stop("Length of 's.names' must equal the number of samples in 'm.obj'")
  }

  feature_df <- as.data.frame(g.features)
  colnames(feature_df)[1] <- "chr"
  if (!(id_column %in% colnames(feature_df))) {
    if (is.null(mapping)) {
      stop("g.features must contain column '", id_column,
           "' or a 'mapping' data.frame must be supplied")
    }
    feature_df[[id_column]] <- mapping[as.character(feature_df$name), "ID"]
  }
  feature_df <- feature_df[!is.na(feature_df[[id_column]]), ]
  feature_df <- feature_df[!duplicated(feature_df[[id_column]]), ]

  promobj <- regionCounts(m.obj, g.features)

  sample_tables <- vector("list", nsamples)
  names(sample_tables) <- s.names
  for (i in seq_len(nsamples)) {
    merged <- merge(
      getData(promobj[[i]]),
      feature_df,
      by.x = c("chr", "start", "end", "strand"),
      by.y = c("chr", "start", "end", "strand")
    )
    if (nrow(merged) == 0) {
      sample_tables[[i]] <- merged
      next
    }
    merged <- merged[merged$coverage >= lo.count, , drop = FALSE]
    if (nrow(merged) == 0) {
      sample_tables[[i]] <- merged
      next
    }
    merged$beta <- merged$numCs / merged$coverage
    merged$m <- log2((merged$numCs + 1) / (merged$numTs + 1))
    merged <- merged[!duplicated(merged[[id_column]]), , drop = FALSE]
    rownames(merged) <- merged[[id_column]]
    sample_tables[[i]] <- merged
  }

  id_lists <- lapply(sample_tables, function(df) df[[id_column]])
  id_lists <- id_lists[sapply(id_lists, length) > 0]
  if (length(id_lists) == 0) {
    empty <- matrix(nrow = 0, ncol = nsamples)
    colnames(empty) <- s.names
    return(list(coverage = empty, beta = empty, m = empty))
  }
  ids <- Reduce(intersect, id_lists)
  if (length(ids) == 0) {
    empty <- matrix(nrow = 0, ncol = nsamples)
    colnames(empty) <- s.names
    return(list(coverage = empty, beta = empty, m = empty))
  }

  build_matrix <- function(col_name) {
    mat <- matrix(nrow = length(ids), ncol = nsamples)
    for (i in seq_len(nsamples)) {
      df <- sample_tables[[i]]
      if (nrow(df) == 0) {
        mat[, i] <- NA_real_
      } else {
        mat[, i] <- df[ids, col_name, drop = TRUE]
      }
    }
    colnames(mat) <- s.names
    rownames(mat) <- ids
    mat
  }

  coverage_mat <- build_matrix("coverage")
  beta_mat <- build_matrix("beta")
  m_mat <- build_matrix("m")

  list(coverage = coverage_mat, beta = beta_mat, m = m_mat)
}

# get methylation data in a genomics range
getMethInGRange <- function(meth, gr) {
  meth.gr <- as(meth, "GRanges")
  hits <- !is.na(findOverlaps(meth.gr, gr, select = "first"))
  return (meth[hits, ])
}

# Read transcript features for a selected set of genes defined in gene_list_file.
read_selected_gene_features <- function(gene_list_file, genome_bed,
                                        up_flank = 1000, down_flank = 1000,
                                        rna_column = "rna",
                                        remove_unusual = FALSE,
                                        unique_prom = FALSE) {
  genes <- read.table(gene_list_file, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE)
  if (!(rna_column %in% colnames(genes))) {
    stop("Column '", rna_column, "' is not present in ", gene_list_file)
  }
  transcripts <- genes[[rna_column]]
  features <- readTranscriptFeatures(genome_bed,
                                     remove.unusual = remove_unusual,
                                     unique.prom = unique_prom,
                                     up.flank = up_flank,
                                     down.flank = down_flank)
  for (slot_name in names(features)) {
    features[[slot_name]] <- features[[slot_name]][
      features[[slot_name]]$name %in% transcripts
    ]
  }
  features
}

# Backwards-compatible wrapper using historic name and default genome bed path.
readSelectedGeneFeatures <- function(gene_list_file, up_flank = 1000,
                                     down_flank = 1000,
                                     genome_bed = file.path("..", "v2", "Ehux_genbank.bed"),
                                     ...) {
  read_selected_gene_features(gene_list_file,
                              genome_bed = genome_bed,
                              up_flank = up_flank,
                              down_flank = down_flank,
                              ...)
}

##
plot_DMR_Gviz <- function(dmrs, dmrTSS, id, methDiff = myDiff.DSS) {
  # All DMR's associated with the gene of id
  dmr_id <- dmrTSS[dmrTSS$ID == id, ]
  # Get GRanges for the DMR's
  dmr_id.gr <- as(dmrs[dmr_id$target.row, ], "GRanges")
  strand(dmr_id.gr) <- "*"
  
  # Get exon structure of gene of id
  rna <- rownames(mapping[mapping$ID == id, ])
  exons <- gene.parts$exons[gene.parts$exons$name == rna, ]
  grlist <- GRangesList("exons" = exons[, 2], "dmr" = dmr_id.gr[, 1])
  gr <- unlist(grlist)
  start <- min(start(range(gr)))
  end <- max(end(range(gr)))
  chr <- as.character(seqnames(gr))[1]   # convert factor to string
  
  
  gtrack <- GenomeAxisTrack()
  #meth_id <- getMethInGRange(methDiff, range(gr))
  meth_id <- getMethInGRange(methDiff, range(dmr_id.gr))
  meth_id.gr <- as(meth_id, "GRanges")
  meth_id.dt <- DataTrack(range=meth_id.gr, data=meth_id.gr$meth.diff, legend = TRUE, start = start, end = end, type = "h", chromosome = chr, name = "Methyl Diff")
  exons.track <- AnnotationTrack(exons, name = "Exons", chromosome = chr, start = start, end = end, id = "exons")
  plotTracks(list(meth_id.dt, gtrack, exons.track), main = paste("Gene", id, chr), sizes=c(1, 0.5, 0.3), from = start, to = end)
}

##
# Plot DMR's associated with a gene of id
plotOneDMR <- function(dmrs, dmrTSS, id, methDiff = myDiff.DSS) {
  dmr_id <- dmrTSS[dmrTSS$ID == id, ]
  dmr_id.gr <- as(dmrs[dmr_id$target.row, ], "GRanges")
  gr <- range(dmr_id.gr)
  strand(gr) <- "*"
  methhits <- getMethInGRange(methDiff, gr)
  matplot(methhits$start, methhits$meth.diff, type="o", pch = 20,col = 'blue',
          xlab=paste("Genomic location", seqnames(gr), sep = " "), 
          ylab="Methylation difference" )
  abline(h=c(-10,0,10),lty=2)
}

plot_DMRs <- function(dmrs) {
  plot(dmrs$meth.diff, dmrs$logFC, type="p", cex = 1, pch = 20, col= "blue", xlab = "Methylation Difference", ylab = "Log Fold Change")
  text(dmrs$meth.diff, dmrs$logFC + 0.25, labels = dmrs$ID, cex = 0.5)
}

# Plot methylation data for a gene
plot_Methyl_Grange <- function(methData, grang) {
  # 
  methhits <- getMethInGRange(methData, grang)
  matplot(methhits$start, methhits$meth.diff, type="o", pch = 20,col = 'blue',
          xlab=paste("Genomic location", methhits$chr[1]), 
          ylab="Methylation difference" )
  abline(h=c(-10,0,10),lty=2)
}

######################################################################
# New modular helpers for readability and efficiency
######################################################################

# Create directories if missing
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# Load key inputs (mapping, gene features, expression table)
load_inputs <- function(config) {
  if (is.null(config$mapping_file) || !file.exists(config$mapping_file)) {
    stop(paste0("Missing mapping_file: ", config$mapping_file))
  }
  if (is.null(config$genome_bed) || !file.exists(config$genome_bed)) {
    stop(paste0("Missing genome_bed: ", config$genome_bed))
  }
  mapping <- read.table(config$mapping_file, row.names = 1)
  colnames(mapping) <- c("gene", "ID")

  gene_parts <- readTranscriptFeatures(config$genome_bed, remove.unusual = FALSE, unique.prom = FALSE)
  gene_parts_up1000 <- readTranscriptFeatures(config$genome_bed, remove.unusual = FALSE,
                                              up.flank = 1000, down.flank = 0, unique.prom = FALSE)
  gene_parts_dn1000 <- readTranscriptFeatures(config$genome_bed, remove.unusual = FALSE,
                                              up.flank = 0, down.flank = 1000, unique.prom = FALSE)

  expr <- NULL
  if (!is.null(config$expr_file) && nzchar(config$expr_file) && file.exists(config$expr_file)) {
    expr <- read.csv(config$expr_file, header = TRUE, row.names = 1, sep = "\t")
  }
  return(list(mapping = mapping,
              gene_parts = gene_parts,
              gene_parts_up1000 = gene_parts_up1000,
              gene_parts_dn1000 = gene_parts_dn1000,
              expr = expr))
}

# Build methylation objects with simple caching
build_methyl_objects <- function(config) {
  ensure_dir(config$out_dir)
  cache_dir <- ifelse(is.null(config$cache_dir), file.path(config$out_dir, "cache"), config$cache_dir)
  ensure_dir(cache_dir)

  meth_rds <- file.path(cache_dir, "meth.rds")
  meth10_rds <- file.path(cache_dir, "meth10.rds")

  if (file.exists(meth_rds)) {
    meth <- readRDS(meth_rds)
  } else {
    methobj <- read_Bismark_coverage(config$data_dir)
    meth <- methylKit::unite(methobj)
    saveRDS(meth, meth_rds)
  }

  if (file.exists(meth10_rds)) {
    meth10 <- readRDS(meth10_rds)
  } else {
    # Recreate methobj if not in memory; read again (cheap I/O compared to compute)
    methobj <- read_Bismark_coverage(config$data_dir)
    meth10 <- unite(filterByCoverage(methobj, lo.count = config$min_cov))
    saveRDS(meth10, meth10_rds)
  }

  # Optional CHG via bismark list
  bismark <- NULL
  if (!isFALSE(config$load_chg)) {
    bismark <- read_Bismark_CpG_CHG(base_dir = config$data_dir)
  }
  return(list(meth = meth, meth10 = meth10, bismark = bismark))
}

# Simple QC plotting wrappers
plot_correlation_heatmap <- function(meth, out_file, meta = NULL) {
  cor_mat <- compute_methyl_correlation(meth)
  png(filename = out_file, width = 2400, height = 1600, res = 600)
  on.exit(dev.off(), add = TRUE)
  pheatmap::pheatmap(cor_mat, annotation = meta, fontsize = 8)
}

plot_pca <- function(meth, out_file, adj_lim = c(0.4, 0.1)) {
  png(filename = out_file, width = 4800, height = 3200, res = 600)
  on.exit(dev.off(), add = TRUE)
  methylKit::PCASamples(meth, adj.lim = adj_lim, scale = FALSE)
}

# Run methylKit differential methylation
run_methylkit_dm <- function(meth, diff_percent = 25, qvalue = 0.01, mc.cores = 1) {
  diff_obj <- calculateDiffMeth(meth, mc.cores = mc.cores)
  sig <- getMethylDiff(diff_obj, difference = diff_percent, qvalue = qvalue)
  return(list(all = diff_obj, sig = sig))
}
