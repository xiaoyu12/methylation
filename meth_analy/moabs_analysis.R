library(methylKit)
library(genomation)
library(stringr)
library(ggplot2)

source("meth_analysis_funcs.R")

config <- list(
  cpg_dmr_file = "dmr_M2_EH1516.G.bed_vs_EH217.G.bed.bed",
  cpg_dmc_file = "dmc_M2_EH1516.G.bed_vs_EH217.G.bed.bed",
  chg_dmr_file = "dmr_M2_EH1516.HG.bed_vs_EH217.HG.bed.bed",
  chg_dmc_file = "dmc_M2_EH1516.HG.bed_vs_EH217.HG.bed.bed",
  alt_cpg_dmr_file = "dmr_M2_moabs1.3.4.bed",
  biominer_gene_file = "../v2/biominer_genes.txt",
  de_annotation_file = "../../gene_expression/Ehux_1516_v_217/Jan2019/output/DE_E1516_vs_E217_FDR0.05.annot.txt",
  output_dir = "moabs_graphs",
  plot_res = 600
)

if (!exists("ensure_dir")) {
  ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }
}

ensure_dir(config$output_dir)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

read_moabs_bed <- function(filename) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required to import MOABS files.")
  }
  if (!file.exists(filename)) {
    stop("File not found: ", filename)
  }
  gr <- rtracklayer::import(filename, format = "bed")
  name_tokens <- stringr::str_sub(gr$name, 5L)
  value_pairs <- stringr::str_split_fixed(name_tokens, "vs", 2)
  meth_diff <- suppressWarnings(as.numeric(value_pairs[, 1]) - as.numeric(value_pairs[, 2]))
  gr$meth.diff <- meth_diff
  strand(gr) <- "*"
  gr
}

readMOABS_DMR_bedfile <- read_moabs_bed

dmr_to_df <- function(gr) {
  df <- as.data.frame(gr)
  if ("seqnames" %in% names(df)) {
    names(df)[names(df) == "seqnames"] <- "chr"
  }
  df
}

select_logfc_column <- function(expr_table) {
  if ("logFC" %in% colnames(expr_table)) return("logFC")
  if ("LFC" %in% colnames(expr_table)) return("LFC")
  stop("Expression coefficient table must contain either 'logFC' or 'LFC'.")
}

extract_overlap_genes <- function(query, subject, mapping) {
  hits <- findOverlaps(query, subject, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(subject[FALSE, ])
  }
  genes_hit <- subject[unique(to(hits)), ]
  genes_hit$ID <- mapping[genes_hit$name, "ID"]
  genes_hit
}

build_overlap_table <- function(dmr, feature, mapping, expr_table, logfc_col, min_width = NULL) {
  if (!is.null(min_width)) {
    dmr <- dmr[width(dmr) >= min_width]
  }
  hits <- findOverlaps(dmr, feature, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(data.frame())
  }
  dmr_idx <- from(hits)
  feat_idx <- to(hits)
  dmr_df <- as.data.frame(dmr[dmr_idx])
  feature_df <- as.data.frame(feature[feat_idx])
  gene_ids <- mapping[feature_df$name, "ID"]
  log_fc <- expr_table[as.character(gene_ids), logfc_col]
  data.frame(
    dmr.seqnames = dmr_df$seqnames,
    dmr.start = dmr_df$start,
    dmr.end = dmr_df$end,
    dmr.width = dmr_df$width,
    meth.diff = dmr_df$meth.diff,
    feature.name = feature_df$name,
    ID = gene_ids,
    logFC = log_fc,
    stringsAsFactors = FALSE
  )
}

subset_dmr_tss <- function(dmrTSS, ids) {
  if (length(ids) == 0) {
    return(dmrTSS[FALSE, , drop = FALSE])
  }
  subset(dmrTSS, ID %in% ids)
}

getMethDiff_for_Ids <- function(methDiff, diffTSS, ids) { # methDiff kept for signature compatibility
  subset_dmr_tss(diffTSS, ids)
}

getDiffExprCoef_for_Ids <- function(ids) {
  expr.coef[as.character(ids), log_fc_column]
}

plotMethyl_and_LFC <- function(methDiffAll, diffTSS, ids, label_points = TRUE) {
  subset_df <- getMethDiff_for_Ids(methDiffAll, diffTSS, ids)
  if (nrow(subset_df) == 0) {
    return(subset_df)
  }
  subset_df$lfc <- expr.coef[as.character(subset_df$ID), log_fc_column]
  plot(subset_df$meth.diff, subset_df$lfc, pch = 20, col = "blue",
       xlab = "Methylation Difference", ylab = "Log Fold Change")
  if (isTRUE(label_points)) {
    text(subset_df$meth.diff, subset_df$lfc + 0.25, labels = subset_df$ID, cex = 0.5)
  }
  subset_df
}

aggregate_DMR <- function(dmrTSS) {
  if (nrow(dmrTSS) == 0) {
    return(data.frame(ID = character(), dist.to.feature = character(), count = integer(),
                      meth.diff = character(), meth.diff.mean = numeric(), logFC = numeric(),
                      stringsAsFactors = FALSE))
  }
  dist_concat <- aggregate(dist.to.feature ~ ID, data = dmrTSS, paste, collapse = ";")
  dist_count <- aggregate(dist.to.feature ~ ID, data = dmrTSS, length)
  meth_concat <- aggregate(meth.diff ~ ID, data = dmrTSS, paste, collapse = ";")
  meth_mean <- aggregate(meth.diff ~ ID, data = dmrTSS, mean)
  merged <- Reduce(function(x, y) merge(x, y, by = "ID"),
                   list(dist_concat, dist_count, meth_concat, meth_mean))
  colnames(merged) <- c("ID", "dist.to.feature", "count", "meth.diff", "meth.diff.mean")
  merged$logFC <- expr.coef[as.character(merged$ID), log_fc_column]
  merged
}

plot_Mean_Meth_Diff_and_LFC <- function(dmr_mean_meth_diff) {
  if (nrow(dmr_mean_meth_diff) == 0) {
    return(dmr_mean_meth_diff)
  }
  valid <- subset(dmr_mean_meth_diff, !is.na(logFC) & !is.na(meth.diff.mean))
  if (nrow(valid) == 0) {
    return(dmr_mean_meth_diff)
  }
  x_limits <- range(c(valid$meth.diff.mean, -50, 50), na.rm = TRUE)
  x_limits <- c(min(-50, x_limits[1] - 5), max(50, x_limits[2] + 5))
  y_limits <- range(c(valid$logFC, -2, 2), na.rm = TRUE)
  y_limits <- c(min(-2, y_limits[1] - 0.25), max(2, y_limits[2] + 0.5))
  plot(valid$meth.diff.mean, valid$logFC, pch = 20, cex = 0.8, col = "blue",
       xlab = "Avg DMR methylation diff", ylab = "Log Fold Change",
       xlim = x_limits, ylim = y_limits)
  abline(h = 0, col = "red")
  abline(v = 0, col = "red")
  abline(h = c(1, -1), lty = 2, col = "gray")
  dmr_mean_meth_diff
}

plot_DMR_and_LFC <- function(methDiffAll, dmrTSS, ids, label_points = TRUE) {
  subset_df <- plotMethyl_and_LFC(methDiffAll, dmrTSS, ids, label_points = label_points)
  mean_df <- aggregate_DMR(subset_df)
  plot_Mean_Meth_Diff_and_LFC(mean_df)
  mean_df
}

build_dmr_tss_table <- function(dmr, annotation, mapping, expr_table, logfc_col) {
  assoc <- getAssociationWithTSS(annotation)
  dmr_df <- as.data.frame(dmr)
  combined <- cbind(dmr_df[assoc$target.row, , drop = FALSE], assoc)
  combined$ID <- mapping[combined$feature.name, "ID"]
  combined$logFC <- expr_table[as.character(combined$ID), logfc_col]
  combined
}

# -----------------------------------------------------------------------------
# Validate required global objects
# -----------------------------------------------------------------------------

required_objects <- c("mapping", "expr.coef", "genes", "gene.parts")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop("Missing required objects in environment: ", paste(missing_objects, collapse = ", "))
}

log_fc_column <- select_logfc_column(expr.coef)

if (!all(c("Gene", "UTR") %in% names(genes))) {
  stop("The 'genes' object must contain 'Gene' and 'UTR' components.")
}

# -----------------------------------------------------------------------------
# Load MOABS DMR/DMC data
# -----------------------------------------------------------------------------

moabs.DMR <- readMOABS_DMR_bedfile(config$cpg_dmr_file)
moabs.DMC <- readMOABS_DMR_bedfile(config$cpg_dmc_file)

moabs.134 <- readMOABS_DMR_bedfile(config$alt_cpg_dmr_file)
overlap_idx <- findOverlaps(moabs.134, moabs.DMR)
if (length(overlap_idx) > 0) {
  message("Replacing CpG DMRs with high-confidence overlaps (", length(unique(from(overlap_idx))), " regions).")
  moabs.DMR <- moabs.134[unique(from(overlap_idx)), ]
}
moabs.DMR.df <- dmr_to_df(moabs.DMR)

moabs.HG.DMC <- readMOABS_DMR_bedfile(config$chg_dmc_file)
moabs.HG.DMR <- readMOABS_DMR_bedfile(config$chg_dmr_file)
moabs.HG.DMR.df <- dmr_to_df(moabs.HG.DMR)
moabs.HG.DMR.l100 <- subset(moabs.HG.DMR.df, width >= 100)

if (nrow(moabs.HG.DMR.l100) >= 2 && exists("methCHG_DSS_mincov10")) {
  try(showOneDMR(moabs.HG.DMR.l100[2, ], methCHG_DSS_mincov10$BSobj), silent = TRUE)
}

moabs.all <- c(moabs.DMR, moabs.HG.DMR)

# -----------------------------------------------------------------------------
# Overlaps with gene features
# -----------------------------------------------------------------------------

DMGs <- extract_overlap_genes(moabs.all, genes$Gene, mapping)
DMGs.CpG <- extract_overlap_genes(moabs.DMR, genes$Gene, mapping)
DMGs.CHG <- extract_overlap_genes(moabs.HG.DMR, genes$Gene, mapping)

DMR.ov.genes <- build_overlap_table(moabs.DMR, genes$Gene, mapping, expr.coef, log_fc_column, min_width = 100)
if (nrow(DMR.ov.genes) > 0) {
  png(filename = file.path(config$output_dir, "DMGs_Methdiff_LFG.png"), res = config$plot_res,
      width = 2400, height = 2400)
  print(ggplot(DMR.ov.genes) +
          geom_point(aes(x = meth.diff, y = logFC), colour = "blue", shape = 20) +
          theme_light())
  dev.off()
}

if (exists("up2000", inherits = FALSE) && exists("scaffolds.gr", inherits = FALSE)) {
  up2000 <- intersect(up2000, scaffolds.gr)
  DMPs <- extract_overlap_genes(moabs.all, up2000, mapping)
  DMPs.CpG <- extract_overlap_genes(moabs.DMR, up2000, mapping)
  DMPs.CHG <- extract_overlap_genes(moabs.HG.DMR, up2000, mapping)
} else {
  warning("Objects 'up2000' and/or 'scaffolds.gr' not found; promoter overlap analysis skipped.")
  DMPs <- DMPs.CpG <- DMPs.CHG <- genes$Gene[FALSE, ]
}

if (exists("up1000") && exists("down1000")) {
  annot.up1000 <- annotateWithFeature(moabs.DMR, up1000)
  genomation::plotTargetAnnotation(annot.up1000)
  annot.down1000 <- annotateWithFeature(moabs.DMR, down1000)
  genomation::plotTargetAnnotation(annot.down1000)
}

# -----------------------------------------------------------------------------
# Gene part annotations
# -----------------------------------------------------------------------------

tss_seqnames <- unique(seqnames(gene.parts$TSSes))
moabs.DMR <- moabs.DMR[seqnames(moabs.DMR) %in% tss_seqnames]
moabs.DMR.Ann <- annotateWithGeneParts(moabs.DMR, gene.parts, intersect.chr = FALSE)
moabs.DMR.TSS <- build_dmr_tss_table(moabs.DMR, moabs.DMR.Ann, mapping, expr.coef, log_fc_column)
moabs.DMR.TSS$logFC <- expr.coef[as.character(moabs.DMR.TSS$ID), log_fc_column]
genomation::plotTargetAnnotation(moabs.DMR.Ann)

moabs.DMR.genes <- annotateWithFeatureFlank(moabs.DMR, genes$Gene, genes$UTR,
                                           feature.name = "Gene", flank.name = "UTR",
                                           intersect.chr = TRUE)
genomation::plotTargetAnnotation(moabs.DMR.genes,
                                 col = c("green", "gray", "white"),
                                 main = "DMRs overlapped with genes")

# -----------------------------------------------------------------------------
# Long DMR inspection
# -----------------------------------------------------------------------------

moabs.longDMR <- subset(moabs.DMR.df, width > 1000)
if (nrow(moabs.longDMR) > 0 && exists("BSobj")) {
  try(showOneDMR(moabs.longDMR[min(8, nrow(moabs.longDMR)), ], BSobj), silent = TRUE)
}

# -----------------------------------------------------------------------------
# Biominer and DE gene analyses
# -----------------------------------------------------------------------------

biomine85 <- read.table(config$biominer_gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
biomine85.DMR.mean <- plot_DMR_and_LFC(moabs.DMR, moabs.DMR.TSS, biomine85$ID)

annot.DE <- read.csv(config$de_annotation_file, header = TRUE, row.names = 1, sep = "\t")
DE_ids <- intersect(moabs.DMR.TSS$ID, rownames(annot.DE))
DE.DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, DE_ids)
DE.DMR.mean <- aggregate_DMR(DE.DMR.TSS)
DE.DMR.mean <- plot_Mean_Meth_Diff_and_LFC(DE.DMR.mean)
DE.DMR.mean$Description <- annot.DE[as.character(DE.DMR.mean$ID), "Description"]
DE.DMR.mean.up <- subset(DE.DMR.mean, logFC > 0)
DE.DMR.mean.down <- subset(DE.DMR.mean, logFC < 0)
write.csv(DE.DMR.mean.up, file = "DE.DMR.mean.up.csv")
write.csv(DE.DMR.mean.down, file = "DE.DMR.mean.down.csv")

png(filename = file.path(config$output_dir, "DMGs_DMR_LFC.png"), width = 7800, height = 5400,
    res = config$plot_res)
plot_DMR_and_LFC(moabs.DMR, moabs.DMR.TSS, DMGs.CpG$ID)
dev.off()

DE.DMR.TSS.promoters <- subset(DE.DMR.TSS, dist.to.feature <= 1000 & dist.to.feature >= -1000)
if (nrow(DE.DMR.TSS.promoters) > 0) {
  plot_DMRs(DE.DMR.TSS.promoters)
}

if (exists("shared.up")) {
  shared_up_ids <- intersect(moabs.DMR.TSS$ID, rownames(shared.up))
  shared_up.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, shared_up_ids)
  shared_up.mean <- aggregate_DMR(shared_up.TSS)
  plot_Mean_Meth_Diff_and_LFC(shared_up.mean)
}

if (exists("carbonic_anhydrases")) {
  carbonic_DMR_ids <- intersect(moabs.DMR.TSS$ID, rownames(carbonic_anhydrases))
  carbonic_DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, carbonic_DMR_ids)
  carbonic_DMR.mean <- aggregate_DMR(carbonic_DMR.TSS)
  plot_Mean_Meth_Diff_and_LFC(carbonic_DMR.mean)
}

EHUX_ids <- moabs.DMR.TSS$ID
EHUX.DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, EHUX_ids)
EHUX.DMR.mean <- aggregate_DMR(EHUX.DMR.TSS)
plot_Mean_Meth_Diff_and_LFC(EHUX.DMR.mean)

# -----------------------------------------------------------------------------
# Locus-specific visualisation
# -----------------------------------------------------------------------------

if (exists("myDiff10")) {
  try(plot_DMR_Gviz(moabs.DMR, moabs.DMR.TSS, "233823", methDiff = myDiff10), silent = TRUE)
}
if (exists("myDiff.DSS")) {
  try(plot_DMR_Gviz(moabs.DMR, moabs.DMR.TSS, "233304", methDiff = myDiff.DSS), silent = TRUE)
  try(plot_DMR_Gviz(moabs.DMR, moabs.DMR.TSS, "193771", methDiff = myDiff.DSS), silent = TRUE)
  try(plotOneDMR(moabs.DMR, moabs.DMR.TSS, "193771", methDiff = myDiff.DSS), silent = TRUE)
}

if (exists("data.raw")) {
  id <- "193771"
  up_len <- 3000
  rna <- rownames(mapping[mapping$ID == id, ])
  TSSes <- gene.parts$TSSes
  exons <- gene.parts$exons
  upstream <- TSSes[TSSes$name == rna[1], ]
  if (length(upstream) > 0) {
    if (as.logical(strand(upstream) == '-')) {
      end(upstream) <- end(upstream) + up_len
    } else {
      start(upstream) <- start(upstream) - up_len
    }
    ex <- exons[exons$name == rna[1], ]
    gene_n_up <- c(upstream, ex)
    gene_n_up <- range(gene_n_up)
    strand(gene_n_up) <- '*'
    methhits <- getMethInGRange(data.raw, gene_n_up)
    matplot(methhits$start, methhits$meth.diff, type = "o", pch = 20, col = 'blue',
            xlab = paste("Genomic location", methhits$chr[1]),
            ylab = "Methylation difference")
    abline(h = c(-10, 0, 10), lty = 2)
  }
}

# -----------------------------------------------------------------------------
# Contingency analyses
# -----------------------------------------------------------------------------

non_DMGs_IDs <- setdiff(genes$Gene$ID, DMGs$ID)
non_DEGs_IDs <- setdiff(genes$Gene$ID, rownames(annot.DE))
DMGs_IDs <- DMGs$ID
DEGS_IDs <- rownames(annot.DE)
n_11 <- length(intersect(DMGs_IDs, DEGS_IDs))
n_12 <- length(intersect(DMGs_IDs, non_DEGs_IDs))
n_21 <- length(intersect(non_DMGs_IDs, DEGS_IDs))
n_22 <- length(intersect(non_DMGs_IDs, non_DEGs_IDs))
DMG_n_DEG <- matrix(c(n_11, n_12, n_21, n_22), nrow = 2, byrow = TRUE)
dimnames(DMG_n_DEG) <- list("DMG" = c("DMG_Y", "DMG_N"),
                           "DEG" = c("DEG_Y", "DEG_N"))
test.DMG_n_DEG <- chisq.test(DMG_n_DEG)

if (exists("DMPs")) {
  DMPs_IDs <- DMPs$ID
  non_DMPs_IDs <- setdiff(genes$Gene$ID, DMPs_IDs)
  n_11 <- length(intersect(DMPs_IDs, DEGS_IDs))
  n_12 <- length(intersect(DMPs_IDs, non_DEGs_IDs))
  n_21 <- length(intersect(non_DMPs_IDs, DEGS_IDs))
  n_22 <- length(intersect(non_DMPs_IDs, non_DEGs_IDs))
  DMP_n_DEG <- matrix(c(n_11, n_12, n_21, n_22), nrow = 2, byrow = TRUE)
  dimnames(DMP_n_DEG) <- list("DMP" = c("DMP_Y", "DMP_N"),
                              "DEG" = c("DEG_Y", "DEG_N"))
  test.DMP_n_DEG <- chisq.test(DMP_n_DEG)
}
