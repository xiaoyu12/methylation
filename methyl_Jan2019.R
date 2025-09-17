library(methylKit)
library(genomation)
library(stringr)
library(Biobase)
library(pheatmap)

source("read_bismark_methyl.R")
source("meth_analysis_funcs.R")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
config <- list(
  data_dir = file.path("data", "methylKit"),
  output_dir = ".",
  genome_bed = file.path("..", "v2", "Ehux_genbank.bed"),
  mapping_file = "genbank_mapping.txt",
  expr_coef_file = "Ehux_217_vs_1516.LFC.csv",
  de_annotation_file = file.path("..", "..", "gene_expression", "Ehux_1516_v_217",
                                 "Jan2019", "output", "DE_E1516_vs_E217.annot.txt"),
  carbonic_list = "biom_carbonic_anhydrases.txt",
  diff_threshold = list(difference = 25, qvalue = 0.01),
  promoter_diff_threshold = list(difference = 20, qvalue = 0.01),
  min_cov = list(default = 3, high = 10, promoter = 100),
  dss_cores = 4,
  chg_cores = 10
)

out_path <- function(filename) file.path(config$output_dir, filename)

# Helper to generate consistent QC outputs for each methylBase object
generate_qc_plots <- function(meth, prefix, meta, adj_lim = c(0.4, 0.1)) {
  getCorrelation(meth, plot = FALSE)
  plot_correlation_heatmap(meth, out_path(paste0(prefix, "_pheatmap.png")), meta)
  clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
  plot_pca(meth, out_path(paste0(prefix, "_PCA.png")), adj_lim = adj_lim)
}

manifest <- default_methyl_manifest()
meta <- build_sample_metadata(manifest)
sample_ids <- manifest$sample_id

inputs <- load_inputs(list(
  mapping_file = config$mapping_file,
  genome_bed = config$genome_bed,
  expr_file = config$expr_coef_file
))

mapping <- inputs$mapping
expr.coef <- inputs$expr
if (is.null(expr.coef)) {
  stop("Expression coefficients file not found: ", config$expr_coef_file)
}

gene.parts <- inputs$gene_parts
gene.parts_up1000 <- inputs$gene_parts_up1000
gene.parts_dn1000 <- inputs$gene_parts_dn1000
gene.parts_up3000 <- readTranscriptFeatures(
  config$genome_bed,
  remove.unusual = FALSE,
  unique.prom = FALSE,
  up.flank = 3000,
  down.flank = 0
)

options(methylation.genome_bed = config$genome_bed)

annot.DE <- read.csv(config$de_annotation_file, header = TRUE, row.names = 1, sep = "\t")
annot.DE.down <- subset(annot.DE, log2FC < 0)
annot.DE.up <- subset(annot.DE, log2FC > 0)

carbonic_anhydrases <- read.csv(config$carbonic_list, header = FALSE, row.names = 1, sep = "\t")

# -----------------------------------------------------------------------------
# Load methylation data
# -----------------------------------------------------------------------------
methobj <- read_Bismark_coverage(
  base_dir = config$data_dir,
  manifest = manifest,
  assembly = "ehux",
  mincov = config$min_cov$default
)

meth <- methylKit::unite(methobj)
meth10 <- methylKit::unite(
  filterByCoverage(methobj, lo.count = config$min_cov$high)
)

# -----------------------------------------------------------------------------
# QC: correlation, clustering, PCA
# -----------------------------------------------------------------------------
meth_cor <- compute_methyl_correlation(meth)
generate_qc_plots(meth, "meth", meta, adj_lim = c(0.4, 0.1))
generate_qc_plots(meth10, "meth_cov10", meta, adj_lim = c(0.3, 0.1))

# -----------------------------------------------------------------------------
# Differential methylation analyses
# -----------------------------------------------------------------------------
myDiff.DSS <- calculateDiffMethDSS(meth, mc.cores = config$dss_cores)
sigDiff.DSS <- getMethylDiff(
  myDiff.DSS,
  difference = config$diff_threshold$difference,
  qvalue = config$diff_threshold$qvalue
)

meth10_diff <- run_methylkit_dm(
  meth10,
  diff_percent = config$diff_threshold$difference,
  qvalue = config$diff_threshold$qvalue
)
myDiff10 <- meth10_diff$all
sigDiff10 <- meth10_diff$sig

sigDiff10.o.gp <- overlap_dmr_gene_features(sigDiff10, gene.parts, expr.coef, mapping)
sigDiff10.o.up3000 <- overlap_dmr_gene_features(sigDiff10, gene.parts_up3000, expr.coef, mapping)

up3000 <- getPromoterRegions(up = 3000, down = 0)
up1000 <- getPromoterRegions(up = 1000, down = 0)

genes.DE_no_vrlp <- read_selected_gene_features(
  "DE_no_overlap.tsv",
  genome_bed = config$genome_bed,
  up_flank = 1000,
  down_flank = 1000
)

# -----------------------------------------------------------------------------
# CHG context analysis
# -----------------------------------------------------------------------------
bismark <- read_Bismark_CpG_CHG(
  base_dir = config$data_dir,
  manifest = manifest,
  assembly = "ehux",
  mincov = config$min_cov$default
)

meth.CHG.cov10 <- methylKit::unite(
  filterByCoverage(bismark$methCHG, lo.count = config$min_cov$high)
)

generate_qc_plots(meth.CHG.cov10, "CHG", meta, adj_lim = c(0.4, 0.1))

CHG.cov10.Diff <- calculateDiffMeth(meth.CHG.cov10, mc.cores = config$chg_cores)
sig.CHG.cov10.Diff <- getMethylDiff(
  CHG.cov10.Diff,
  difference = config$diff_threshold$difference,
  qvalue = config$diff_threshold$qvalue
)

# -----------------------------------------------------------------------------
# Feature-level methylation summaries
# -----------------------------------------------------------------------------
m <- getFeatureMethyl(methobj, gene.parts$promoters, sample_ids, mapping = mapping)
m_up1000 <- getFeatureMethyl(methobj, gene.parts_up1000$promoters, sample_ids, mapping = mapping)
m_dn1000 <- getFeatureMethyl(methobj, gene.parts_dn1000$promoters, sample_ids, mapping = mapping)

m.EH1516.beta <- m$beta[, 1:2, drop = FALSE]
m.EH217.beta <- m$beta[, 3:5, drop = FALSE]

m.DMR.beta <- m$beta[
  abs(rowMedians(m.EH1516.beta) - rowMedians(m.EH217.beta)) >= 0.2,
  ,
  drop = FALSE
]

annot.DE.DMR <- annot.DE[intersect(rownames(annot.DE), rownames(m.DMR.beta)), ]

m10 <- getFeatureMethyl(
  filterByCoverage(methobj, lo.count = config$min_cov$high),
  gene.parts$promoters,
  sample_ids,
  mapping = mapping
)
m10_up1000 <- getFeatureMethyl(
  filterByCoverage(methobj, lo.count = config$min_cov$high),
  gene.parts_up1000$promoters,
  sample_ids,
  mapping = mapping
)

m_exon <- getFeatureMethyl(methobj, gene.parts$exons, sample_ids, mapping = mapping)
m_intron <- getFeatureMethyl(methobj, gene.parts$introns, sample_ids, mapping = mapping)

# -----------------------------------------------------------------------------
# Carbonic anhydrase genes: expression vs fold-change
# -----------------------------------------------------------------------------
carbonic_expr <- expr.coef[rownames(carbonic_anhydrases), ]
png(filename = out_path("carbonic_lfc.png"), width = 1600, height = 1600, res = 200)
plot(
  carbonic_expr$baseMean,
  carbonic_expr$LFC,
  pch = 20,
  xlab = "Avg Expression",
  ylab = "logFC",
  col = "blue"
)
carbonic_de <- subset(carbonic_expr, abs(LFC) > 1)
text(
  carbonic_de$baseMean + runif(1, min = -2.5, max = 2.5),
  carbonic_de$LFC + 0.5,
  labels = rownames(carbonic_de),
  cex = 0.5
)
points(carbonic_de$baseMean, carbonic_de$LFC, pch = 20, col = "red")
abline(h = 0)
abline(h = c(1, -1), lty = 2, col = "gray")
dev.off()

# -----------------------------------------------------------------------------
# Promoter-level differential methylation
# -----------------------------------------------------------------------------
promobj <- regionCounts(methobj, gene.parts$promoters)
meth.prom <- methylKit::unite(
  filterByCoverage(promobj, lo.count = config$min_cov$promoter)
)
getCorrelation(meth.prom, plot = FALSE)
meth.prom.DSS <- calculateDiffMeth(meth.prom)
meth.prom.DSS.20p <- getMethylDiff(
  meth.prom.DSS,
  difference = config$promoter_diff_threshold$difference,
  qvalue = config$promoter_diff_threshold$qvalue
)

promoters <- as.data.frame(gene.parts$promoters)
promoters$ID <- mapping[promoters$name, "ID"]
promoters <- promoters[!duplicated(promoters$ID), ]
colnames(promoters)[1] <- "chr"

meth.prom.DSS.data <- merge(
  getData(meth.prom.DSS),
  promoters,
  by.x = c("chr", "start", "end", "strand"),
  by.y = c("chr", "start", "end", "strand")
)

meth.prom.DSS.20p <- merge(
  getData(meth.prom.DSS.20p),
  promoters,
  by.x = c("chr", "start", "end", "strand"),
  by.y = c("chr", "start", "end", "strand")
)

png(filename = out_path("meth_promoter_diff_vs_lfc.png"), width = 1600, height = 1600, res = 200)
plot(
  meth.prom.DSS.20p$meth.diff,
  expr.coef[as.character(meth.prom.DSS.20p$ID), ]$LFC,
  pch = 20,
  xlab = "Promoter Methyl Diff",
  ylab = "LFC",
  col = "blue"
)
text(
  meth.prom.DSS.20p$meth.diff,
  expr.coef[as.character(meth.prom.DSS.20p$ID), ]$LFC + 0.5,
  labels = as.character(meth.prom.DSS.20p$ID),
  cex = 0.5
)
abline(h = 0)
abline(v = 0)
dev.off()
