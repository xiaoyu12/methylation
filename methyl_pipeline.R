suppressPackageStartupMessages({
  library(methylKit)
  library(genomation)
  library(pheatmap)
  library(BiocParallel)
})

source("read_bismark_methyl.R")
source("meth_analysis_funcs.R")

config <- list(
  data_dir = "./data/methylKit/",
  genome_bed = "../v2/Ehux_genbank.bed",
  mapping_file = "./genbank_mapping.txt",
  expr_file = "../../gene_expression/Ehux_1516_v_217/Jan2019/output/DE_E1516_vs_E217.annot.txt",
  out_dir = ".",
  cache_dir = "./cache",
  min_cov = 10,
  diff_percent = 25,
  qvalue = 0.01,
  promoter_up = 1000,
  promoter_down = 0,
  mc_cores = max(1, parallel::detectCores() - 1),
  load_chg = TRUE
)

register(MulticoreParam(workers = config$mc_cores))

# Load inputs (mapping, gene features, expression)
inputs <- load_inputs(config)

# Build methylation objects (cached)
objs <- build_methyl_objects(config)

# Metadata for plotting
strains <- factor(c(rep("EH1516",2), rep("EH217", 3)))
meta <- data.frame(strains)
rownames(meta) <- slot(objs$meth, "sample.ids")

# QC plots
plot_correlation_heatmap(objs$meth, file.path(config$out_dir, "meth_pheatmap.png"), meta)
plot_pca(objs$meth, file.path(config$out_dir, "meth_PCA.png"))

# Differential methylation (CpG)
dm <- run_methylkit_dm(objs$meth10, diff_percent = config$diff_percent, qvalue = config$qvalue, mc.cores = config$mc_cores)

# Overlap with gene features
sigDiff10 <- dm$sig
ann_default <- overlap_DMR_GeneFeatures2(sigDiff10, inputs$gene_parts, inputs$expr, inputs$mapping)
ann_up3000 <- overlap_DMR_GeneFeatures2(sigDiff10, readTranscriptFeatures(config$genome_bed, remove.unusual = FALSE, up.flank = 3000, down.flank = 0, unique.prom = FALSE), inputs$expr, inputs$mapping)

# Save results
ensure_dir(config$out_dir)
saveRDS(list(dm = dm, ann_default = ann_default, ann_up3000 = ann_up3000), file.path(config$out_dir, "results.rds"))

message("Pipeline completed. Results saved to results.rds")
