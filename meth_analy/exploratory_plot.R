# Exploratory methylation analysis helper script
# --------------------------------------------------------------
# Expected objects in the global environment before sourcing:
#   - promobj: BSseq object with promoter methylation calls
#   - promoters: GenomicRanges describing promoter regions
#   - sample.ids: character vector naming the samples in promobj
#   - conditions: factor (or coercible vector) defining sample groups
#   - cpgi: list containing CpG island annotations (at least $CpGi)
#   - gene.parts: output from annotateWithGeneParts()
#   - mapping: lookup table with at least a column named "ID"
#   - methobj: BSseq object with genome-wide methylation calls
# --------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
})

source("meth_analysis_funcs.R")

# ---- Validation ------------------------------------------------------------
required_objects <- c(
  "promobj", "promoters", "sample.ids", "conditions",
  "cpgi", "gene.parts", "mapping", "methobj"
)
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "The following objects must exist before running this script: ",
    paste(missing_objects, collapse = ", ")
  )
}

if (!is.factor(conditions)) {
  conditions <- factor(conditions)
}

# Ensure sample naming is consistent across objects
promoter_meth <- methyl_to_data_frame(promobj, promoters)
stopifnot(
  ncol(promoter_meth) == length(sample.ids),
  identical(colnames(promoter_meth), colnames(promoter_meth))
)
colnames(promoter_meth) <- sample.ids

if (length(conditions) != length(sample.ids)) {
  stop("Length of `conditions` (", length(conditions),
       ") must match number of samples (", length(sample.ids), ")")
}

# ---- Visualization ---------------------------------------------------------
plot_palette <- grDevices::rainbow(length(sample.ids))

plot_promoter_density <- function(meth_matrix, palette) {
  densities <- apply(meth_matrix, 2, density)
  max_y <- max(vapply(densities, function(d) max(d$y), numeric(1)))

  plot(
    densities[[1]], col = palette[1], ylim = c(0, max_y),
    main = "Promoter methylation density", xlab = "Methylation M-values"
  )
  if (length(densities) > 1) {
    for (i in 2:length(densities)) {
      lines(densities[[i]], col = palette[i])
    }
  }
  legend(
    "topright", inset = 0.05, legend = colnames(meth_matrix),
    fill = palette, cex = 0.7
  )
}

plot_promoter_density(promoter_meth, plot_palette)

plotMDS(
  x = as.matrix(promoter_meth), top = 1000,
  col = plot_palette, dim = c(1, 2)
)

# ---- Differential analysis (promoters) -------------------------------------
design <- model.matrix(~ conditions)
rownames(design) <- sample.ids

fit <- lmFit(as.matrix(promoter_meth), design)
fit <- eBayes(fit)

limma_summary <- summary(decideTests(fit))
top_promoter_table <- topTable(fit, number = Inf, coef = 2)
promoter_dm <- subset(top_promoter_table, adj.P.Val < 0.05)

message(
  "Promoter differential methylation: ",
  nrow(promoter_dm), " regions (adj.P.Val < 0.05)"
)

# ---- CpG island analysis ---------------------------------------------------
if (!"CpGi" %in% names(cpgi)) {
  stop("cpgi list must contain a component named 'CpGi'.")
}

cpg_islands <- cpgi$CpGi
if (length(cpg_islands$name) == 0L && length(gene.parts$TSSes$name) > 0L) {
  cpg_islands$name <- gene.parts$TSSes$name
}

if (!is.null(mapping) && is.null(cpg_islands$ID)) {
  common_names <- intersect(cpg_islands$name, rownames(mapping))
  if (length(common_names) == 0) {
    warning("No overlapping names between CpG islands and mapping table.")
  } else {
    cpg_islands$ID <- mapping[cpg_islands$name, "ID"]
  }
}

cpg_island_meth <- getFeatureMethyl(methobj, cpg_islands, sample.ids)
cpg_fit <- lmFit(as.matrix(cpg_island_meth$m), design)
cpg_fit <- eBayes(cpg_fit)

cpg_summary <- summary(decideTests(cpg_fit))
cpgi_table <- topTable(cpg_fit, number = Inf, coef = 2)
cpg_island_dm <- subset(cpgi_table, adj.P.Val < 0.05)

message(
  "CpG island differential methylation: ",
  nrow(cpg_island_dm), " regions (adj.P.Val < 0.05)"
)

# ---- Output objects --------------------------------------------------------
analysis_outputs <- list(
  palette = plot_palette,
  promoter_density = promoter_meth,
  promoter_fit = fit,
  promoter_summary = limma_summary,
  promoter_dm = promoter_dm,
  cpg_fit = cpg_fit,
  cpg_summary = cpg_summary,
  cpg_dm = cpg_island_dm
)

message("Analysis complete. Outputs available in `analysis_outputs`.")
