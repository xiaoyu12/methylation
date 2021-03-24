library(methylKit)
library(genomation)
library(stringr)
library(Biobase)

# Analyze DMRs related to gene features and expression profiles
# DMR: DMR dataframe. It has a field named "meth.diff"
# features: Genomic features such as promoters, cpgi etc.
# gene_exps: dataframe of differentially expressed genes with a field "logFC"
# selected_ids: consider a subset of DMRs for given selected_ids 
DMR_gene_analysis <- function(DMR, features, gene_exps, selected_ids = NULL, ignore.strand = TRUE) {
  DMR.gr <- makeGRangesFromDataFrame(DMR, keep.extra.columns = TRUE, ignore.strand = ignore.strand)
  DMR.Ann <- annotateWithGeneParts(DMR.gr, features, intersect.chr = TRUE)
  DMR.TSS <- getAssociationWithTSS(DMR.Ann)
  genomation::plotTargetAnnotation(DMR.Ann)
  DMR.TSS <- cbind(DMR[DMR.TSS$target.row, ], DMR.TSS)
  # mapping feature name to gene ID's
  DMR.TSS$ID <- mapping[DMR.TSS$feature.name, 2]
  DMR.TSS$logFC <- gene_exps[as.character(DMR.TSS$ID), ]$logFC
  if(!is.null(selected_ids)) {
    ids <- intersect(DMR.TSS$ID, selected_ids)
    DMR.TSS.sel <- DMR.TSS[which(DMR.TSS$ID %in% ids), ]
  }
  newList <- list("ann" = DMR.Ann, "tss" = DMR.TSS, "sel" = DMR.TSS.sel)
  return (newList)
}

methCHG_DSS_mincov10$DMR$meth.diff <- -100 * methCHG_DSS_mincov10$DMR$diff.Methy
dmr_CHG <- DMR_gene_analysis(methCHG_DSS_mincov10$DMR, gene.parts, expr.coef, selected_ids = rownames(annot.DE))
plot(dmr_CHG$sel$meth.diff, dmr_CHG$sel$logFC, type="p", cex = 1, pch = 20, col= "blue", xlab = "Methylation Difference", ylab = "Log Fold Change")
text(dmr_CHG$sel$meth.diff, dmr_CHG$sel$logFC + 0.25, labels = dmr_CHG$sel$ID, cex = 0.5)


# List of CpG coverage files
file.list <- list("EH1516C.merged_CpG_evidence.cov",
                  "EH1516D.merged_CpG_evidence.cov",
                  "EH217A.merged_CpG_evidence.cov",
                  "EH217B.merged_CpG_evidence.cov",
                  "EH217C.merged_CpG_evidence.cov")


# Read Ehux gene structures
gene.parts <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, unique.prom = FALSE)
# Promoters defined as up.flank = 1000, down.flank = 0
gene.parts_up1000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, 
                                            up.flank = 1000, down.flank = 0, unique.prom = FALSE)
# Promoters defined as up.flank = 0, down.flank = 1000
gene.parts_dn1000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, 
                                            up.flank = 0, down.flank = 1000, unique.prom = FALSE)

sample.ids = c("EH1516C", "EH1516D", "EH217A", "EH217B", "EH217C")
# Read methylation data
methobj <- methRead(file.list, sample.id = as.list(sample.ids),
                    assembly = "ehux", pipeline = "bismarkCoverage", treatment = c(0, 0, 1, 1, 1),
                    context = "CpG", mincov = 3, header = FALSE)

#Unite methylRawList to a single table
meth <- unite(methobj)
rm(methobj)

myDiff.DSS <- calculateDiffMethDSS(meth, mc.cores = 4)

getCorrelation(meth,plot=FALSE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, scale = FALSE)

meth10 <- unite(filterByCoverage(methobj, lo.count = 10))
getCorrelation(meth10, plot = FALSE)
clusterSamples(meth10, dist="correlation", plot=TRUE)
PCASamples(meth10, scale = FALSE)
myDiff10 <- calculateDiffMeth(meth10)
#myDiff15 <- calculateDiffMeth(unite(filterByCoverage(methobj, lo.count = 15)))

meth.CpG.cov10 <- unite(methCpG, destrand=TRUE)
PCASamples(meth.CpG.cov10, scale = FALSE)
getCorrelation(meth.CpG.cov10, plot = FALSE)

meth.CHG.cov10 <- unite(filterByCoverage(methCHG, lo.count = 10))
getCorrelation(meth.CHG.cov10, plot = FALSE)
clusterSamples(meth.CHG.cov10, dist="correlation", plot=TRUE)
PCASamples(meth.CHG.cov10, scale = FALSE)

# Read mapping from rna* to gene ID
mapping <- read.table("../v2/genbank_mapping.txt", row.names = 1)
colnames(mapping) <- c("gene", "ID")

m <- getFeatureMethyl(methobj, gene.parts$promoters, sample.ids)
m_up1000 <- getFeatureMethyl(methobj, gene.parts_up1000$promoters, sample.ids)
m_dn1000 <- getFeatureMethyl(methobj, gene.parts_dn1000$promoters, sample.ids)

m.EH1516.beta <- m$beta[, 1:2]
m.EH217.beta <- m$beta[, 3:5]

m.DMR.beta <- m$beta[which(abs(rowMedians(m.EH1516.beta) - rowMedians(m.EH217.beta)) >= 0.2), ]

# read Ehux gene annotations
annotations <- read.csv("../v2/ehux_JGI_best_b2g.txt", header=TRUE, row.names=1, sep="\t")
annot.DE <- read.csv("../../gene_expression/Ehux_1516_v_217/Jan2019/output/DE_E1516_vs_E217.annot.txt", 
                     header = TRUE, row.names = 1, sep="\t")
#annot.DE <- read.csv("../../gene_expression/Ehux_1516_v_217/Jan2019/output_fdr0.05//DE_E1516_vs_E217.annot.txt", 
#                     header = TRUE, row.names = 1, sep="\t")
annot.DE.down <- subset(annot.DE, log2FC < 0)
annot.DE.up <- subset(annot.DE, log2FC > 0)

# Read list of carbonic anhydrases
carbonic_anhydrases <- read.csv("biom_carbonic_anhydrases.txt", header = FALSE, row.names = 1, sep = "\t")
# Plot Average expression and expression LFC of carbonic anhydrases genes
carbonic_expr <- expr.coef[rownames(carbonic_anhydrases), ]
plot(carbonic_expr$baseMean, carbonic_expr$LFC, pch=20, xlab = "Avg Expression", ylab="logFC", col="blue")
carbonic_de = subset(carbonic_expr, abs(LFC) > 1)
text(carbonic_de$baseMean+runif(1, min=-2.5, max=2.5), carbonic_de$LFC+0.5, labels = rownames(carbonic_de), cex=0.5)
points(carbonic_de$baseMean, carbonic_de$LFC, pch=20, , col="red")
abline(h = 0)
#abline(v = 0)
abline(h = 1, lty = 2, col="gray")
abline(h = -1, lty=2, col="gray")

annot.DE.DMR <- annot.DE[intersect(rownames(annot.DE), rownames(m.DMR.beta)), ]

m10 <- getFeatureMethyl(filterByCoverage(methobj, lo.count = 10), gene.parts$promoters, sample.ids)
m10_up1000 <- getFeatureMethyl(filterByCoverage(methobj, lo.count = 10), gene.parts_up1000$promoters, sample.ids)

m_exon <- getFeatureMethyl(methobj, gene.parts$exons, sample.ids)
m_intron <- getFeatureMethyl(methobj, gene.parts$introns, sample.ids)


# Diff mehtylated in promoter regions
promobj <- regionCounts(methobj, gene.parts$promoters)
meth.prom <- unite(filterByCoverage(promobj, lo.count = 100))
getCorrelation(meth.prom, plot = FALSE)
meth.prom.DSS <- calculateDiffMeth(meth.prom)
meth.prom.DSS.20p <- getMethylDiff(meth.prom.DSS, difference = 10, qvalue = 0.01)
promoters <- as.data.frame(gene.parts$promoters)
promoters$ID = mapping[promoters$name, "ID"]
promoters = promoters[!duplicated(promoters$ID), ]    # remove duplicates due to alternative splicing
colnames(promoters)[1] <- "chr"
meth.prom.DSS.data <- merge(getData(meth.prom.DSS), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
meth.prom.DSS.20p <- merge(getData(meth.prom.DSS.20p), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
plot(meth.prom.DSS.20p$meth.diff, expr.coef[as.character(meth.prom.DSS.20p$ID), ]$LFC, pch=20, xlab="Promoter Methyl Diff", ylab="LFC", col="blue")
text(meth.prom.DSS.20p$meth.diff, expr.coef[as.character(meth.prom.DSS.20p$ID), ]$LFC+0.5, labels = as.character(meth.prom.DSS.20p$ID), cex=0.5)
abline(h=0)
abline(v=0)
# make tiles of size 100
tiles100 = tileMethylCounts(filterByCoverage(methobj, lo.count = 10), win.size = 100, step.size = 100)
