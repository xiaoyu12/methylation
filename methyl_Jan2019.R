library(methylKit)
library(genomation)
library(stringr)
library(Biobase)
#source("read_bismark_methyl.R")

methobj <- read_Bismark_coverage()
#Unite methylRawList to a single table
meth <- methylKit::unite(methobj)
rm(methobj)

# Add a column of gene IDs to the Granges or dataframe
addMapID <- function(ranges) {
  mapping <- read.table("../v2/genbank_mapping.txt", row.names = 1)
  colnames(mapping) <- c("gene", "ID")
  ranges$ID = mapping[ranges$name, "ID"]
  ranges <- as.data.frame(ranges)
  ranges = ranges[!duplicated(ranges$ID), ]    # remove duplicates due to alternative splicing
  # convert back to Granges
  ranges <- as(ranges, "GRanges")
  return(ranges)
}

# Get Granges object representing regions upstream and downstream from TSSes with given len
getPromoterRegions <- function(up = 1000, down = 1000) {
  gene_structs <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE,
                                         up.flank = up, down.flank = down, unique.prom = FALSE)
  promoters <- addMapID(gene_structs$promoter)
  return(promoters)
}

up3000 <- getPromoterRegions(up = 3000, down = 0)
# get methylation data in the upstream 3000 regions
meth.up3000 <- selectByOverlap(meth, up3000)

up1000 <- getPromoterRegions(up = 1000, down = 0)
meth.up1000 <- selectByOverlap(meth, up1000)

down1000 <- getPromoterRegions(up = 0, down = 1000)

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



# Read Ehux gene structures
gene.parts <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, unique.prom = FALSE)
# Promoters defined as up.flank = 1000, down.flank = 0
gene.parts_up1000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, 
                                            up.flank = 1000, down.flank = 0, unique.prom = FALSE)
# Promoters defined as up.flank = 0, down.flank = 1000
gene.parts_dn1000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, 
                                            up.flank = 0, down.flank = 1000, unique.prom = FALSE)

readSelectedGenes <- function(gene_list_file, up_flank = 1000, down_flank = 1000) {
  # read gene_list into a data frame
  tmp_file = "./Ehux_tmp.bed"
  gene_list = read.table(gene_list_file, sep = "\t", header = TRUE, as.is = TRUE)
  gene_bed = read.table("../v2/Ehux_genbank.bed", sep="\t", header = FALSE, as.is = TRUE)
  gene_bed = gene_bed[gene_bed[, 4] %in% gene_list$rna, ]
  write.table(gene_bed, tmp_file, sep = "\t", row.names = FALSE, col.names =FALSE, quote = FALSE)
  features <- readTranscriptFeatures(tmp_file, remove.unusual = FALSE, unique.prom = FALSE,
                         up.flank = up_flank, down.flank = down_flank)
  file.remove(tmp_file)
  return(features)
}

genes.DE_no_vrlp <- readSelectedGenes("DE_no_overlap.tsv")


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

meth.CpG.cov10 <- methylKit::unite(methCpG, destrand=TRUE)
PCASamples(meth.CpG.cov10, scale = FALSE)
getCorrelation(meth.CpG.cov10, plot = FALSE)

meth.CHG.cov10 <- methylKit::unite(filterByCoverage(methCHG, lo.count = 10))
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
