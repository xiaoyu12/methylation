library(methylKit)
library(genomation)
library(stringr)
library(Biobase)
library(pheatmap)
source("read_bismark_methyl.R")
source("meth_analysis_funcs.R")

# Read Bismark coverage file for CpG methylation
methobj <- read_Bismark_coverage("./data/methylKit/")
# Unite methylRawList to a single table
# meth has a minimum coverage of 3 for each sample
meth <- methylKit::unite(methobj)
rm(methobj)

strains <- factor(c(rep("EH1516",2), rep("EH217", 3)))
#sample.ids = c("EH1516B", "EH1516C", "EH217A", "EH217B", "EH217C")
meta <- data.frame(strains, row.names = slot(meth, "sample.ids"))

# calculate the correlation of sample meth profiles
# It should get the same results as the getCorrelation() function in methylKit
# but it returns a correlation matrix that can be used for pheatmap plot
getMethCorrelation <- function(meth) {
  x <- getData(meth)
  x.cov <- x[, slot(meth, "coverage.index")]
  x.C <- x[, slot(meth, "numCs.index")]
  x <- x.C / x.cov
  colnames(x) <- slot(meth, "sample.ids")
  meth.cor <- cor(x)
  return (meth.cor)
}

# Cluster samples based on correlations between CpG methylations
getCorrelation(meth,plot=FALSE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

# Plot correlation
meth.cor <- getMethCorrelation(meth)
png(filename = "meth_pheatmap.png", width=2400, height=1600, res=600)
pheatmap(meth.cor, annotation = meta)
dev.off()

# Plot PCA analysis
png(filename = "meth_PCA.png", width=4800, height=3200, res=600)
methylKit::PCASamples(meth, adj.lim=c(0.4, 0.1), scale = FALSE)
dev.off()

# Plot the methylation profile on scaffold_1 (chr1)
# use a 10-kb window
# make tiles of size 100
#tiles100 = tileMethylCounts(filterByCoverage(methobj, lo.count = 10), win.size = 100, step.size = 100)
tiles10K = tileMethylCounts(methobj, win.size = 10000, step.size = 10000)
meth.chr1 <- subset(getData(tiles10K), chr == "chr1")
#sum.chr1 <- matrix(rep(0, times=len*5), ncol=5)
#for (i in 1:len) {
#  meth.win <- subset(meth.chr1, start < i*10000 & start >= (i-1)*10000)
#}
ggplot(meth.chr1) + geom_smooth(mapping=aes(x = , y = ))

# calculate Differential Methylation with DSS
myDiff.DSS <- calculateDiffMethDSS(meth, mc.cores = 4)
sigDiff.DSS <- getMethylDiff(myDiff.DSS, difference = 25, qvalue = 0.01)

# meth10 has a min coverage of 10 for each sample
meth10 <- unite(filterByCoverage(methobj, lo.count = 10))
getCorrelation(meth10, plot = FALSE)
clusterSamples(meth10, dist="correlation", plot=TRUE)
methylKit::PCASamples(meth10, adj.lim=c(0.3, 0.1), scale = FALSE)
myDiff10 <- calculateDiffMeth(meth10)
myDiff10.DSS <- calculateDiffMethDSS(meth10, mc.cores = 24)
sigDiff10 <- getMethylDiff(myDiff10, difference = 25, qvalue = 0.01)
sigDiff10.DSS <- getMethylDiff(myDiff10.DSS, difference = 25, qvalue = 0.01)

# Read mapping from rna* to gene ID
mapping <- read.table("../v2/genbank_mapping.txt", row.names = 1)
colnames(mapping) <- c("gene", "ID")

# Read Ehux gene structures
# By default, promoters are defined as 1000 bp up and down stream from the TSSes 
gene.parts <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, unique.prom = FALSE)

# Get overlap of sigDiff10 and default gene parts. (up=1000, down=1000)
sigDiff10.o.gp <- overlap_DMR_GeneFeatures(sigDiff10, gene.parts, expr.coef)

# Read Ehux gene structures, using up 3000 and down 0 for promoter regions
gene.parts.up3000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE,
                       up.flank = 3000, down.flank = 0, unique.prom = FALSE)

# Get overlaps of sigDiff10 and gene.parts with up 3000 promoters
sigDiff10.o.up3000 <- overlap_DMR_GeneFeatures(sigDiff10, gene.parts.up3000, expr.coef)

up3000 <- getPromoterRegions(up = 3000, down = 0)
# get methylation data in the upstream 3000 regions
#meth.up3000 <- selectByOverlap(meth, up3000)

up1000 <- getPromoterRegions(up = 1000, down = 0)
#meth.up1000 <- selectByOverlap(meth, up1000)

#down1000 <- getPromoterRegions(up = 0, down = 1000)

#Read the features of only genes in a list file
readSelectedGeneFeatures <- function(gene_list_file, up_flank = 1000, down_flank = 1000) {
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

genes.DE_no_vrlp <- readSelectedGeneFeatures("DE_no_overlap.tsv")


#meth.CpG.cov10 <- methylKit::unite(methCpG)
#PCASamples(meth.CpG.cov10, scale = FALSE)
#getCorrelation(meth.CpG.cov10, plot = FALSE)

#*************************************************************************#
#* Analyze CHG context
#*************************************************************************#
# Read Bismark CpG and CHG methylation data
bismark <- read_Bismark_CpG_CHG(base_dir = "~/ESSD/data/methylKit/")

meth.CHG.cov10 <- methylKit::unite(filterByCoverage(bismark$methCHG, lo.count = 10))
CHG.cor <- getMethCorrelation(meth.CHG.cov10)
# Draw pheatmap
png(filename = "CHG_pheatmap.png", width=2400, height=1600, res=600)
pheatmap(CHG.cor, annotation = meta)
dev.off()
clusterSamples(meth.CHG.cov10, dist="correlation", plot=TRUE)
png(filename = "CHG_PCA.png", width=4800, height=3200, res=600)
PCASamples(meth.CHG.cov10, adj.lim=c(0.4, 0.1), scale = FALSE)
dev.off()

# calculate DMC in CHG context
CHG.cov10.Diff <- calculateDiffMeth(meth.CHG.cov10, mc.cores = 10)
sig.CHG.cov10.Diff <- getMethylDiff(CHG.cov10.Diff, difference = 25, qvalue = 0.01)

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
meth.prom.DSS.20p <- getMethylDiff(meth.prom.DSS, difference = 20, qvalue = 0.01)
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

