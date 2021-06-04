library(methylKit)
library(genomation)
library(stringr)
source("meth_analysis_funcs.R")

# read MOABS DMR from a file
readMOABS_DMR_bedfile <- function(filename) {
  library(rtracklayer)
  dmr <- import(filename, format = "bed")
  # calculate methyl difference
  name <- dmr$name
  name <- substring(name, 5)
  vals <- str_split_fixed(name, "vs", 2)
  # meth.diff = EH217 - EH1516
  meth.diff <- as.numeric(vals[, 1]) - as.numeric(vals[, 2])
  values(dmr) <- cbind(values(dmr), DataFrame(meth.diff))
  return(dmr)
}

# Read MOBS DMR file
moabs.DMR <- readMOABS_DMR_bedfile("dmr_M2_EH1516.G.bed_vs_EH217.G.bed.bed")
# Find the overlaps of DMRs and genes upstream regions
annot.up1000 <- annotateWithFeature(moabs.DMR, up1000)
genomation::plotTargetAnnotation(annot.up1000)
annot.down1000 <- annotateWithFeature(moabs.DMR, down1000)
genomation::plotTargetAnnotation(annot.down1000)

moabs.DMR.Ann <- annotateWithGeneParts(moabs.DMR, gene.parts, intersect.chr = TRUE)
moabs.DMR.TSS <- getAssociationWithTSS(moabs.DMR.Ann)
moabs.DMR.TSS$ID <- mapping[moabs.DMR.TSS$feature.name, 2]
moabs.DMR.TSS$logFC <- expr.coef[as.character(moabs.DMR.TSS$ID), ]$logFC
genomation::plotTargetAnnotation(moabs.DMR.Ann)
# Read cpg islands info
# The returned cpgi seem to be regions upstream from TSS
genes = readFeatureFlank("../v2/Ehux_genbank.bed", feature.flank.name=c("Gene","UTR"))
moabs.DMR.cpgi <- annotateWithFeatureFlank(as(moabs.DMR, "GRanges"), cpgi$CpGi,cpgi$shores,feature.name="CpGi",flank.name="shores", intersect.chr = TRUE)
genomation::plotTargetAnnotation(moabs.DMR.cpgi,col=c("green","gray","white"),main="differential methylation annotation")
# convert DMR to a dataframe
moabs.DMR.df <- as.data.frame(moabs.DMR)

cpg_i <- cpgi$CpGi
cpg_i$name = gene.parts$TSSes$name  
cpg_i$ID <- mapping[cpg_i$name, "ID"]
moabs.DMR.in.cgpi <- moabs.DMR.TSS[which(slot(moabs.DMR.cpgi, "members")[, 1] > 0), ] # not correct because moabs.DMR.cpgi has longer length than moabs.DMR.TSS
moabs.DMR.in.cgpi$meth.diff <- as.data.frame(moabs.DMR)[moabs.DMR.in.cgpi$target.row, "meth.diff"]
plot(moabs.DMR.in.cgpi$meth.diff, moabs.DMR.in.cgpi$logFC, cex=1, pch = 20, 
     col= "blue", xlab = "Methylation Difference", ylab = "Log Fold Change")
# get methylation Diff data for a set of ids
getMethDiff_for_Ids <- function(methDiff, diffTSS, ids) {
  diffTSS.sel <- diffTSS[diffTSS$ID %in% ids, ]
  row <- diffTSS[diffTSS$ID %in% ids, ]$target.row
  methDiff.sel <- data.frame(methDiff[row, ])
  methDiff.sel <- cbind(methDiff.sel, diffTSS.sel)
  return(methDiff.sel)
}

# get Log fold change in gene expression values for a list of transcript ids in form of "rna***"
getDiffExprCoef_for_Ids <- function(ids) {
  gids <- as.character(ids)
  data <- expr.coef[gids, ]$logFC
  return (data)
}

# plot MOABS methyl diff and LFC for selected genes
plotMethyl_and_LFC <- function(methDiffAll, diffTSS, ids) {
  methDiff <- getMethDiff_for_Ids(methDiffAll, diffTSS, ids)
  mids <- as.vector(methDiff$ID)
  lfc <- getDiffExprCoef_for_Ids(mids)
  methDiff <- cbind(methDiff, lfc)
  
  plot(methDiff$meth.diff, methDiff$lfc, type="p", cex = 1, pch = 20, col= "blue", xlab = "Methylation Difference", ylab = "Log Fold Change")
  text(methDiff$meth.diff, methDiff$lfc + 0.25, labels = methDiff$ID, cex = 0.5)
  return(methDiff)
}

# Compute aggreated distance and meth.diff to TSS features
aggregate_DMR <- function(dmrTSS) {
  dmr_genes <- aggregate(dist.to.feature~ID, data=dmrTSS, paste, collapse=";")
  dmr_genes <- merge(dmr_genes, aggregate(dist.to.feature~ID, data=dmrTSS, length), by="ID")
  dmr_genes <- merge(dmr_genes, aggregate(meth.diff~ID, data=dmrTSS, paste, collapse=";"), by="ID")
  dmr_genes <- merge(dmr_genes, aggregate(meth.diff~ID, data=dmrTSS, mean, collapse=";"), by="ID")
  colnames(dmr_genes) <- c("ID", "dist.to.feature", "count", "meth.diff", "meth.diff.mean")
  return(dmr_genes)
} 

# Plot aggreated DMR.TSS mean meth.diff vs. LFC
plot_Mean_Meth_Diff_and_LFC <- function(dmr_mean_meth_diff) {
    dmr_mean_meth_diff$LFC <- expr.coef[as.character(dmr_mean_meth_diff$ID), ]$LFC
    t <- subset(dmr_mean_meth_diff, !is.na(LFC))
    plot(t$meth.diff.mean, t$LFC, pch=20, col="blue", xlab='Avg DMR methylation diff', ylab='Log Fold Change',
         xlim=c(min(-50, min(t$meth.diff.mean)-5), max(50, max(t$meth.diff.mean)+5)),
         ylim=c(min(-2, min(t$LFC)-0.25), max(2, max(t$LFC)+0.5)))
    text(t$meth.diff.mean, t$LFC+0.25, labels = t$ID, cex=0.5)
    abline(h=0, col='red')
    abline(v=0, col='red')
    abline(h = 1, lty = 2, col="gray")
    abline(h = -1, lty=2, col="gray")
    return(t)
}

# Plot LogFC vs. methylation diff for given set of genes
plot_DMR_and_LFC <- function(methDiffAll, dmrTSS, ids) { 
  sel.ids <- intersect(dmrTSS$ID, ids)
  sel.TSS <- plotMethyl_and_LFC(methDiffAll, dmrTSS, sel.ids)
  sel.mean <- aggregate_DMR(sel.TSS)
  sel.mean <- plot_Mean_Meth_Diff_and_LFC(sel.mean)
  return(sel.mean)
}


biomine85 <- read.table("../v2/biominer_genes.txt", header=TRUE, sep="\t", as.is = TRUE)
#ids <- intersect(moabs.DMR.TSS$ID, biomine85$ID)
#biomine85.DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, ids)
#biomine85.DMR.mean <- aggregate_DMR(biomine85.DMR.TSS)
#plot_Mean_Meth_Diff_and_LFC(biomine85.DMR.mean)
biomine85.DMR.mean <- plot_DMR_and_LFC(moabs.DMR, moabs.DMR.TSS, biomine85$ID)

# Plot LogFC vs. methylation diff for differentially expressed (DE) genes
annot.DE <- read.csv("../../gene_expression/Ehux_1516_v_217/Jan2019/output/DE_E1516_vs_E217.annot.txt", 
                     header = TRUE, row.names = 1, sep="\t")
DE_ids <- intersect(moabs.DMR.TSS$ID, rownames(annot.DE))
DE.DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, DE_ids)
DE.DMR.mean <- aggregate_DMR(DE.DMR.TSS)
DE.DMR.mean <- plot_Mean_Meth_Diff_and_LFC(DE.DMR.mean)
DE.DMR.mean$Description <- annot.DE[as.character(DE.DMR.mean$ID), "Description"]
DE.DMR.mean.up <- subset(DE.DMR.mean, logFC > 0)
DE.DMR.mean.down <- subset(DE.DMR.mean, logFC < 0)
#write.csv(DE.DMR.mean.up, file = "DE.DMR.mean.up.csv")

# Retrieve the DMRs overlapped with promoter regions
# s4@field is used to extract a field from a S4 object
DE.DMR.TSS.promoters <- subset(DE.DMR.TSS, dist.to.feature <= 1000 & dist.to.feature >= -1000)
plot_DMRs(DE.DMR.TSS.promoters)

# Plot LFC vs. methylation diff for shared up genes
shared_up_ids <- intersect(moabs.DMR.TSS$ID, rownames(shared.up))
shared_up.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, shared_up_ids)
shared_up.mean <- aggregate_DMR(shared_up.TSS)
plot_Mean_Meth_Diff_and_LFC(shared_up.mean)

# Plot LFC vs. methylation diff for carbonic anhydrases
carbonic_DMR_ids <- intersect(moabs.DMR.TSS$ID, rownames(carbonic_anhydrases))
carbonic_DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, carbonic_DMR_ids)
carbonic_DMR.mean <- aggregate_DMR(carbonic_DMR.TSS)
plot_Mean_Meth_Diff_and_LFC(carbonic_DMR.mean)

ehux.DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, moabs.DMR.TSS$ID)
ehux.DMR.mean <- aggregate_DMR(ehux.DMR.TSS)
ehux.DMR.mean <- plot_Mean_Meth_Diff_and_LFC(ehux.DMR.mean)

t <- plot_DMR_Gviz(moabs.DMR, moabs.DMR.TSS, "233823", methDiff = myDiff10)
plot_DMR_Gviz(moabs.DMR, moabs.DMR.TSS, "233304", methDiff = myDiff.DSS)
plot_DMR_Gviz(moabs.DMR, moabs.DMR.TSS, "193771", methDiff = myDiff.DSS)

plotOneDMR(moabs.DMR, moabs.DMR.TSS, "193771", methDiff = myDiff.DSS)

id = "193771"
up_len = 3000
rna = rownames(mapping[mapping$ID == id, ])
TSSes = gene.parts$TSSes
exons = gene.parts$exons
upstream = TSSes[TSSes$name == rna[1], ]
if (as.logical(strand(upstream) == '-')) {
  end(upstream) = end(upstream) + up_len
} else {
  start(upstream) = start(upstream) - up_len
}
ex = exons[exons$name == rna[1], ]
#gene_n_up = GRangesList(upstream, ex)
gene_n_up = c(upstream, ex)
gene_n_up = range(gene_n_up)
strand(gene_n_up) = '*'
methhits <- getMethInGRange(data.raw, gene_n_up)
matplot(methhits$start, methhits$meth.diff, type="o", pch = 20,col = 'blue',
        xlab=paste("Genomic location", methhits$chr[1]), 
        ylab="Methylation difference" )
abline(h=c(-10,0,10),lty=2)


