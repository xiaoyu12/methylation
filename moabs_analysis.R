library(methylKit)
library(genomation)
library(stringr)
source("meth_analysis_funcs.R")

# read MOABS DMR from a file (backup function if not in meth_analysis_funcs.R)
if (!exists("readMOABS_DMR_bedfile")) {
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
}

# Read MOABS DMR and DMC files
moabs.DMR <- readMOABS_DMR_bedfile("dmr_M2_EH1516.G.bed_vs_EH217.G.bed.bed")
moabs.DMC <- readMOABS_DMR_bedfile("dmc_M2_EH1516.G.bed_vs_EH217.G.bed.bed")
strand(moabs.DMR) = "*"
strand(moabs.DMC) = "*"

# change strand to *


# Read MOABS DMRs from moabs1.3.4 analysis
moabs.134 <- readMOABS_DMR_bedfile("dmr_M2_moabs1.3.4.bed")
ov <- findOverlaps(moabs.134, moabs.DMR)
length(unique(from(ov)))
length(unique(to(ov)))
# Only keep the overlapped ones. Others are susceptible based on manual inspection
moabs.DMR <- moabs.134[unique(from(ov)), ]
# convert DMR to a dataframe
moabs.DMR.df <- as.data.frame(moabs.DMR)
colnames(moabs.DMR.df)[1] <- "chr"


# Read MOABS CHG DMRs 
# DMRs are computed using converted bismark files 
# Data are from the directory /bigdata/xiaoyu/bisulfite/Jan2019/moabs_bismark on the server
moabs.HG.DMC <- readMOABS_DMR_bedfile("dmc_M2_EH1516.HG.bed_vs_EH217.HG.bed.bed")
moabs.HG.DMR <- readMOABS_DMR_bedfile("dmr_M2_EH1516.HG.bed_vs_EH217.HG.bed.bed")
moabs.HG.DMR.df <- as.data.frame(moabs.HG.DMR)
colnames(moabs.HG.DMR.df)[1] <- "chr"
moabs.HG.DMR.l100 <- subset(moabs.HG.DMR.df, width >= 100)
showOneDMR(moabs.HG.134.l100[2, ], methCHG_DSS_mincov10$BSobj)

# combine MOABS CpG and CHG DMRs into a genomicrange list
moabs.all <- c(moabs.DMR, moabs.HG.DMR)

# Find overlaps of DMRs and genes
# Removed previous incorrect referrence to CPGi as Ehux_genbank.bed doesn't have that info
# It only has gene structures. 

ov <- findOverlaps(moabs.all, genes$Gene, ignore.strand=TRUE)
DMGs <- genes$Gene[unique(to(ov)), ]
DMGs$ID <- mapping[DMGs$name, ]$ID
ov <- findOverlaps(moabs.DMR, genes$Gene, ignore.strand=TRUE)
DMGs.CpG <- genes$Gene[unique(to(ov)), ]
DMGs.CpG$ID <- mapping[DMGs.CpG$name,]$ID

# Plot DMR values and gene LFC
ov <- findOverlaps(moabs.DMR[width(moabs.DMR) >= 100], genes$Gene, ignore.strand=TRUE)
DMR.ov.genes <- as.data.frame(ov) 
DMR.ov.genes <- cbind(DMR.ov.genes, as.data.frame(moabs.DMR[from(ov)]))
DMR.ov.genes <- cbind(DMR.ov.genes, as.data.frame(genes$Gene[to(ov)]))
colnames(DMR.ov.genes)[8] <- "dmr.name"
colnames(DMR.ov.genes)[3:6] <- c("dmr.seqnames", "dmr.start", "dmr.end", "dmr.width")
DMR.ov.genes[7] <- NULL
DMR.ov.genes$ID <- mapping[DMR.ov.genes$name, ]$ID
DMR.ov.genes$logFC <- expr.coef[as.character(DMR.ov.genes$ID), "logFC"]
png(filename = "moabs_graphs/DMGs_Methdiff_LFG.png", res=600)
ggplot(DMR.ov.genes) + geom_point(aes(x=meth.diff, y=logFC), color="blue", shape=20) + theme_light() 
dev.off()

ov <- findOverlaps(moabs.HG.DMR, genes$Gene, ignore.strand=TRUE)
DMGs.CHG <- genes$Gene[unique(to(ov)), ]
DMGs.CHG$ID <- mapping[DMGs.CHG$name, ]$ID

# extend gene regions upstream 1000 bp
genes_up1000 <- resize(genes, width(genes) + 1000L, fix = "end")

#write.csv(DMGs.CpG$ID, file="DMGs.CpG.csv", quote=FALSE)
#write.csv(DMGs.CHG$ID, file="DMGs.CHG.csv", quote=FALSE)
#write.csv(DMGs$ID, file="DMGs.csv", quote = FALSE)

# Find the overlaps of DMRs and genes upstream regions
up2000 <- intersect(up2000, scaffolds.gr)
ov <- findOverlaps(moabs.all, up2000, ignore.strand=TRUE)
DMPs <- up2000[unique(to(ov)), ]
ov <- findOverlaps(moabs.DMR, up2000, ignore.strand=TRUE)
DMPs.CpG <- up2000[unique(to(ov)), ]
ov <- findOverlaps(moabs.HG.DMR, up2000, ignore.strand=TRUE)
DMPs.CHG <- up2000[unique(to(ov)), ]
#write.csv(DMPs.CpG$ID, file="DMPs.CpG.csv", quote=FALSE)
#write.csv(DMPs.CHG$ID, file="DMPs.CHG.csv", quote=FALSE)
#write.csv(DMPs$ID, file="DMPs.csv", quote = FALSE)

annot.up1000 <- annotateWithFeature(moabs.DMR, up1000)
genomation::plotTargetAnnotation(annot.up1000)
annot.down1000 <- annotateWithFeature(moabs.DMR, down1000)
genomation::plotTargetAnnotation(annot.down1000)

moabs.DMR.Ann <- annotateWithGeneParts(moabs.DMR, gene.parts, intersect.chr = TRUE)
moabs.DMR.TSS <- getAssociationWithTSS(moabs.DMR.Ann)
moabs.DMR.TSS$ID <- mapping[moabs.DMR.TSS$feature.name, 2]
moabs.DMR.TSS$logFC <- expr.coef[as.character(moabs.DMR.TSS$ID), ]$logFC
genomation::plotTargetAnnotation(moabs.DMR.Ann)

# MOABS DMR overlapping with genes and flanking regions
moabs.DMR.genes <- annotateWithFeatureFlank(moabs.DMR, genes$Gene,genes$UTR,
                                           feature.name="Gene",flank.name="UTR", intersect.chr = TRUE)
genomation::plotTargetAnnotation(moabs.DMR.genes,
                                 col=c("green","gray","white"),main="DMRs overlapped with genes")


#moabs.DMR.in.genes <- moabs.DMR.TSS[which(slot(moabs.DMR.genes, "members")[, 1] > 0), ] # not correct because moabs.DMR.cpgi has longer length than moabs.DMR.TSS
#moabs.DMR.in.genes$meth.diff <- as.data.frame(moabs.DMR)[moabs.DMR.in.genes$target.row, "meth.diff"]
#plot(moabs.DMR.in.genes$meth.diff, moabs.DMR.in.genes$logFC, cex=1, pch = 20, 
#     col= "blue", xlab = "Methylation Difference", ylab = "Log Fold Change")

# Plot longest MOABS DMRs
moabs.longDMR <- subset(moabs.DMR.df, width > 1000)
colnames(moabs.longDMR)[1] <- "chr"
showOneDMR(moabs.longDMR[8, ], BSobj)

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
    plot(t$meth.diff.mean, t$LFC, pch=20, cex = 0.8, col="blue", xlab='Avg DMR methylation diff', ylab='Log Fold Change',
         xlim=c(min(-50, min(t$meth.diff.mean)-5), max(50, max(t$meth.diff.mean)+5)),
         ylim=c(min(-2, min(t$LFC)-0.25), max(2, max(t$LFC)+0.5)))
    #text(t$meth.diff.mean, t$LFC+0.25, labels = t$ID, cex=0.5)
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
annot.DE <- read.csv("../../gene_expression/Ehux_1516_v_217/Jan2019/output/DE_E1516_vs_E217_FDR0.05.annot.txt", 
                     header = TRUE, row.names = 1, sep="\t")
DE_ids <- intersect(moabs.DMR.TSS$ID, rownames(annot.DE))
DE.DMR.TSS <- plotMethyl_and_LFC(moabs.DMR, moabs.DMR.TSS, DE_ids)
DE.DMR.mean <- aggregate_DMR(DE.DMR.TSS)
DE.DMR.mean <- plot_Mean_Meth_Diff_and_LFC(DE.DMR.mean)
DE.DMR.mean$Description <- annot.DE[as.character(DE.DMR.mean$ID), "Description"]
DE.DMR.mean.up <- subset(DE.DMR.mean, LFC > 0)
DE.DMR.mean.down <- subset(DE.DMR.mean, LFC < 0)
#write.csv(DE.DMR.mean.up, file = "DE.DMR.mean.up.csv")

# Plot the LFC and Avg DMR methylation difference for DMGs
png(filename = "moabs_graphs/DMGs_DMR_LFC.png", width=7800, height=5400, res=600)
plot_DMR_and_LFC(moabs.DMR, moabs.DMR.TSS, DMGs.CpG$ID)
dev.off()

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

# Plot MOABS DMR near a selected gene
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

# Find the number of overlapped genes between DMGs and DEGs
non_DMGs_IDs <- setdiff(genes$Gene$ID, DMGs$ID)
non_DEGs_IDs <- setdiff(genes$Gene$ID, rownames(annot.DE))
DMGs_IDs <- DMGs$ID
DEGS_IDs <- rownames(annot.DE)
n_11 <- length(intersect(DMGs_IDs, DEGS_IDs))
n_12 <- length(intersect(DMGs_IDs, non_DEGs_IDs))
n_21 <- length(intersect(non_DMGs_IDs, DEGS_IDs))
n_22 <- length(intersect(non_DMGs_IDs, non_DEGs_IDs))
DMG_n_DEG <- matrix(c(n_11, n_12, n_21, n_22), nrow=2, byrow = TRUE) 
dimnames(DMG_n_DEG) <- list("DMG" = c("DMG_Y", "DMG_N"), 
                       "DEG" = c("DEG_Y", "DEG_N"))
# Perform Chi-square test of independence
test.DMG_n_DEG <- chisq.test(DMG_n_DEG)

DMPs_IDs <- DMPs$ID
non_DMPs_IDs <- setdiff(genes$Gene$ID, DMPs_IDs)
n_11 <- length(intersect(DMPs_IDs, DEGS_IDs))
n_12 <- length(intersect(DMPs_IDs, non_DEGs_IDs))
n_21 <- length(intersect(non_DMPs_IDs, DEGS_IDs))
n_22 <- length(intersect(non_DMPs_IDs, non_DEGs_IDs))
DMP_n_DEG <- matrix(c(n_11, n_12, n_21, n_22), nrow=2, byrow = TRUE) 
dimnames(DMP_n_DEG) <- list("DMP" = c("DMP_Y", "DMP_N"), 
                            "DEG" = c("DEG_Y", "DEG_N"))
# Perform Chi-square test of independence
test.DMP_n_DEG <- chisq.test(DMP_n_DEG)