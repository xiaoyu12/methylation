######################################################################
# Use the bioconductor DSS for differential methylation analysis
######################################################################
library('DSS')
require(bsseq)
source("meth_analysis_funcs.R")

# create a BSSeq object from a methyRawList 
# methobj: methylBase created by methylKit unite function
# g1 and g2: 1st and 2nd group vectors of sample names
makeBSobj<- function(meth, g1, g2, mincov = 3) {
  meth.list <- list()
  df <- getData(meth)
  for (i in 1:5) {
    print(i)
    df_i <- df[, c(1, 2, 3*i+2, 3*i+3)]
    colnames(df_i) <- c("chr", "pos", "N", "X")
    meth.list[[i]] <-df_i
  } 
  
  # Create BSseq object
  BSobj <- makeBSseqData(meth.list, c(g1, g2))
  
  rm(meth.list)
  rm(df_i)
  rm(df)
  
  return(BSobj)
}

# Use DSS to calcualte differentially methylated sites and regions 
calcDMR_DSS <- function (BSobj, g1, g2, delta = 0.1, 
                         minlen = 6, minCG = 3, p.threshold=0.01,
                         cpu.cores = 20) {

  #pa <- SnowParam(workers = cpu.cores, type = "SOCK", progressbar = TRUE)
  pa <- MulticoreParam(workers = cpu.cores, progressbar = TRUE)
  dmlTest <- DMLtest(BSobj, group1 = g1, group2 = g2, BPPARAM=pa)
  dmcs <- callDML(dmlTest, delta=0.25, p.threshold=1e-5)
  dmrs <- callDMR(dmlTest, delta = delta, minlen = minlen, 
                  minCG = minCG, p.threshold = p.threshold)
  newList <- list("BSobj" = BSobj, "dml" = dmlTest, "DMR" = dmrs, "DMC"=dmcs)
  return (newList)
}

#save(methCHG_DSS_mincov10, file="methCHG_DSS_mincov10.RData")
# Read Ehux gene structures, using up 3000 and down 0 for promoter regions
gene.parts.up3000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE,
                                            up.flank = 3000, down.flank = 0, unique.prom = FALSE)
gene.parts.up2000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE,
                                            up.flank = 2000, down.flank = 0, unique.prom = FALSE)

# Get overlaps of sigDiff10 and gene.parts with up 3000 promoters
sigDiff10.o.up3000 <- overlap_DMR_GeneFeatures(sigDiff10, gene.parts.up3000, expr.coef)

meth.CHG.cov10.gr <- makeGRangesFromDataFrame(methCHG_DSS_mincov10$DMR, keep.extra.columns = TRUE)
methCHG_DSS_mincov10.Ann <- annotateWithGeneParts(meth.CHG.cov10.gr, gene.parts.up3000)
genomation::plotTargetAnnotation(methCHG_DSS_mincov10.Ann)

methCHG_DSS_mincov10$DMR$meth.diff <- -100 * methCHG_DSS_mincov10$DMR$diff.Methy
dmr_CHG <- overlap_DMR_GeneFeatures(methCHG_DSS_mincov10$DMR, gene.parts, expr.coef, 
                                    selected_ids = rownames(annot.DE))
plot(dmr_CHG$sel$meth.diff, dmr_CHG$sel$logFC, type="p", cex = 1, pch = 20, col= "blue", xlab = "Methylation Difference", ylab = "Log Fold Change")
text(dmr_CHG$sel$meth.diff, dmr_CHG$sel$logFC + 0.25, labels = dmr_CHG$sel$ID, cex = 0.5)

#write BSobj to a txt file
write_BSobj <- function(BSobj, file_name="meth_raw_avg_data.txt") {
  # raw methylation data
  data.raw = as.data.frame(getMeth(BSobj, type='raw'))
  cov.raw = as.data.frame(getCoverage(BSobj))
  grang = as.data.frame(granges(BSobj))
  data.raw = cbind(grang, data.raw)
  library(tidyr)
  # convert NA in methylation data to 0
  data.raw[is.na(data.raw)] <- 0
  #data.raw <- drop_na(data.raw)
  cov.raw$EH1516 = (cov.raw$EH1516C + cov.raw$EH1516D)
  cov.raw$EH1516[cov.raw$EH1516 == 0] <- 1
  cov.raw$EH217 = (cov.raw$EH217A + cov.raw$EH217B + cov.raw$EH217C)
  cov.raw$EH217[cov.raw$EH217 == 0] <- 1
  # number of methylated C's
  data.raw$EH1516 = (cov.raw$EH1516C*data.raw$EH1516C + 
                       cov.raw$EH1516D*data.raw$EH1516D)
  data.raw$EH217 = (cov.raw$EH217A*data.raw$EH217A + 
                      cov.raw$EH217B*data.raw$EH217B + 
                      cov.raw$EH217C*data.raw$EH217C)
  # divide the number of methylated C's with total C's
  data.raw$EH1516 = data.raw$EH1516 / cov.raw$EH1516
  data.raw$EH217 = data.raw$EH217 / cov.raw$EH217
  data.raw$meth.diff = data.raw$EH217 - data.raw$EH1516
  data.out <- data.raw[, c("seqnames", "start", "end", "EH1516", "EH217")]
  colnames(data.out) <- c("V1",	"V2",	"V3",	"beta.1516", "beta.217")
  write.table(data.out, file=file_name, sep="\t", row.names=FALSE, quote=FALSE)
}



#load methylKit methobj into BSSeq data structure
meth.list <- list() 
# Only consider sites with a min coverage of 3 in all samples.
df <- getData(meth)
#df <- getData(meth.CpG.cov10)
for (i in 1:5) {
  print(i)
  df_i <- df[, c(1, 2, 3*i+2, 3*i+3)]
  colnames(df_i) <- c("chr", "pos", "N", "X")
  meth.list[[i]] <-df_i
} 
# BSobj for CpG context
BSobj <- makeBSseqData(meth.list, c('EH1516B', 'EH1516C', 'EH217A', 'EH217B', 'EH217C'))

snow <- SnowParam(workers = 24, type = "SOCK", progressbar = TRUE)

# calculate DE methylation sites and regions
# smoothed results aren't good. We use the unsmoothed data 
dmlTest <- DMLtest(BSobj, group1 = c('EH1516B', 'EH1516C'), 
                   group2 = c('EH217A', 'EH217B', 'EH217C'),
                   BPPARAM = snow)

dmls <- callDML(dmlTest, p.threshold = 0.001)
# meth.diff is the methylation rate difference (EH217 - EH1516)
dmls$meth.diff = (dmls$mu2 - dmls$mu1)*100
library(tibble)
dmls = add_column(dmls, dmls$pos+1, .after="pos")
colnames(dmls)[2:3] = c("start", "end")

#dmrs <- callDMR(dmlTest, delta = 0.1, minlen = 6, minCG = 2, p.threshold = 0.05)
# more stringent criteria for DMRs
dmrs <- callDMR(dmlTest, minlen = 10);
# filter out dmrs with nCG density < 0.02
dmrs <- dmrs %>% filter (nCG / length > 0.02)

# BSobj for CHG context
meth.list <- list()
df <- getData(meth.CHG.cov3)
for (i in 1:5) {
  print(i)
  df_i <- df[, c(1, 2, 3*i+2, 3*i+3)]
  colnames(df_i) <- c("chr", "pos", "N", "X")
  meth.list[[i]] <-df_i
} 
BSobj.CHG <- makeBSseqData(meth.list, c('EH1516B', 'EH1516C', 'EH217A', 'EH217B', 'EH217C'))
dmlTest.CHG <- DMLtest(BSobj.CHG, group1 = c('EH1516B', 'EH1516C'), 
                       group2 = c('EH217A', 'EH217B', 'EH217C'),
                       BPPARAM = snow)
# Find the difference between two dataframes that represent genomic regions
diffRegions <- function(df1, df2) {
  gr1 = as(df1, "GRanges")
  gr2 = as(df2, "GRanges")
  o <- findOverlaps(gr1, gr2)
  hits <- unique(queryHits(o))
  nohits <- setdiff(1:length(gr1), hits)
  return(df1[nohits, ])
}


dmrs$meth.diff <- -100 * dmrs$diff.Methy     # meth.diff = Methy2 - Meth1
dmrs.gr <- makeGRangesFromDataFrame(dmrs, ignore.strand = TRUE, keep.extra.columns = TRUE)
dmrs.Ann <- annotateWithGeneParts(dmrs.gr, gene.parts)
dmrs.TSS <- getAssociationWithTSS(dmrs.Ann)
dmrs.TSS$ID <- mapping[dmrs.TSS$feature.name, 2]
dmrs.TSS$LFC <- expr.coef[as.character(dmrs.TSS$ID), ]$LFC
dmrs.TSS$logFC <- expr.coef[as.character(dmrs.TSS$ID), ]$logFC
genomation::plotTargetAnnotation(dmrs.Ann)

# Compare DMRs to MOABS results
dmrs.nohits <- diffRegions(dmrs, moabs.DMR.df)
moabs.nohits <- diffRegions(moabs.DMR.df, dmrs)
showOneDMR(dmrs.nohits[1, ], BSobj)

# Plot one DMR
showOneDMR(dmrs["221935", ], BSobj)

dmrs.long <- subset(dmrs, length > 1000)
showOneDMR(dmrs.long[1, ], BSobj)


dmrs.CHG <- callDMR(dmlTest.CHG, minlen = 10)
dmls.CHG <- callDML(dmlTest.CHG, p.threshold = 0.001)
dmls.CHG$meth.diff = (dmls.CHG$mu2 - dmls.CHG$mu1)*100
dmls.CHG = add_column(dmls.CHG, dmls.CHG$pos, .after="pos")
colnames(dmls.CHG)[2:3] = c("start", "end")

# filter out dmrs with nCG density < 0.02
dmrs.CHG <- dmrs.CHG %>% filter (nCG / length > 0.02)

methCHG_DSS_mincov10$BSobj <- makeBSobj(methCHG, g1 = c('EH1516B', 'EH1516C'),
                                        g2 = c('EH217A', 'EH217B', 'EH217C'), 
                                        mincov=3)
methCHG_DSS_mincov10 <- calcDMR_DSS(methCHG_DSS_mincov10$BSobj, 
                                    c('EH1516B', 'EH1516C'), 
                                    c('EH217A', 'EH217B', 'EH217C'),
                                    cpu.cores = 20)
dmrs.CHG.nohits <- diffRegions(dmrs.CHG, methCHG_DSS_mincov10$DMR)
showOneDMR(dmrs.CHG.nohits[1, ], BSobj.CHG)
t <- diffRegions(moabs.HG.DMR.df, dmrs.CHG)

dmrs.CHG.Ann <- annotateWithGeneParts(as(dmrs.CHG, "GRanges"), gene.parts)
par(mfrow=c(1, 1))
genomation::plotTargetAnnotation(dmrs.CHG.Ann)
moabs.HG.DMR.Ann = annotateWithGeneParts(moabs.HG.DMR, gene.parts)
genomation::plotTargetAnnotation(moabs.HG.DMR.Ann)

# raw methylation data
write_BSobj(BSobj)



