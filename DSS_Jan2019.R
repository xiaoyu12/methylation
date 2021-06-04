######################################################################
# Use the bioconductor DSS for differential methylation analysis
######################################################################
library('DSS')
require(bsseq)

# calculate DMRs using DSS package
# methobj: methyRawList created by methylKit package
# g1 and g2: 1st and 2nd group vectors of sample names
calcDMR_DSS <- function(methobj, g1, g2, mincov = 1, delta = 0.1, minlen = 6, minCG = 2, p.threshold=0.05) {
  meth.list <- list()
  i = 1
  for(mobj in methobj) {
    df <- getData(mobj)
    df <- df[which(df$coverage >= mincov), c(1, 2, 5, 6)]
    colnames(df) <-  c("chr", "pos", "N", "X")
    meth.list[[i]] <- df
    i <- i+1
  }
  # Create BSseq object
  BSobj <- makeBSseqData(meth.list, c(g1, g2))
  dmlTest <- DMLtest(BSobj, group1 = g1, group2 = g2)
  dmrs <- callDMR(dmlTest, delta = delta, minlen = minlen, minCG = minCG, p.threshold = p.threshold)
  newList <- list("BSobj" = BSobj, "DMR" = dmrs)
  return(newList)
}

methobj_DSS_mincov10 <- calcDMR_DSS(methobj, c('EH1516_r1', 'EH1516_r2'), c('EH217_r1', 'EH217_r2', 'EH217_r3'), mincov=10)
methCHG_DSS_mincov10 <- calcDMR_DSS(methCHG, c('EH1516_r1', 'EH1516_r2'), c('EH217_r1', 'EH217_r2', 'EH217_r3'), mincov=10)
save(methCHG_DSS_mincov10, file="methCHG_DSS_mincov10.RData")
meth.CHG.cov10.gr <- makeGRangesFromDataFrame(methCHG_DSS_mincov10$DMR, keep.extra.columns = TRUE)
methCHG_DSS_mincov10.Ann <- annotateWithGeneParts(meth.CHG.cov10.gr, gene.parts)
genomation::plotTargetAnnotation(methCHG_DSS_mincov10.Ann)

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
  data.raw$EH1516 = (cov.raw$EH1516C*data.raw$EH1516C + cov.raw$EH1516D*data.raw$EH1516D)
  data.raw$EH217 = (cov.raw$EH217A*data.raw$EH217A + cov.raw$EH217B*data.raw$EH217B + cov.raw$EH217C*data.raw$EH217C)
  # divide the number of methylated C's with total C's
  data.raw$EH1516 = data.raw$EH1516 / cov.raw$EH1516
  data.raw$EH217 = data.raw$EH217 / cov.raw$EH217
  data.raw$meth.diff = data.raw$EH217 - data.raw$EH1516
  data.out <- data.raw[, c("seqnames", "start", "end", "EH1516", "EH217")]
  colnames(data.out) <- c("V1",	"V2",	"V3",	"beta.1516", "beta.217")
  write.table(data.out, file=file_name, sep="\t", row.names=FALSE, quote=FALSE)
}

# write smoothed BSobj to a file
write_BSobj_sm <- function(BSobj.sm, file_name) {
  data.sm = as.data.frame(getMeth(BSobj.sm))
  grang = as.data.frame(granges(BSobj.sm))
  data.sm <- cbind(grang, data.sm)
  data.sm$EH1516 = (data.sm$EH1516C + data.sm$EH1516D)/2
  data.sm$EH217 = (data.sm$EH217A + data.sm$EH217B + data.sm$EH217C)/3
  data.sm$meth.diff = data.sm$EH217 - data.sm$EH1516
  data.out <- data.sm[, c("seqnames", "start", "end", "EH1516", "EH217")]
  colnames(data.out) <- c("V1",	"V2",	"V3",	"beta.1516", "beta.217")
  write.table(data.out, file=file_name, sep="\t", row.names=FALSE, quote=FALSE)
}

#load methylKit methobj into BSSeq data structure
meth.list <- list() 
i = 1
for(mobj in methobj) {
  df <- getData(mobj)
  df <- df[, c(1, 2, 5, 6)]
  colnames(df) <-  c("chr", "pos", "N", "X")
  meth.list[[i]] <- df
  i <- i+1
}
BSobj <- makeBSseqData(meth.list, c('EH1516C', 'EH1516D', 'EH217A', 'EH217B', 'EH217C'))
#BSSmooth the methylation data
#BSobj.sm <- BSmooth(BSobj,ns = 30, maxGap = 10^6,
#                    BPPARAM = MulticoreParam(workers = parallel::detectCores()-2, 
#                                        progressbar = TRUE))
BSobj.sm <- BSmooth(BSobj, ns=50, h=500, BPPARAM = MulticoreParam(workers = 24, progressbar = TRUE))

BSobj.sm20 <- BSmooth(BSobj, ns=20, h=200, maxGap = 10^6,
                      BPPARAM = MulticoreParam(workers = parallel::detectCores()-2, 
                                               progressbar = TRUE))
BSobj.sm30 <- BSmooth(BSobj, ns=30, h=200, maxGap = 10^6,
                      BPPARAM = MulticoreParam(workers = 24, 
                                               progressbar = TRUE))
#BSobj.sm <- BSobj.sm20

# calculate DE methylation sites and regions
dmlTest <- DMLtest(BSobj, group1 = c('EH1516C', 'EH1516D'), group2 = c('EH217A', 'EH217B', 'EH217C'))
dmlTest.sm <- DMLtest(BSobj, group1 = c('EH1516C', 'EH1516D'), group2 = c('EH217A', 'EH217B', 'EH217C'), smoothing = TRUE)
dmls <- callDML(dmlTest, p.threshold = 0.001)
dmls.sm <- callDML(dmlTest.sm, p.threshold = 0.001)
dmrs <- callDMR(dmlTest, delta = 0.1, minlen = 6, minCG = 2, p.threshold = 0.05)
dmrs.sm <- callDMR(dmlTest.sm, delta = 0.1, minlen = 6, minCG = 2, p.threshold = 0.05)

plot_Methyl_Grange(data.raw, as(dmrs.sm[2, ], "GRanges"))

dmrs.sm$meth.diff <- -100 * dmrs.sm$diff.Methy     # meth.diff = Met
hy2 - Meth1
dmrs.sm.gr <- makeGRangesFromDataFrame(dmrs.sm, ignore.strand = TRUE, keep.extra.columns = TRUE)
dmrs.sm.Ann <- annotateWithGeneParts(dmrs.sm.gr, gene.parts)
dmrs.sm.TSS <- getAssociationWithTSS(dmrs.sm.Ann)
dmrs.sm.TSS$ID <- mapping[dmrs.sm.TSS$feature.name, 2]
dmrs.sm.TSS$LFC <- expr.coef[as.character(dmrs.sm.TSS$ID), ]$LFC
dmrs.sm.TSS$logFC <- expr.coef[as.character(dmrs.sm.TSS$ID), ]$logFC
plotTargetAnnotation(dmrs.sm.Ann)

dmrs$meth.diff <- -100 * dmrs$diff.Methy     # meth.diff = Methy2 - Meth1
dmrs.gr <- makeGRangesFromDataFrame(dmrs, ignore.strand = TRUE, keep.extra.columns = TRUE)
dmrs.Ann <- annotateWithGeneParts(dmrs.gr, gene.parts)
dmrs.TSS <- getAssociationWithTSS(dmrs.Ann)
dmrs.TSS$ID <- mapping[dmrs.TSS$feature.name, 2]
dmrs.TSS$LFC <- expr.coef[as.character(dmrs.TSS$ID), ]$LFC
dmrs.TSS$logFC <- expr.coef[as.character(dmrs.TSS$ID), ]$logFC
plotTargetAnnotation(dmrs.Ann)

# Plot one DMR
showOneDMR(dmrs.sm["1645", ], BSobj)
showOneDMR(dmrs["221935", ], BSobj)

write_BSobj_sm(BSobj.sm, file_name = "meth_bsmooth_avg_data.txt")
write_BSobj_sm(BSobj.sm20, file_name="meth_bsm20_avg_data.txt")
write_BSobj_sm(BSobj.sm30, file_name="meth_bsm30_avg_data.txt")

# raw methylation data
write_BSobj(BSobj)



