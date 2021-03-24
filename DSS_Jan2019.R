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

meth.list <- list() 
i = 1
for(mobj in methobj) {
  df <- getData(mobj)
  df <- df[, c(1, 2, 5, 6)]
  colnames(df) <- c("chr", "pos", "N", "X")
  meth.list[[i]] <- df
  i <- i+1
}
BSobj <- makeBSseqData(meth.list, c('EH1516C', 'EH1516D', 'EH217A', 'EH217B', 'EH217C'))
dmlTest <- DMLtest(BSobj, group1 = c('EH1516C', 'EH1516D'), group2 = c('EH217A', 'EH217B', 'EH217C'))
dmlTest.sm <- DMLtest(BSobj, group1 = c('EH1516C', 'EH1516D'), group2 = c('EH217A', 'EH217B', 'EH217C'), smoothing = TRUE)
dmls <- callDML(dmlTest, p.threshold = 0.001)
dmls.sm <- callDML(dmlTest.sm, p.threshold = 0.001)
dmrs <- callDMR(dmlTest, delta = 0.1, minlen = 6, minCG = 2, p.threshold = 0.05)
dmrs.sm <- callDMR(dmlTest.sm, delta = 0.1, minlen = 6, minCG = 2, p.threshold = 0.05)

dmrs.sm$meth.diff <- -100 * dmrs.sm$diff.Methy     # meth.diff = Methy2 - Meth1
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
region <- data.frame(chr = "chr67", start = 429356, end = 433135, length=3780)
showOneDMR(region, BSobj)
showOneDMR(dmrs["221935", BSobj])

# Plot LogFC vs. methylation diff for differentially expressed (DE) genes
biomine85.dmrs.mean <- plot_DMR_and_LFC(dmrs.gr, dmrs.TSS, biomine85$ID)
