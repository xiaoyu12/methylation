##########################################################################
# This file contains code/data not used in the final analysis
#
# It was cut/pasted from other files
#########################################################################

# methCpG isn't used because it's a duplicate of methobj and meth
methCpG <- bismark$methCpG
meth.CpG.cov10 <- methylKit::unite(methCpG)
#PCASamples(meth.CpG.cov10, scale = FALSE)
#getCorrelation(meth.CpG.cov10, plot = FALSE)
#CpG10K <- tileMethylCounts(bismark$methCpG, win.size = 10000, step.size = 10000)
#CpG10K <- unite(CpG10K)
#CpG10K.chr1 <- subset(getData(CpG10K), chr == "chr1")
#CpG10K.chr1$Ehx1516 <- (CpG10K.chr1$numCs1+CpG10K.chr1$numCs2) / 
#  (CpG10K.chr1$coverage1+CpG10K.chr1$coverage2)
#CpG10K.chr1$Ehx217 <- (CpG10K.chr1$numCs3 + CpG10K.chr1$numCs4 + CpG10K.chr1$numCs5) / 
#  (CpG10K.chr1$coverage3 + CpG10K.chr1$coverage4 + CpG10K.chr1$coverage5)
#ggplot(CpG10K.chr1) + geom_line(mapping=aes(x = end/10000, y = Ehx217))


# smoothed BSobj doesn't produce good results
dmlTest.sm <- DMLtest(BSobj, group1 = c('EH1516B', 'EH1516C'), 
                      group2 = c('EH217A', 'EH217B', 'EH217C'), 
                      smoothing = TRUE, BPPARAM=snow)
dmls.sm <- callDML(dmlTest.sm, p.threshold = 0.001)
dmrs.sm <- callDMR(dmlTest.sm, delta = 0.1, minlen = 6, minCG = 2, p.threshold = 0.05)
# Use bsseq dmrFinder to locate DMRs
BSobj.sm.stat <-  BSmooth.tstat(BSobj.sm, group1 = c('EH1516C', 'EH1516D'), 
                                group2 = c('EH217A', 'EH217B', 'EH217C'),
                                mc.cores = 20)
dmr.bsseq <- dmrFinder(BSobj.sm.stat, na.rm=TRUE)

dmrs.sm$meth.diff <- -100 * dmrs.sm$diff.Methy     # meth.diff = Met
hy2 - Meth1
dmrs.sm.gr <- makeGRangesFromDataFrame(dmrs.sm, ignore.strand = TRUE, keep.extra.columns = TRUE)
dmrs.sm.Ann <- annotateWithGeneParts(dmrs.sm.gr, gene.parts)
dmrs.sm.TSS <- getAssociationWithTSS(dmrs.sm.Ann)
dmrs.sm.TSS$ID <- mapping[dmrs.sm.TSS$feature.name, 2]
dmrs.sm.TSS$LFC <- expr.coef[as.character(dmrs.sm.TSS$ID), ]$LFC
dmrs.sm.TSS$logFC <- expr.coef[as.character(dmrs.sm.TSS$ID), ]$logFC
plotTargetAnnotation(dmrs.sm.Ann)
plot_Methyl_Grange(data.raw, as(dmrs.sm[2, ], "GRanges"))

#BSSmooth the methylation data
BSobj.sm <- BSmooth(BSobj,ns = 30, maxGap = 10^6,
                    BPPARAM = MulticoreParam(workers = 20, progressbar = TRUE))
BSobj.sm <- BSmooth(BSobj, ns=50, h=500, BPPARAM = MulticoreParam(workers = 24, progressbar = TRUE))

BSobj.sm20 <- BSmooth(BSobj, ns=20, h=200, maxGap = 10^6,
                      BPPARAM = MulticoreParam(workers = parallel::detectCores()-2, 
                                               progressbar = TRUE))
BSobj.sm30 <- BSmooth(BSobj, ns=30, h=200, maxGap = 10^6,
                      BPPARAM = MulticoreParam(workers = 24, 
                                               progressbar = TRUE))
#BSobj.sm <- BSobj.sm20

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

write_BSobj_sm(BSobj.sm, file_name = "meth_bsmooth_avg_data.txt")
write_BSobj_sm(BSobj.sm20, file_name="meth_bsm20_avg_data.txt")
write_BSobj_sm(BSobj.sm30, file_name="meth_bsm30_avg_data.txt")

# DSS object is named BSobj
#methobj_DSS_mincov10 <- calcDMR_DSS(methobj, c('EH1516B', 'EH1516C'), c('EH217A', 'EH217B', 'EH217C'), mincov=10)

DM_count_by_region <- function(region.gr, min_len=500) {
  dmc.count = as.data.frame(region.gr)
  dmc.count$len = dmc.count$end - dmc.count$start + 1
  
  dmc.count$onoff_dmc = countOverlaps(region.gr, meth.dmc.gr, type="any", ignore.strand=TRUE)
  # count overlaps with moabs DMCs
  dmc.count$moabs = countOverlaps(region.gr, moabs.DMC, type = "any", ignore.strand=TRUE)
  
  # Count promoter overlapped with moabs.DMRs
  dmc.count$dmrs = countOverlaps(region.gr, moabs.DMR, type="any", ignore.strand=TRUE)
  
  # Count promoter overlapped with moabs CHG DMCs
  dmc.count$hg.dmcs = countOverlaps(region.gr, moabs.HG.DMC, type="any", ignore.strand=TRUE)
  dmc.count$hg.dmrs = countOverlaps(region.gr, moabs.HG.DMR, type="any", ignore.strand=TRUE)
  
  # separate hyper and hypo-methylated DMCs
  dmc.count$dmc.hypo = countOverlaps(region.gr, moabs.DMC.hypo, type="any", ignore.strand=TRUE)
  dmc.count$dmc.hyper = countOverlaps(region.gr, moabs.DMC.hyper, type="any", ignore.strand=TRUE)
  
  # filter out regions shorter than min_len
  dmc.count = dmc.count %>% filter(len >= min_len)
  
  return(dmc.count)
}
