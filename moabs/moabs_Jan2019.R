library(methylKit)
library(genomation)
library(stringr)

# List of moabs coverage bed files
moabs.bed.list <- list ("./EH1516_r1.moabs.bed",
                        "./EH1516_r2.moabs.bed",
                        "./EH217_r1.moabs.bed",
                        "./EH217_r2.moabs.bed",
                        "./EH217_r3.moabs.bed")

moabsobj <- methRead(moabs.bed.list, 
                     sample.id = list("EH1516_r1", "EH1516_r2", "EH217_r1", "EH217_r2", "EH217_r3"),
                     assembly = "Ehux", 
                     pipeline = list(fraction = FALSE, chr.col=1, start.col=2, end.col=3,
                                    coverage.col=5, freqC.col=4, strand.col=7),
                     treatment = c(1, 1, 0, 0, 0), mincov = 3, header=TRUE)

moabs <- unite(filterByCoverage(moabsobj, lo.count = 10))
getCorrelation(moabs, plot=FALSE)
clusterSamples(moabs, dist="correlation", method="ward", plot=TRUE)
moabsDiff <- calculateDiffMeth(moabs)
