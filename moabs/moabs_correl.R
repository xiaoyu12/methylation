library(methylKit)
library(genomation)
library(stringr)

#bismark.files <- list("./EH1516_r1.bismark.txt", "./EH1516_r2.bismark.txt", "./EH217_r1.bismark.txt",
#                      "./EH217_r2.bismark.txt", "./EH217_r3.bismark.txt")

bismark.files <- list("../EH1516C.merged_CpG_evidence.cov.gz",
                        "../EH1516D.merged_CpG_evidence.cov.gz",
                        "../EH217A.merged_CpG_evidence.cov.gz",
                        "../EH217B.merged_CpG_evidence.cov.gz",
                        "../EH217C.merged_CpG_evidence.cov.gz"
                        )

moabs.bed.files <- list("D:/xiaoyu/Documents/moabs_rudy/EH1516_r1.bismark.txt",
                        "D:/xiaoyu/Documents/moabs_rudy/EH1516_r2.bismark.txt",
                        "D:/xiaoyu/Documents/moabs_rudy/EH217_r1.bismark.txt",
                        "D:/xiaoyu/Documents/moabs_rudy/EH217_r2.bismark.txt",
                        "D:/xiaoyu/Documents/moabs_rudy/EH217_r3.bismark.txt")

all.list <- list ("../EH1516C.merged_CpG_evidence.cov.gz",
                  "../EH1516D.merged_CpG_evidence.cov.gz",
                  "../EH217A.merged_CpG_evidence.cov.gz",
                  "../EH217B.merged_CpG_evidence.cov.gz",
                  "../EH217C.merged_CpG_evidence.cov.gz",
                  "D:/xiaoyu/Documents/moabs_rudy/EH1516_r1.bismark.txt",
                  "D:/xiaoyu/Documents/moabs_rudy/EH1516_r2.bismark.txt",
                  "D:/xiaoyu/Documents/moabs_rudy/EH217_r1.bismark.txt",
                  "D:/xiaoyu/Documents/moabs_rudy/EH217_r2.bismark.txt",
                  "D:/xiaoyu/Documents/moabs_rudy/EH217_r3.bismark.txt")

samples <- list("EH1516_r1", "EH1516_r2", "EH217_r1", "EH217_r2", "EH217_r3")

moabsobj <- methRead(moabs.bed.files, 
                     sample.id = samples,
                     assembly = "Ehux", 
                     pipeline = "bismarkCoverage",
                     treatment = c(1, 1, 0, 0, 0), mincov = 10, header=TRUE)

methobj <- methRead(bismark.files, sample.id = list("EH1516C", "EH1516D", "EH217A", "EH217B", "EH217C"),
                    assembly = "Ehux", pipeline = "bismarkCoverage", #list(fraction = TRUE, chr.col=1, start.col=2, end.col=3, coverage.col=5, freqC.col=4, strand.col=7), 
                    treatment = c(1, 1, 0, 0, 0),
                    context = "CpG", mincov = 10, header = TRUE)

moabs <- unite(filterByCoverage(moabsobj, lo.count = 10))
meth <- unite(filterByCoverage(methobj, lo.count = 10))
getCorrelation(moabs, plot=FALSE)
getCorrelation(meth, plot=FALSE)

#allobj <- methobj
#for (i in 1:5) {
#  allobj[[i+5]] <- moabsobj[[i]]
#}

allobj <- methRead(all.list, 
                   sample.id = list("EH1516C", "EH1516D", "EH217A", "EH217B", "EH217C", "EH1516_r1", "EH1516_r2", "EH217_r1", "EH217_r2", "EH217_r3"),
                   assembly = "Ehux", 
                   pipeline = "bismarkCoverage",
                   treatment = c(1, 1, 0, 0, 0, 1, 1, 0, 0, 0), mincov = 10, header=TRUE)

all <- unite(filterByCoverage(allobj, lo.count = 10))

getCorrelation(all, plot = FALSE)

clusterSamples(all, dist="correlation", method="ward", plot=TRUE)
