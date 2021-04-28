library(methylKit)
library(genomation)
library(stringr)
library(Biobase)

# List of CpG coverage files
file.list <- list("../v2/1516-FCH.merged_CpG_evidence2.cov.gz",
                  "../v2/217-FCH.merged_CpG_evidence2.cov.gz",
                  "EH1516C.merged_CpG_evidence.cov.gz",
                  "EH1516D.merged_CpG_evidence.cov.gz",
                  "EH217A.merged_CpG_evidence.cov.gz",
                  "EH217B.merged_CpG_evidence.cov.gz",
                  "EH217C.merged_CpG_evidence.cov.gz")


# Read Ehux gene structures
gene.parts <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, unique.prom = FALSE)
# Promoters defined as up.flank = 1000, down.flank = 0
gene.parts_up1000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, 
                                            up.flank = 1000, down.flank = 0, unique.prom = FALSE)
# Promoters defined as up.flank = 0, down.flank = 1000
gene.parts_dn1000 <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, 
                                            up.flank = 0, down.flank = 1000, unique.prom = FALSE)

sample.ids = c("EH1516-FCH", "EH217-FCH", "EH1516C", "EH1516D", "EH217A", "EH217B", "EH217C")
# Read methylation data
methobj <- methRead(file.list, sample.id = as.list(sample.ids),
                    assembly = "ehux", pipeline = "bismarkCoverage", treatment = c(1, 0, 1, 1, 0, 0, 0),
                    context = "CpG", mincov = 3, header = FALSE)

#Unite methylRawList to a single table
meth <- unite(methobj)
rm(methobj)

myDiff.DSS <- calculateDiffMethDSS(meth, mc.cores=6)

getCorrelation(meth,plot=FALSE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, scale = FALSE)

meth10 <- unite(filterByCoverage(methobj, lo.count = 10))
getCorrelation(meth10, plot = FALSE)
clusterSamples(meth10, dist="correlation", plot=TRUE)
PCASamples(meth10, scale = FALSE)

# Read mapping from rna* to gene ID
mapping <- read.table("../v2/genbank_mapping.txt", row.names = 1)
colnames(mapping) <- c("gene", "ID")


# Get methylation values for gene annotated features: e.g. promoters, exons ...
getFeatureMethyl <- function(m.obj, g.features, s.names) {
  nsamples <- length(m.obj)
  promoters <- as.data.frame(g.features)
  promoters$ID = mapping[promoters$name, "ID"]
  promoters = promoters[!duplicated(promoters$ID), ]    # remove duplicates due to alternative splicing
  colnames(promoters)[1] <- "chr"
  
  promobj <- regionCounts(m.obj, g.features)

  meth.prom <- list()
  for (i in 1:nsamples) {
    t <- merge(getData(promobj[[i]]), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
    t <- t[!duplicated(t$ID), ]    # Remove duplicates due to alternative splicing
    t <- t[t$coverage >= 100, ]
    t$beta <- t$numCs / t$coverage
    t$m <- log2((t$numCs+1)/(t$numTs+1))
    rownames(t) <- t$ID
    meth.prom[[i]] <- t
  }
  
  ids <- meth.prom[[1]]$ID
  for (i in 2:nsamples) {
    ids <- intersect(ids, meth.prom[[i]]$ID)
  }
  mtx.coverage <- meth.prom[[1]][as.character(ids), "coverage"]
  mtx.m <- meth.prom[[1]][as.character(ids), "m"]
  mtx.beta <- meth.prom[[1]][as.character(ids), "beta"]
  for (i in 2:nsamples) {
    mtx.coverage <- cbind(mtx.coverage, meth.prom[[i]][as.character(ids), "coverage"])
    mtx.m <- cbind(mtx.m, meth.prom[[i]][as.character(ids), "m"])
    mtx.beta <- cbind(mtx.beta, meth.prom[[i]][as.character(ids), "beta"])
  }
  colnames(mtx.coverage) <- s.names
  colnames(mtx.m) <- s.names
  colnames(mtx.beta) <- s.names
  rownames(mtx.coverage) <- as.character(ids)
  rownames(mtx.m) <- as.character(ids)
  rownames(mtx.beta) <-as.character(ids)
  
  return (list("coverage"=mtx.coverage, "beta"=mtx.beta, "m"=mtx.m))
}


m <- getFeatureMethyl(methobj, gene.parts$promoters, sample.ids)
m_up1000 <- getFeatureMethyl(methobj, gene.parts_up1000$promoters, sample.ids)
m_dn1000 <- getFeatureMethyl(methobj, gene.parts_dn1000$promoters, sample.ids)

m.EH1516.beta <- m$beta[, 3:4]
m.EH217.beta <- m$beta[, 5:7]

m.DMR.beta <- m$beta[which(abs(rowMedians(m.EH1516.beta) - rowMedians(m.EH217.beta)) >= 0.2), ]

# read Ehux gene annotations
annotations <- read.csv("../v2/ehux_JGI_best_b2g.txt", header=TRUE, sep="\t")
annot.DE <- read.csv("DE_E1516_vs_E217.annot.txt", header = TRUE, row.names = 1, sep="\t")

annot.DE.DMR <- annot.DE[intersect(rownames(annot.DE), rownames(m.DMR.beta)), ]

m10 <- getFeatureMethyl(filterByCoverage(methobj, lo.count = 10), gene.parts$promoters, sample.ids)
m10_up1000 <- getFeatureMethyl(filterByCoverage(methobj, lo.count = 10), gene.parts_up1000$promoters, sample.ids)

m_exon <- getFeatureMethyl(methobj, gene.parts$exons, sample.ids)
m_intron <- getFeatureMethyl(methobj, gene.parts$introns, sample.ids)

meth.prom <- unite(filterByCoverage(promobj, lo.count = 100))
getCorrelation(meth.prom, plot = FALSE)
meth.prom.DSS <- calculateDiffMeth(meth.prom)
meth.prom.DSS.20p <- getMethylDiff(meth.prom.DSS, difference = 20, qvalue = 0.01)
promoters <- as.data.frame(gene.parts$promoters)
promoters$ID = mapping[promoters$name, "ID"]
promoters = promoters[!duplicated(promoters$ID), ]    # remove duplicates due to alternative splicing
colnames(promoters)[1] <- "chr"
meth.prom.DSS <- merge(getData(meth.prom.DSS), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
meth.prom.DSS.20p <- merge(getData(meth.prom.DSS.20p), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
