############################################################################
# This script analyzes CpG methylation data of E. hux
############################################################################
library(methylKit)
library(genomation)
source("meth_analysis_funcs.R")

#
#1. Read methylation data from bismark output
#

# Bismark CpG coverage files. 
# CpG methylations from both strands are combined in those files. Example lines:
# chr1    69      70      0.000000        0       3
# chr1    131     132     0.000000        0       5
file.list <- list("EH1516C.merged_CpG_evidence.cov",
                  "EH1516D.merged_CpG_evidence.cov",
                  "EH217A.merged_CpG_evidence.cov",
                  "EH217B.merged_CpG_evidence.cov",
                  "EH217C.merged_CpG_evidence.cov")

# Read methylation data
methobj <- methRead(file.list, sample.id = as.list(sample.ids),
                    assembly = "ehux", pipeline = "bismarkCoverage", treatment = c(0, 0, 1, 1, 1),
                    context = "CpG", mincov = 3, header = FALSE)

# Read Ehux gene structures
gene.parts <- readTranscriptFeatures("../v2/Ehux_genbank.bed", remove.unusual = FALSE, unique.prom = FALSE)

m <- getFeatureMethyl(methobj, gene.parts$promoters, sample.ids)


# Save the data in the workspace
save(list = ls(all.names = TRUE), file = "CpG_analysis.RData", envir = )