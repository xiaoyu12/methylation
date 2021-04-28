###################################################################
# This script merges CpG methylation data of the same E. hux sample
###################################################################

# Read CpG evidence files
EH1516C <- read.table("EH1516C.merged_CpG_evidence.cov", sep = "\t")
EH1516D <- read.table("EH1516D.merged_CpG_evidence.cov", sep = "\t")
EH1516 <- merge(EH1516C, EH1516D, by=c("V1", "V2", "V3"), all=FALSE, sort=FALSE)

EH217A <- read.table("EH217A.merged_CpG_evidence.cov", sep="\t")
EH217B <- read.table("EH217B.merged_CpG_evidence.cov", sep="\t")
EH217C <- read.table("EH217C.merged_CpG_evidence.cov", sep="\t")
EH217AB <- merge(EH217A, EH217B, by=c("V1", "V2", "V3"), all=FALSE, sort=FALSE)
EH217 <- merge(EH217AB, EH217C, by=c("V1", "V2", "V3"), all=FALSE, sort=FALSE)

# Sum up methylated and total C sites to calcaute the average methylation rate
EH1516$methc <- EH1516$V5.x + EH1516$V5.y
EH1516$totc <- EH1516$V5.x +  EH1516$V6.x + EH1516$V5.y + EH1516$V6.y
EH1516$beta <- EH1516$methc / EH1516$totc

EH217$methc <- EH217$V5 + EH217$V5.x + EH217$V5.y
EH217$totc <- EH217$V5 + EH217$V5.x + EH217$V5.y + EH217$V6 + EH217$V6.x + EH217$V6.y
EH217$beta <- EH217$methc / EH217$totc

#write.table(EH1516[, c("V1", "V2", "V3", "methc", "totc", "beta")], 
#            file="EH1516_meth_merged.tsv",
#            sep="\t", row.names = FALSE, quote = FALSE)

#write.table(EH217[, c("V1", "V2", "V3", "methc", "totc", "beta")],
#            file="EH217_meth_merged.tsv",
#            sep="\t", row.names = FALSE, quote = FALSE)


# merge EH1516 and EH217 methylation data 
EH1516_tot = EH1516[, c("V1", "V2", "V3", "methc", "totc", "beta")]
EH217_tot = EH217[, c("V1", "V2", "V3", "methc", "totc", "beta")]
EH_all = merge(EH1516_tot, EH217_tot, by=c("V1", "V2", "V3"), all=FALSE, sort=FALSE)
colnames(EH_all)[4:9] = c("methc.1516", "totc.1516", "beta.1516", 
                          "methc.217", "totc.217", "beta.217")

write.table(EH_all, file="Ehux_meth_merged.tsv", sep="\t", row.names=FALSE, quote=FALSE)



