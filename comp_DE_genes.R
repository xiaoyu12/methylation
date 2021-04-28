# Compare DE genes in 217 vs. 1516 against DE genes in 1516 9S vs. 9

# DE genes in strain 217 vs. 1516
head(annot.DE)

annot.9v9s <- read.csv("../../gene_expression/ehux/output_genbank/DE_9-217_vs_9S-217_0.1.annot.txt", 
                        header = TRUE, row.names = 1, sep="\t")
#annot.9v9s <- read.csv("../../gene_expression/ehux/output_genbank/DE_9-217_vs_9S-217_0.05.annot.txt", 
#                       header = TRUE, row.names = 1, sep="\t")
annot.9v9s.down <- subset(annot.9v9s, log2FC < 0)
annot.9v9s.up <- subset(annot.9v9s, log2FC > 0)

shared_ids <- intersect(rownames(annot.9v9s), rownames(annot.DE))
shared <- annot.9v9s[shared_ids, ]
shared  <- cbind(annot.DE[shared_ids, ]$log2FC, shared)
colnames(shared)[1] <- "log2FC_217vs1516"
cor(shared$log2FC_217vs1516, shared$log2FC)

plot(shared$log2FC_217vs1516, shared$log2FC, type="p", pch=20, cex=0.5, col="blue", 
     xlab = "LFC strain 217 vs. 1516", ylab = 'LFC 9S vs. 9 ')
abline(h=0, col='black')
abline(v=0, col='black')
t <- subset(shared, (log2FC <0 & log2FC_217vs1516 >0) |  (log2FC >0 & log2FC_217vs1516 <0))
points(t$log2FC_217vs1516, t$log2FC, col='red', pch=20, cex=0.5)

# Select genes in quartet 1
shared.up <- subset(shared, log2FC >0 & log2FC_217vs1516 >0)
shared.down <- subset(shared,log2FC < 0 & log2FC_217vs1516 < 0)
write.csv(shared.up, file="DE_shared_up.csv")

