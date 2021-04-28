# Read expression matrix and calculate expression changes using edgeR
library(limma)
library(edgeR)
library(DESeq2)
mat <- read.table("EHX_1516_v_217_Jan2019.txt", row.names = 1, header = TRUE, sep = "\t")
conditions <- factor(c(rep("EH1516", 3), rep("EH217", 3)))
design <- model.matrix(~conditions)
DGE <- DGEList(mat)
DGE <- calcNormFactors(DGE)
v <- voom(DGE, design = design, plot = TRUE)
fit <- eBayes(lmFit(v, design))
expr.coef <- topTable(fit, number = Inf, coef = 2, sort.by = "P")

#DESeq2 analysis
conditions = data.frame(conditions)
ddsTable <- DESeqDataSetFromMatrix(countData = mat, colData = conditions, design = ~conditions)
ddsTable <- estimateSizeFactors(ddsTable)
dds <- DESeq(ddsTable)
res <- results(dds)
res$padj[is.na(res$padj)]  <- 1
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
expr.coef <- cbind(res[rownames(expr.coef), ]$log2FoldChange, expr.coef)
colnames(expr.coef)[1] <- "LFC" # logFC from DESeq2 analysis
expr.coef <- cbind(res[rownames(expr.coef), ]$baseMean, expr.coef)
colnames(expr.coef)[1] <- "baseMean"
