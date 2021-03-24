library(limma)
source("meth_analysis_funcs.R")
# Density plot of methylation in promtoers
prom.m <- methyl_to_data_frame(promobj, promoters)
colnames(prom.m) <- sample.ids

colramp = rainbow(10)
# plot density of methy M-values
plot(density(prom.m[, 1]), col=colramp[1], ylim=c(0, 0.8),
     xlab = 'Methylation M-values')
for(i in 2:5) {
  lines(density(prom.m[, i]), col=colramp[i])
}

legend("topright", inset=0.05, sample.ids, fill=colramp, cex=0.6)

plotMDS(prom.m, top=1000, col=colramp, dim=c(1,2))

# Differentially analysis using limma
#conditions = data.frame(conditions=factor(c(rep("EH1516", 2), rep("EH217", 3))))
#rownames(conditions) <- colnames(prom.m)
conditions <- factor(c(rep("EH1516", 2), rep("EH217", 3)))
design <- model.matrix(~conditions)
fit <- lmFit(prom.m, design)
fit <- eBayes(fit)
summary(decideTests(fit))
t <- topTable(fit, number = Inf, coef=2)
prom.dm <- subset(t, adj.P.Val < 0.05) # Differential methylated promoter regions

# Calculate methylations of CPG islands and shores
#cpgi.Ann <- annotateWithGeneParts(cpgi[[1]], gene.parts)
#cpgi.TSS <- getAssociationWithTSS(cpgi.Ann)
cpg_i <- cpgi$CpGi
cpg_i$name = gene.parts$TSSes$name    # Add a gene name column
cpg_i$ID = mapping[cpg_i$name, "ID"]

m.cpg_i <- getFeatureMethyl(methobj, cpg_i, sample.ids)
fit <- lmFit(m.cpg_i$m, design)
fit <- eBayes(fit)
summary(decideTests(fit)) 
t <- topTable(fit, number = Inf, coef=2)
cpg_i.dm <- subset(t, adj.P.Val < 0.05)

#cpg_s <- as.data.frame(cpgi[[2]])
