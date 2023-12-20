library(tidyverse)
# A DMC is defined as a off site increase to methylation level > 90%
# or an on site decreases to methylation level < 10%

# meth has a minimum coverage of 3 for each sample
meth$Ehx1516 <- (meth$numCs1+meth$numCs2) / (meth$coverage1+meth$coverage2)
meth$Ehx217 <- (meth$numCs3 + meth$numCs4 + meth$numCs5) / 
  (meth$coverage3 + meth$coverage4 + meth$coverage5)
meth.dmc <- getData(meth) %>% filter(Ehx1516 <= 0.1 & Ehx217 >= 0.9)
meth.dmc <- rbind(meth.dmc, getData(meth) %>% filter(Ehx1516 >= 0.9 & Ehx217 <= 0.1))
meth.dmc <- merge(meth.dmc, getData(myDiff.DSS), by=c("chr", "start", "end"), sort=FALSE)

#up_down_2000 <- GRangesList(up2000, down2000)
#names(up_down_2000) <- c("up2000", "down2000")
# Find percentage of dmc/r with genomic feature regions
# dm: dmc/dmr
#overlapFeatureRegions <- function(dm)

# Plot DMCs annotated by gene regions
plot_DMCs_by_Regions <- function(DMC.gr, filename) {
  dmc.genes <- annotateWithFeature(DMC.gr, genes$Gene, feature.name = "gene regions")
  
  dmc.exons <- annotateWithFeature(DMC.gr, gene.parts$exons, feature.name = "exons")
  
  dmc.introns <- annotateWithFeature(DMC.gr, gene.parts$introns, feature.name = "introns")
  
  dmc.up2000 <- annotateWithFeature(DMC.gr, up2000, feature.name = "up2000")
  dmc.down2000 <- annotateWithFeature(DMC.gr, down2000, feature.name = "down2000")
  dmc.promoters <- annotateWithFeature(DMC.gr, gene.parts$promoters, 
                                             feature.name = "promoters")
  
  png(filename = filename, width=6600, height=3000, res=600)
  par(mfrow=c(2, 3))
  genomation::plotTargetAnnotation(dmc.genes, col=c("green","gray","white"),
                                   main="overlapped with gene regions")
  genomation::plotTargetAnnotation(dmc.exons, col=c("green","gray","white"),
                                   main="overlapped with exons")
  genomation::plotTargetAnnotation(dmc.introns, col=c("green","gray","white"),
                                   main="overlapped with introns")
  genomation::plotTargetAnnotation(dmc.up2000, col=c("green","gray","white"),
                                   main="overlapped with 2000bp upstream from TSS")
  genomation::plotTargetAnnotation(dmc.down2000, col=c("green","gray","white"),
                                   main="overlapped with 2000bp downstream from TSS")
  genomation::plotTargetAnnotation(dmc.promoters, col=c("green","gray","white"),
                                   main="overlapped with promoters\n(1000bp up and downstream around TSS)")
  dev.off()
}

plot_DMCs_by_Regions(moabs.DMC, "moabs_graphs/moabs_dmc_by_regions.png")
plot_DMCs_by_Regions(as(dmls, "GRanges"), "dss_dmc_by_regions.png")

# Get the DMC counts in different regions in DEGs vs non-DEGs, and up-DEGs vs down-DEGs
# The DMCs input can come from different programs, e.g. MOABS, DSS, and on-n-off
DMC_count_in_region <- function(DMCs, region.gr, min_len=100) {
  dmc.count = as.data.frame(region.gr)
  dmc.count$len = dmc.count$end - dmc.count$start + 1
  dmc.gr = as(DMCs, "GRanges")
  dmc.count$dmcs = countOverlaps(region.gr, dmc.gr, type = "any", ignore.strand=TRUE)
  dmc.hypo = subset(dmc.gr, meth.diff < 0)
  dmc.hyper = subset(dmc.gr, meth.diff > 0)
  dmc.count$hypo = countOverlaps(region.gr, dmc.hypo, type = "any", ignore.strand=TRUE)
  dmc.count$hyper = countOverlaps(region.gr, dmc.hyper, type = "any", ignore.strand=TRUE)
  
  # filter out regions shorter than min_len
  dmc.count = dmc.count %>% filter(len >= min_len)
  
  return(dmc.count)
}

# Two sample test of DEGs vs. non-DEGs
DEG_nonDEG_test <- function(dmc.count) {
  dmc.count.DEG = dmc.count %>% filter (ID %in% rownames(annot.DE))
  dmc.count.nonDEG = dmc.count %>% filter (!(ID %in% rownames(annot.DE)))
  
  x = 1000 * dmc.count.DEG$dmcs / dmc.count.DEG$len
  y = 1000 * dmc.count.nonDEG$dmcs / dmc.count.nonDEG$len
  return(t.test(x, y))
  
  # tests = list()
  # results = list()
  # # two sample t test of on/off DMCs
  # x = 1000*dmc.count.DEG$onoff_dmc / dmc.count.DEG$len
  # y = 1000*dmc.count.nonDEG$onoff_dmc / dmc.count.nonDEG$len
  # tests <- append(tests, "on/off DMCs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # # two sample t test of MOABS CpG DMCs
  # x = 1000*dmc.count.DEG$moabs / dmc.count.DEG$len
  # y = 1000*dmc.count.nonDEG$moabs / dmc.count.nonDEG$len
  # tests <- append(tests, "MOABS CpG DMCs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # # two sample t test of MOABS CpG DMRs
  # x = dmc.count.DEG$dmrs
  # y = dmc.count.nonDEG$dmrs
  # tests <- append(tests, "MOABS CpG DMRs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # # two sample t test of MOABS CHG DMCs
  # x = 1000*dmc.count.DEG$hg.dmcs / dmc.count.DEG$len
  # y = 1000*dmc.count.nonDEG$hg.dmcs / dmc.count.nonDEG$len
  # tests <- append(tests, "MOABS CHG DMCs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # # two sample t test of MOABS CHG DMRs
  # x = 1000*dmc.count.DEG$hg.dmrs / dmc.count.DEG$len
  # y = 1000*dmc.count.nonDEG$hg.dmrs / dmc.count.nonDEG$len
  # tests <- append(tests, "MOABS CHG DMRs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # return(list("tests"=tests, "results"=results))
}


# Two sample test of up-regulated vs. down-regulated DEGs
DEG_up_down_test <- function(dmc.count) {
  # select promoters.dmc.count subset in DEGs
  dmc.count.DEG.up = dmc.count %>% filter (ID %in% rownames(annot.DE.up))
  dmc.count.DEG.down = dmc.count %>% filter (ID %in% rownames(annot.DE.down))
  
  x = 1000*dmc.count.DEG.up$dmcs / dmc.count.DEG.up$len
  y = 1000*dmc.count.DEG.down$dmcs / dmc.count.DEG.down$len
  return(t.test(x, y))
  
  # tests = list()
  # results = list()
  # 
  # # Two sample tests of on-off DMCs
  # x = 1000*dmc.count.DEG.up$onoff_dmc / dmc.count.DEG.up$len
  # y = 1000*dmc.count.DEG.down$onoff_dmc / dmc.count.DEG.down$len
  # tests <- append(tests, "On-off DMCs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # # Two sample tests between DE.up and DE.down
  # x = 1000*dmc.count.DEG.up$moabs / dmc.count.DEG.up$len
  # y = 1000*dmc.count.DEG.down$moabs / dmc.count.DEG.down$len
  # tests <- append(tests, "MOABS CpG DMCs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # # Two sample test of hyper and hypo DMCs in up and down DEGs
  # x = 1000*dmc.count.DEG.up$dmc.hypo / dmc.count.DEG.up$len
  # y = 1000*dmc.count.DEG.down$dmc.hypo / dmc.count.DEG.down$len
  # tests <- append(tests, "hypo CpG DMCs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # x = 1000*dmc.count.DEG.up$dmc.hyper / dmc.count.DEG.up$len
  # y = 1000*dmc.count.DEG.down$dmc.hyper / dmc.count.DEG.down$len
  # tests <- append(tests, "hypo CpG DMCs")
  # results[length(results)+1] <- list(t.test(x, y))
  # 
  # return(list("tests"=tests, "results"=results))
  
}

# Test for promoter regions
promoters.dmc.count = DMC_count_in_region(moabs.DMC, promoters.gr)
DEG_nonDEG_test(promoters.dmc.count)
DEG_up_down_test(promoters.dmc.count)

promoters.dmc.count$dmrs <- countOverlaps(as(promoters.dmc.count, "GRanges"), moabs.DMR, type="any", ignore.strand=TRUE)
genes.dmc.count$dmrs <- countOverlaps(as(genes.dmc.count, "GRanges"), moabs.DMR, type="any", ignore.strand=TRUE)
# find the number of DEGs and non-DEGs that overlap with DMRs
dmc.count.DEG = genes.dmc.count %>% filter (ID %in% rownames(annot.DE))
dmc.count.nonDEG = genes.dmc.count %>% filter (!(ID %in% rownames(annot.DE)))
sum(dmc.count.DEG$dmrs > 0) / nrow(dmc.count.DEG)
sum(dmc.count.DEG$dmr ) / sum(dmc.count.DEG$dmrs > 0)
sum(dmc.count.nonDEG$dmrs > 0) / nrow(dmc.count.nonDEG)
sum(dmc.count.nonDEG$dmr ) / sum(dmc.count.nonDEG$dmr > 0)

dmc.count.DEG = promoters.dmc.count %>% filter (ID %in% rownames(annot.DE))
dmc.count.nonDEG = promoters.dmc.count %>% filter (!(ID %in% rownames(annot.DE)))
sum(dmc.count.DEG$dmrs > 0) / nrow(dmc.count.DEG)
sum(dmc.count.nonDEG$dmrs > 0) / nrow(dmc.count.nonDEG)

# CHG context
promoters.dmc.chg.count = DMC_count_in_region(as(dmls.CHG, "GRanges"), promoters.gr) 
DEG_nonDEG_test(promoters.dmc.chg.count)
DEG_up_down_test(promoters.dmc.chg.count)

# separate hyper and hypo-methylated DMCs
moabs.DMC.hypo = subset(moabs.DMC, meth.diff < 0)
moabs.DMC.hyper = subset(moabs.DMC, meth.diff > 0)
moabs.HG.DMC.hypo = subset(moabs.HG.DMC, meth.diff < 0)
moabs.HG.DMC.hyper = subset(moabs.HG.DMC, meth.diff > 0)


# Test for gene regions
genes.dmc.count = DMC_count_in_region(moabs.DMC, as(gene_regions, "GRanges"))
DEG_nonDEG_test(genes.dmc.count)
DEG_up_down_test(genes.dmc.count)
genes.dmc.count$DEG = genes.dmc.count$ID %in% rownames(annot.DE)
genes.dmc.count$density = 1000 * genes.dmc.count$dmcs / genes.dmc.count$len
# plot a histogram
x = genes.dmc.count %>% filter(ID %in% rownames(annot.DE))
y = genes.dmc.count %>% filter(!(ID %in% rownames(annot.DE)))
ggplot(genes.dmc.count, aes(x=density, color=DEG)) +
  geom_histogram(fill="white", position="dodge")

# Consider genes with at least 3 DMCs in the gene regions
genes.DE.3dmcs = genes.dmc.count %>% filter(ID %in% rownames(annot.DE) & dmcs >= 3)
ggplot(genes.DE.3dmcs, aes(x=density, color=DEG)) +
  geom_histogram(fill="white", position="dodge", binwidth = 1)
genes.nonDE.3dmcs = genes.dmc.count %>% filter(!(ID %in% rownames(annot.DE)) & dmcs >= 3)
ggplot(genes.nonDE.3dmcs, aes(x=density, color=DEG)) +
  geom_histogram(fill="white", position="dodge", binwidth = 1)
# calculate average methylation changes of DMCs for genes with at least 3DMCs
x = findOverlaps(as(genes.DE.3dmcs, "GRanges"), moabs.DMC, type = "any")
x = as.data.frame(x)
genes.DE.3dmcs$dmc_mean = 0
for(i in 1:nrow(genes.DE.3dmcs)) {
  genes.DE.3dmcs[i, "dmc_mean"] = mean(moabs.DMC[x[x$queryHits == i, ]$subjectHits]$meth.diff)
}
x = as.data.frame(findOverlaps(as(genes.nonDE.3dmcs, "GRanges"), moabs.DMC, type = "any"))
genes.nonDE.3dmcs$dmc_mean = 0
for(i in 1:nrow(genes.nonDE.3dmcs)) {
  genes.nonDE.3dmcs[i, "dmc_mean"] = mean(moabs.DMC[x[x$queryHits == i, ]$subjectHits]$meth.diff)
}

genes.dmc.count = DMC_count_in_region(as(dmls, "GRanges"), as(gene_regions, "GRanges"))
DEG_nonDEG_test(genes.dmc.count)

r = DEG_nonDEG_test(genes.dmc.count)
r$results
r = DEG_up_down_test(genes.dmc.count)
r$results

# Test for up2000 regions
#up2000.dmc.count = DM_count_by_region(up2000, min_len = 500)
up2000.dmc.count = DMC_count_in_region(moabs.DMC, up2000)
r = DEG_nonDEG_test(up2000.dmc.count)
r$results

# Test for down2000 regions
#down2000.dmc.count = DM_count_by_region(down2000, min_len = 500)
down2000.dmc.count = DMC_count_in_region(moabs.DMC, down2000)
r = DEG_nonDEG_test(down2000.dmc.count)
r$results

DE_nonDEG_conf_interval_by_regions <- function(dmc.count, region="promoters", sel="moabs", conf.int=0.95) {
  dmc.count.DEG = dmc.count %>% filter (ID %in% rownames(annot.DE))
  dmc.count.nonDEG = dmc.count %>% filter (!(ID %in% rownames(annot.DE)))
  data.DEG = dmc.count.DEG[, sel] / dmc.count.DEG$len * 1000
  data.nonDEG = dmc.count.nonDEG[, sel] / dmc.count.nonDEG$len * 1000
  
  d = data.frame("subset"="DEGs", "region" = region, "mean" = mean(data.DEG), 
                 "lci" = t.test(data.DEG, conf.level=conf.int)$conf.int[1], 
                 "uci" = t.test(data.DEG, conf.level=conf.int)$conf.int[2]
  )
  
  d1 = data.frame("subset"="non-DEGs", "region" = region, "mean" = mean(data.nonDEG), 
                  "lci" = t.test(data.nonDEG, conf.level=conf.int)$conf.int[1], 
                  "uci" = t.test(data.nonDEG, conf.level=conf.int)$conf.int[2]
  )
  d = rbind(d, d1)
  return(d)
}

DE_up_down_conf_interval_by_regions <- function(dmc.count, region="promoters", sel="moabs", conf.int=0.95) {
  dmc.count.DEG = dmc.count %>% filter (ID %in% rownames(annot.DE))
  dmc.count.DEG.up = dmc.count %>% filter (ID %in% rownames(annot.DE.up))
  dmc.count.DEG.down = dmc.count %>% filter (ID %in% rownames(annot.DE.down))
  dmc.count.nonDEG = dmc.count %>% filter (!(ID %in% rownames(annot.DE)))
  data.DEG = dmc.count.DEG[, sel] / dmc.count.DEG$len * 1000
  data.DEG.up = dmc.count.DEG.up[, sel] / dmc.count.DEG.up$len * 1000
  data.DEG.down = dmc.count.DEG.down[, sel] / dmc.count.DEG.down$len * 1000
  data.nonDEG = dmc.count.nonDEG[, sel] / dmc.count.nonDEG$len * 1000
  
  d = data.frame("subset"="up DEGs", "region" = region, "mean" = mean(data.DEG.up), 
                 "lci" = t.test(data.DEG.up, conf.level=conf.int)$conf.int[1], 
                 "uci" = t.test(data.DEG.up, conf.level=conf.int)$conf.int[2]
  )
  d1 = data.frame("subset"="down DEGs", "region" = region, "mean" = mean(data.DEG.down), 
                 "lci" = t.test(data.DEG.down, conf.level=conf.int)$conf.int[1], 
                 "uci" = t.test(data.DEG.down, conf.level=conf.int)$conf.int[2]
  )
  d2 = data.frame("subset"="non-DEGs", "region" = region, "mean" = mean(data.nonDEG), 
                  "lci" = t.test(data.nonDEG, conf.level=conf.int)$conf.int[1], 
                  "uci" = t.test(data.nonDEG, conf.level=conf.int)$conf.int[2]
  )
  d = rbind(d, d1, d2)
  return(d)
}


CpG_DMC_conf <- DE_nonDEG_conf_interval_by_regions(promoters.dmc.count, region="promoters", conf.int = 0.99)
CpG_DMC_conf <- rbind(CpG_DMC_conf, 
                      DE_nonDEG_conf_interval_by_regions(genes.dmc.count, region="gene regions", conf.int = 0.99))
CpG_DMC_conf <- rbind(CpG_DMC_conf, 
                      DE_nonDEG_conf_interval_by_regions(up2000.dmc.count, region = "up2000", conf.int = 0.99))
CpG_DMC_conf <- rbind(CpG_DMC_conf, 
                      DE_nonDEG_conf_interval_by_regions(down2000.dmc.count, region = "down2000", conf.int = 0.99))

CpG_DMC_conf <- DE_up_down_conf_interval_by_regions(promoters.dmc.count, region="promoters", conf.int = 0.99)
CpG_DMC_conf <- rbind(CpG_DMC_conf, 
                      DE_up_down_conf_interval_by_regions(genes.dmc.count, region="gene regions", conf.int = 0.99))
CpG_DMC_conf <- rbind(CpG_DMC_conf, 
                      DE_up_down_conf_interval_by_regions(up2000.dmc.count, region = "up2000", conf.int = 0.99))
CpG_DMC_conf <- rbind(CpG_DMC_conf, 
                      DE_up_down_conf_interval_by_regions(down2000.dmc.count, region = "down2000", conf.int = 0.99))
png(filename = "region_analysis/DMC_DE_up_down_Region_conf0.99.png", 
    width=4800, height=3600, res=600)
ggplot(CpG_DMC_conf) + geom_point(aes(x=region, y=mean, color=subset), position = position_dodge(width=0.2)) +
  geom_errorbar(aes(x=region, ymin=lci, ymax=uci, color=subset), position = "dodge", width=0.2) +
  ggtitle("DEG vs non-DEG CpG DMC count per 1000bp by Regions") +
  labs(x = "Region", y = "DMC count per 1000bp") +
  theme_light() + #scale_y_continuous(limits=c(10, 20), n.breaks = 10) +
  theme(axis.title=element_text(size=12), plot.title=element_text(size=14, hjust=0.5)) +
  theme(legend.position = c(0.8, 0.8), legend.direction = "vertical", legend.text = element_text(size=10))
dev.off()

CHG_DMC_conf <- DE_nonDEG_conf_interval_by_regions(promoters.dmc.count, region="promoters", 
                                                   sel = "hg.dmcs", conf.int = 0.99)
CHG_DMC_conf <- rbind(CHG_DMC_conf, 
                      DE_nonDEG_conf_interval_by_regions(genes.dmc.count, region="gene regions", 
                                                         sel = "hg.dmcs", conf.int = 0.99))
CHG_DMC_conf <- rbind(CHG_DMC_conf, 
                      DE_nonDEG_conf_interval_by_regions(up2000.dmc.count, region = "up2000", 
                                                         sel = "hg.dmcs", conf.int = 0.99))
CHG_DMC_conf <- rbind(CHG_DMC_conf, 
                      DE_nonDEG_conf_interval_by_regions(down2000.dmc.count, region = "down2000", 
                                                         sel = "hg.dmcs", conf.int = 0.99))

png(filename = "region_analysis/CHG_Region_conf0.99.png", 
    width=4800, height=3600, res=600)
ggplot(CHG_DMC_conf) + geom_point(aes(x=region, y=mean, color=subset), position = position_dodge(width=0.2)) +
  geom_errorbar(aes(x=region, ymin=lci, ymax=uci, color=subset), position = "dodge", width=0.2) +
  ggtitle("DEG vs non-DEG CHG DMC count per 1000bp by Regions") +
  labs(x = "Region", y = "DMC count per 1000bp") +
  theme_light() + #scale_y_continuous(limits=c(10, 20), n.breaks = 10) +
  theme(axis.title=element_text(size=12), plot.title=element_text(size=14, hjust=0.5)) +
  theme(legend.position = c(0.7, 0.9), legend.direction = "horizontal", legend.text = element_text(size=10))
dev.off()

# Calculate the number of DMCs per 1000 bp in different regions: genes, exons, introns, promoters, up2000, down2000
calc_DMC_region_rate <- function(dmcs) {
  gene_len_perc <- sum(width(GenomicRanges::reduce(genes$Gene))) / genome_len *100
  exons_len_perc <- sum(width(GenomicRanges::reduce(gene.parts$exons))) / genome_len * 100
  introns_len_perc <- sum(width(GenomicRanges::reduce(gene.parts$introns))) / genome_len * 100
  promoters_len_perc <- sum(width(GenomicRanges::reduce(
    GenomicRanges::intersect(gene.parts$promoters, scaffolds.gr, ignore.strand=TRUE)))) / genome_len * 100
  up2000_len_perc <- sum(width(GenomicRanges::reduce(
    GenomicRanges::intersect(up2000, scaffolds.gr, ignore.strand=TRUE)))) / genome_len * 100
  down2000_len_perc <- sum(width(GenomicRanges::reduce(
    GenomicRanges::intersect(down2000, scaffolds.gr, ignore.strand=TRUE)))) / genome_len * 100
  
  dmc.genes <- annotateWithFeature(dmcs, genes$Gene, feature.name = "gene regions")
  dmc.exons <- annotateWithFeature(dmcs, gene.parts$exons, feature.name = "exons")
  dmc.introns <- annotateWithFeature(dmcs, gene.parts$introns, feature.name = "introns")
  dmc.up2000 <- annotateWithFeature(dmcs, up2000, feature.name = "up2000")
  dmc.down2000 <- annotateWithFeature(dmcs, down2000, feature.name = "down2000")
  dmc.promoters <- annotateWithFeature(dmcs, gene.parts$promoters, feature.name = "promoters")
  
  return ( (length(dmcs) * 1000 / genome_len) * c( 
    slot(dmc.up2000, "annotation")[1] / up2000_len_perc,
    slot(dmc.promoters, "annotation")[1] / promoters_len_perc,
    slot(dmc.genes, "annotation")[1]  / gene_len_perc ,
    slot(dmc.exons, "annotation")[1]  / exons_len_perc,
    slot(dmc.introns, "annotation")[1] / introns_len_perc,
    slot(dmc.down2000, "annotation")[1] / down2000_len_perc))
}

# CpG context

# Relative ratio of DMCs in features vs. feature sizes
dmc.regions.rat <- data.frame(region=c("up2000", "promoters", "genes", "exons", "introns", "down2000"), 
                              CpG=0, CHG=0)
rownames(dmc.regions.rat) <- dmc.regions.rat$region
# calc the DMC rate per 1000 bp
dmc.regions.rat$CpG = calc_DMC_region_rate(moabs.DMC)

moabs.dmc.Ann <- annotateWithGeneParts(moabs.DMC, gene.parts)
genomation::plotTargetAnnotation(moabs.dmc.Ann)

# CHG context

#moabs.HG.Ann <- annotateWithGeneParts(moabs.HG.DMC, gene.parts)

dmc.regions.rat$CHG <- calc_DMC_region_rate(moabs.HG.DMC)
  
rownames(dmc.regions.rat) <- seq(1,6)
region_order <- dmc.regions.rat$region
png(filename = "moabs_graphs/DMC_region_rate.png", width=6600, 
    height=3000, res=600)
ggplot(dmc.regions.rat) + geom_point(aes(x = factor(region, level = region_order), y=CpG, colour="CpG")) +
  geom_line(aes(x = seq(1,6), y=CpG, colour="CpG")) +
  geom_point(aes(x = factor(region, level = region_order), y=CHG, colour="CHG")) +
  geom_line(aes(x = seq(1,6), y=CHG, colour="CHG")) +
  scale_colour_manual("", breaks = c("CpG", "CHG"), values = c("CpG"= "red", "CHG" = "blue")) +
  #geom_hline(yintercept=1.0, linetype="dashed") +
  theme_light() +
  scale_y_continuous(limits=c(0.5, 2.5), n.breaks = 7) +
  labs(x = "Feature Region", y = "DMC Rate per 1000 bp") +
  theme(axis.title=element_text(size=12), plot.title=element_text(size=12, hjust=0.5)) +
  theme(legend.position = c(0.8, 0.85), legend.direction = "horizontal", legend.text = element_text(size=12))
dev.off()
#ggsave("moabs_graphs/DMC_region_rate.png", dpi=600)

png(filename = "moabs_graphs/moabs_HG_by_regions.png", width=6600, height=3000, res=600)
par(mfrow=c(2, 3))
genomation::plotTargetAnnotation(moabs.HG.genes, col=c("blue","gray","white"),
                                 main="overlapped with gene regions")
genomation::plotTargetAnnotation(moabs.HG.exons, col=c("blue","gray","white"),
                                 main="overlapped with exons")
genomation::plotTargetAnnotation(moabs.HG.introns, col=c("blue","gray","white"),
                                 main="overlapped with introns")
genomation::plotTargetAnnotation(moabs.HG.up2000, col=c("blue","gray","white"),
                                 main="overlapped with 2000bp upstream from TSS")
genomation::plotTargetAnnotation(moabs.HG.down2000, col=c("blue","gray","white"),
                                 main="overlapped with 2000bp downstream from TSS")
genomation::plotTargetAnnotation(moabs.HG.promoters, col=c("blue","gray","white"),
                                 main="overlapped with promoters\n(1000bp up and downstream around TSS)")
dev.off()

# Plot the ratio of  DMC numbers / scaffold lengths
scaffolds_100$CpG_DMC <- 0
moabs.DMC.df = as.data.frame(moabs.DMC)
colnames(moabs.DMC.df)[1] = "chr"
for (i in 1:nrow(scaffolds_100)) {
  nDMC <-  moabs.DMC.df %>% filter(chr == scaffolds_100[i, ]$chr) %>% nrow()
  scaffolds_100[i, ]$CpG_DMC <- nDMC
}
scaffolds_100$CHG_DMC <- 0
moabs.HG.DMC.df = as.data.frame(moabs.HG.DMC)
colnames(moabs.HG.DMC.df)[1] = "chr"
for (i in 1:nrow(scaffolds_100)) {
  nDMC <- moabs.HG.DMC.df %>% filter(chr == scaffolds_100[i, ]$chr) %>% nrow()
  scaffolds_100[i, ]$CHG_DMC <- nDMC
}

png(filename = "moabs_graphs/DMC_Scaffold_Ratio.png", width=6600, height=3000, res=600)
ggplot(scaffolds_100) + geom_line(mapping=aes(x=seq(1, 100), y = CpG_DMC/len*1000, colour="CpG")) +
  geom_line(mapping=aes(x=seq(1, 100), y = CHG_DMC/len*1000, colour="CHG")) +
  geom_hline(aes(yintercept=mean(CpG_DMC/len)*1000), color="red", linetype="dashed") +
  scale_colour_manual("", breaks = c("CpG", "CHG"), values = c("CpG"= "red", "CHG" = "blue")) + 
  geom_hline(aes(yintercept=mean(CHG_DMC/len)*1000), color="blue", linetype="dashed") +
  theme_light() + #  ggtitle("Numbers of DMCs per 1000 bp by Scaffolds") +
  labs(x = "Scaffolds", y = "Avg Num of DMCs per 1000 bp") +
  scale_y_continuous(limits=c(0, 9), n.breaks = 7) +
  theme(axis.title=element_text(size=12), plot.title=element_text(size=14, hjust=0.5)) +
  theme(legend.position = c(0.8, 0.85), legend.direction = "horizontal", legend.text = element_text(size=12))
dev.off()
#ggsave("moabs_graphs/DMC_Scaffold_Ratio.png", dpi=600)

# Statistics of DMCs in the first 100 scaffolds
x = scaffolds_100$CpG_DMC / scaffolds_100$len * 1000
summary(x)
y = scaffolds_100$CHG_DMC / scaffolds_100$len * 1000
summary(y)
cor(x, y)
cor.test(x, y)

