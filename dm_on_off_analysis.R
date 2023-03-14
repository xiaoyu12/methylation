# A DMC is defined as a off site increase to methylation level > 90%
# or an on site decreases to methylation level < 10%
#meth10.dmc <- data.frame(matrix(ncol=length(meth10), nrow = 0))
#colnames(meth10.dmc) <- colnames(meth10)
meth10$Ehx1516 <- (meth10$numCs1+meth10$numCs2) / (meth10$coverage1+meth10$coverage2)
meth10$Ehx217 <- (meth10$numCs3 + meth10$numCs4 + meth10$numCs5) / 
  (meth10$coverage3 + meth10$coverage4 + meth10$coverage5)
meth10.dmc <- getData(meth10) %>% filter(Ehx1516 <= 0.1 & Ehx217 >= 0.9)
meth10.dmc <- rbind(meth10.dmc, getData(meth10) %>% filter(Ehx1516 >= 0.9 & Ehx217 <= 0.1))

meth10.dmc <- merge(meth10.dmc, getData(myDiff10), by=c("chr", "start", "end"), sort=FALSE)
meth10.dmc.Ann <- annotateWithGeneParts(as(meth10.dmc, "GRanges"), gene.parts)
genomation::plotTargetAnnotation(meth10.dmc.Ann)

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
meth.dmc.genes <- annotateWithFeature(as(meth.dmc, "GRanges"), genes$Gene, feature.name = "gene regions")

meth.dmc.exons <- annotateWithFeature(as(meth.dmc, "GRanges"), gene.parts$exons, feature.name = "exons")

meth.dmc.introns <- annotateWithFeature(as(meth.dmc, "GRanges"), gene.parts$introns, feature.name = "introns")

meth.dmc.up2000 <- annotateWithFeature(as(meth.dmc, "GRanges"), up2000, feature.name = "up2000")
meth.dmc.down2000 <- annotateWithFeature(as(meth.dmc, "GRanges"), down2000, feature.name = "down2000")
meth.dmc.promoters <- annotateWithFeature(as(meth.dmc, "GRanges"), gene.parts$promoters, 
                                          feature.name = "promoters")

genome_len <- sum(scaffolds$len)
# the percentages of genome are in the gene regions
scaffolds.gr <- scaffolds[, 1:2]
scaffolds.gr[, 2] <-1 
scaffolds.gr[, 3] <- scaffolds$len
colnames(scaffolds.gr) <- c("chr", "start", "end")
scaffolds.gr <- as(scaffolds.gr, "GRanges")

promoters.gr = as(promoters, 'GRanges')
# Get the promoter regions overlapped with the genome
o = findOverlaps(promoters.gr, scaffolds.gr)
grl1 = split(promoters.gr[queryHits(o)], 1:length(o)) # You can't mendoapply on a GRanges object
grl2 = split(scaffolds.gr[subjectHits(o)], 1:length(o))
# set the promoters min and max to the min and max of genome scaffold
foo = function(x, y) {
      rv = x
       start(rv) = max(start(x), start(y))
       end(rv) = min(end(x), end(y))
       return(rv)
}
promoters.gr = unlist(mendoapply(foo, grl1, y=grl2))

#promoters.gr = GenomicRanges::findOverlaps(promoters.gr, scaffolds.gr, ignore.strand=TRUE)
# count the number of DMCs in each promoter region
#promoters.dmc.count = countOverlaps(promoters.gr, meth.dmc.gr, type="any", ignore.strand=TRUE)
promoters.dmc.count = as.data.frame(promoters.gr)
n_dmc = countOverlaps(promoters.gr, meth.dmc.gr, type="any", ignore.strand=TRUE)
# count overlaps with moabs DMCs
moabs = countOverlaps(promoters.gr, moabs.DMC, type = "any", ignore.strand=TRUE)
#
promoters.dmc.count = cbind(promoters.dmc.count, n_dmc)
promoters.dmc.count = cbind(promoters.dmc.count, moabs)
promoters.dmc.count$len = promoters.dmc.count$end - promoters.dmc.count$start + 1
# Count promoter overlapped with moabs.DMRs
dmrs = countOverlaps(promoters.gr, moabs.DMR, type="any", ignore.strand=TRUE)
promoters.dmc.count = cbind(promoters.dmc.count, dmrs)
# Count promoter overlapped with moabs CHG DMCs
hg.dmcs = countOverlaps(promoters.gr, moabs.HG.DMC, type="any", ignore.strand=TRUE)
promoters.dmc.count$hg.dmcs = hg.dmcs
hg.dmrs = countOverlaps(promoters.gr, moabs.HG.DMR, type="any", ignore.strand=TRUE)
promoters.dmc.count$hg.dmrs = hg.dmrs

# select promoters.dmc.count subset in DEGs
#promoters.dmc.count = promoters.dmc.count %>% filter (len > 1000)
promoters.dmc.count.DEG = promoters.dmc.count %>% filter (ID %in% rownames(annot.DE))
promoters.dmc.count.nonDEG = promoters.dmc.count %>% filter (!(ID %in% rownames(annot.DE)))
promoters.dmc.count.DEG.up = promoters.dmc.count %>% filter (ID %in% rownames(annot.DE.up))
promoters.dmc.count.DEG.down = promoters.dmc.count %>% filter (ID %in% rownames(annot.DE.down))
# two sample t test of on/off DMCs
x = 1000*promoters.dmc.count.DEG$n_dmc / promoters.dmc.count.DEG$len
y = 1000*promoters.dmc.count.nonDEG$n_dmc / promoters.dmc.count.nonDEG$len
t.test(x, y)
# two sample t test of MOABS CpG DMCs
x = 1000*promoters.dmc.count.DEG$moabs / promoters.dmc.count.DEG$len
y = 1000*promoters.dmc.count.nonDEG$moabs / promoters.dmc.count.nonDEG$len
t.test(x, y)
# two sample t test of MOABS CpG DMRs
x = promoters.dmc.count.DEG$dmrs
y = promoters.dmc.count.nonDEG$dmrs
t.test(x,y)
# two sample t test of MOABS CHG DMCs
x = 1000*promoters.dmc.count.DEG$hg.dmcs / promoters.dmc.count.DEG$len
y = 1000*promoters.dmc.count.nonDEG$hg.dmcs / promoters.dmc.count.nonDEG$len
t.test(x, y)
# two sample t test of MOABS CHG DMRs
x = 1000*promoters.dmc.count.DEG$hg.dmrs / promoters.dmc.count.DEG$len
y = 1000*promoters.dmc.count.nonDEG$hg.dmrs / promoters.dmc.count.nonDEG$len
t.test(x, y)

# Two sample tests between DE.up and DE.down
x = 1000*promoters.dmc.count.DEG.up$moabs / promoters.dmc.count.DEG.up$len
y = 1000*promoters.dmc.count.DEG.down$moabs / promoters.dmc.count.DEG.down$len
t.test(x, y)

x = 1000*promoters.dmc.count.DEG.up$n_dmc / promoters.dmc.count.DEG.up$len
y = 1000*promoters.dmc.count.DEG.down$n_dmc / promoters.dmc.count.DEG.down$len
t.test(x, y)

gene_len_perc <- sum(width(GenomicRanges::reduce(genes$Gene))) / genome_len *100
exons_len_perc <- sum(width(GenomicRanges::reduce(gene.parts$exons))) / genome_len * 100
introns_len_perc <- sum(width(GenomicRanges::reduce(gene.parts$introns))) / genome_len * 100
promoters_len_perc <- sum(width(GenomicRanges::reduce(
  GenomicRanges::intersect(gene.parts$promoters, scaffolds.gr, ignore.strand=TRUE)))) / genome_len * 100
up2000_len_perc <- sum(width(GenomicRanges::reduce(
  GenomicRanges::intersect(up2000, scaffolds.gr, ignore.strand=TRUE)))) / genome_len * 100
down2000_len_perc <- sum(width(GenomicRanges::reduce(
  GenomicRanges::intersect(down2000, scaffolds.gr, ignore.strand=TRUE)))) / genome_len * 100

# Relative ratio of DMCs in features vs. feature sizes
dmc.regions.rat <- data.frame(region=c("genes", "exons", "introns", "promoters", "up2000", "down2000"), 
                              CpG=0, CHG=0)
rownames(dmc.regions.rat) <- dmc.regions.rat$region
dmc.regions.rat["genes", ]$CpG <- slot(meth.dmc.genes, "annotation")[1] / gene_len_perc
dmc.regions.rat["exons", ]$CpG<- slot(meth.dmc.exons, "annotation")[1] / exons_len_perc
dmc.regions.rat["introns", ]$CpG<- slot(meth.dmc.introns, "annotation")[1] / introns_len_perc
dmc.regions.rat["promoters", ]$CpG<- slot(meth.dmc.promoters, "annotation")[1] / promoters_len_perc
dmc.regions.rat["up2000", ]$CpG<- slot(meth.dmc.up2000, "annotation")[1] / up2000_len_perc
dmc.regions.rat["down2000", ]$CpG <- slot(meth.dmc.down2000, "annotation")[1] / down2000_len_perc

meth.dmc.Ann <- annotateWithGeneParts(as(meth.dmc, "GRanges"), gene.parts.up2000)

png(filename = "meth_dmc_by_regions.png", width=6600, height=3000, res=600)
par(mfrow=c(2, 3))
genomation::plotTargetAnnotation(meth.dmc.genes, col=c("green","gray","white"),
                                 main="overlapped with gene regions")
genomation::plotTargetAnnotation(meth.dmc.exons, col=c("green","gray","white"),
                                 main="overlapped with exons")
genomation::plotTargetAnnotation(meth.dmc.introns, col=c("green","gray","white"),
                                 main="overlapped with introns")
genomation::plotTargetAnnotation(meth.dmc.up2000, col=c("green","gray","white"),
                                 main="overlapped with 2000bp upstream from TSS")
genomation::plotTargetAnnotation(meth.dmc.down2000, col=c("green","gray","white"),
                                 main="overlapped with 2000bp downstream from TSS")
genomation::plotTargetAnnotation(meth.dmc.promoters, col=c("green","gray","white"),
                                 main="overlapped with promoters\n(1000bp up and downstream around TSS)")
dev.off()

CHG10.dmc <- getData(meth.CHG.cov10) %>% filter(Ehx1516 <= 0.1 & Ehx217 >= 0.9)
CHG10.dmc <- rbind(CHG10.dmc, getData(meth.CHG.cov10) %>% filter(Ehx1516 >= 0.9 & Ehx217 <= 0.1))

CHG10.dmc.genes <- annotateWithFeature(as(CHG10.dmc, "GRanges"), genes$Gene, feature.name = "gene regions")

CHG10.dmc.exons <- annotateWithFeature(as(CHG10.dmc, "GRanges"), gene.parts$exons, feature.name = "exons")

CHG10.dmc.introns <- annotateWithFeature(as(CHG10.dmc, "GRanges"), gene.parts$introns, feature.name = "introns")

CHG10.dmc.up2000 <- annotateWithFeature(as(CHG10.dmc, "GRanges"), up2000, feature.name = "up2000")
CHG10.dmc.down2000 <- annotateWithFeature(as(CHG10.dmc, "GRanges"), down2000, feature.name = "down2000")
CHG10.dmc.promoters <- annotateWithFeature(as(CHG10.dmc, "GRanges"), gene.parts$promoters, 
                                           feature.name = "promoters")

CHG10.dmc.Ann <- annotateWithGeneParts(as(CHG10.dmc, "GRanges"), gene.parts.up2000)

dmc.regions.rat["genes", ]$CHG <- slot(CHG10.dmc.genes, "annotation")[1] / gene_len_perc
dmc.regions.rat["exons", ]$CHG<- slot(CHG10.dmc.exons, "annotation")[1] / exons_len_perc
dmc.regions.rat["introns", ]$CHG<- slot(CHG10.dmc.introns, "annotation")[1] / introns_len_perc
dmc.regions.rat["promoters", ]$CHG<- slot(CHG10.dmc.promoters, "annotation")[1] / promoters_len_perc
dmc.regions.rat["up2000", ]$CHG<- slot(CHG10.dmc.up2000, "annotation")[1] / up2000_len_perc
dmc.regions.rat["down2000", ]$CHG <- slot(CHG10.dmc.down2000, "annotation")[1] / down2000_len_perc
rownames(dmc.regions.rat) <- seq(1,6)
region_order <- dmc.regions.rat$region
#ggplot(dmc.regions.rat) + geom_point(aes(x = factor(region, level = region_order), y=CpG))
ggplot(dmc.regions.rat) + geom_point(aes(x = factor(region, level = region_order), y=CpG, colour="CpG")) +
  geom_line(aes(x = seq(1,6), y=CpG, colour="CpG")) +
  geom_point(aes(x = factor(region, level = region_order), y=CHG, colour="CHG")) +
  geom_line(aes(x = seq(1,6), y=CHG, colour="CHG")) +
  scale_colour_manual("", breaks = c("CpG", "CHG"), values = c("CpG"= "red", "CHG" = "blue")) +
  geom_hline(yintercept=1.0, linetype="dashed") +
  theme_light() +
  labs(x = "Feature Region", y = "Relative DMC Rate vs. Feature Length") 
  #theme(axis.title=element_text(size=14), plot.title=element_text(size=18, hjust=0.5)) +
  #theme(legend.position = c(0.8, 0.85), legend.direction = "horizontal", legend.text = element_text(size=14))

ggsave("DMC_region_rate.png", dpi=600)

png(filename = "CHG10_dmc_by_regions.png", width=6600, height=3000, res=600)
par(mfrow=c(2, 3))
genomation::plotTargetAnnotation(CHG10.dmc.genes, col=c("blue","gray","white"),
                                 main="overlapped with gene regions")
genomation::plotTargetAnnotation(CHG10.dmc.exons, col=c("blue","gray","white"),
                                 main="overlapped with exons")
genomation::plotTargetAnnotation(CHG10.dmc.introns, col=c("blue","gray","white"),
                                 main="overlapped with introns")
genomation::plotTargetAnnotation(CHG10.dmc.up2000, col=c("blue","gray","white"),
                                 main="overlapped with 2000bp upstream from TSS")
genomation::plotTargetAnnotation(CHG10.dmc.down2000, col=c("blue","gray","white"),
                                 main="overlapped with 2000bp downstream from TSS")
genomation::plotTargetAnnotation(CHG10.dmc.promoters, col=c("blue","gray","white"),
                                 main="overlapped with promoters\n(1000bp up and downstream around TSS)")
dev.off()

# Plot the ratio of  DMC numbers / scaffold lengths
scaffolds_100$CpG_DMC <- 0
for (i in 1:nrow(scaffolds_100)) {
  nDMC <- meth.dmc %>% filter(chr == scaffolds_100[i, ]$chr) %>% nrow()
  scaffolds_100[i, ]$CpG_DMC <- nDMC
}
scaffolds_100$CHG_DMC <- 0
for (i in 1:nrow(scaffolds_100)) {
  nDMC <- CHG10.dmc %>% filter(chr == scaffolds_100[i, ]$chr) %>% nrow()
  scaffolds_100[i, ]$CHG_DMC <- nDMC
}

ggplot(scaffolds_100) + geom_line(mapping=aes(x=seq(1, 100), y = CpG_DMC/len*1000, colour="CpG")) +
  geom_line(mapping=aes(x=seq(1, 100), y = CHG_DMC/len*1000, colour="CHG")) +
  geom_hline(aes(yintercept=mean(CpG_DMC/len)*1000), color="red", linetype="dashed") +
  scale_colour_manual("", breaks = c("CpG", "CHG"), values = c("CpG"= "red", "CHG" = "blue")) +
  geom_hline(aes(yintercept=mean(CHG_DMC/len)*1000), color="blue", linetype="dashed") +
  theme_light() + #  ggtitle("Numbers of DMCs per 1000 bp by Scaffolds") +
  labs(x = "Scaffolds", y = "Avg Num of DMCs per 1000 bp") +
  scale_y_continuous(limits=c(0, 1.2), n.breaks = 7) +
  theme(axis.title=element_text(size=14), plot.title=element_text(size=18, hjust=0.5)) +
  theme(legend.position = c(0.8, 0.85), legend.direction = "horizontal", legend.text = element_text(size=14))

ggsave("DMC_Scaffold_Ratio.png", dpi=600)


## Identify DMRs using a sliding window approach 
## A DMR contains at least 5(?) DMCs in a fixed-length window
meth.dmc.gr <- as(meth.dmc, "GRanges")
dmr_min_len <- 300
dmr_min_cnt <- 3
meth.dmc.flank <- flank(meth.dmc.gr, -dmr_min_len) # 1000 bp window to the right
ovp <- findOverlaps(meth.dmc.flank, meth.dmc.gr)
meth.dmc.flank$nDMC <- countQueryHits(ovp)
# keep the window with at least 5 DMCs
meth.dmr <- meth.dmc.flank %>% as.data.frame() %>% filter(nDMC >= dmr_min_cnt)
meth.dmr <- as(meth.dmr, "GRanges")
# merge overlapped DMRs
meth.dmr <- GenomicRanges::reduce(meth.dmr)
ovp <- findOverlaps(meth.dmr, meth.dmc.gr)
meth.dmr$nDMC <- countQueryHits(ovp)

# meth.dmr.genes <- annotateWithFeature(as(meth.dmr,"GRanges"), genes$Gene, 
#                                       feature.name = "gene regions")
# 
# meth.dmc.exons <- annotateWithFeature(as(meth.dmc, "GRanges"), gene.parts$exons, feature.name = "exons")
# 
# meth.dmc.introns <- annotateWithFeature(as(meth.dmc, "GRanges"), gene.parts$introns, feature.name = "introns")
# 
# meth.dmr.up2000 <- annotateWithFeature(as(meth.dmr, "GRanges"), up2000, feature.name = "up2000")
# meth.dmc.down2000 <- annotateWithFeature(as(meth.dmc, "GRanges"), down2000, feature.name = "down2000")
# meth.dmc.promoters <- annotateWithFeature(as(meth.dmc, "GRanges"), gene.parts$promoters, 
#                                           feature.name = "promoters")

#meth.dmr.o.gp <- overlap_DMR_GeneFeatures(as.data.frame(meth.dmr), gene.parts, expr.coef)

CHG10.dmc.gr <- as(CHG10.dmc, "GRanges")
CHG10.dmc.flank <- flank(CHG10.dmc.gr, -1000)
ovp <-findOverlaps(CHG10.dmc.flank, CHG10.dmc.gr)
CHG10.dmc.flank$nDMC <- countQueryHits(ovp)
# keep the window with at least 5 DMCs
CHG10.dmr <- CHG10.dmc.flank %>% as.data.frame() %>% filter(nDMC >= 5)
CHG10.dmr <- as(CHG10.dmr, "GRanges")
# merge overlapped DMRs
CHG10.dmr <- GenomicRanges::reduce(CHG10.dmr)
ovp <- findOverlaps(CHG10.dmr, CHG10.dmc.gr)
CHG10.dmr$nDMC <- countQueryHits(ovp)

CHG10.dmr.o.gp <- overlap_DMR_GeneFeatures(as.data.frame(CHG10.dmr), gene.parts, expr.coef)
