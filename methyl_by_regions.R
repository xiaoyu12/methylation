# Get Granges object representing regions upstream and downstream from TSSes with given len
getPromoterRegions <- function(up = 1000, down = 1000, bed_file="./Ehux_genbank.bed") {
  library(genomation)
  gene_structs <- readTranscriptFeatures(bed_file, remove.unusual = FALSE,
                                         up.flank = up, down.flank = down, 
                                         unique.prom = FALSE)
  #promoters <- addMapID(gene_structs$promoter, mapping)
  ranges <- gene_structs$promoter
  ranges$ID = mapping[ranges$name, "ID"]
  ranges <- as.data.frame(ranges)
  ranges = ranges[!duplicated(ranges$ID), ]    # remove duplicates due to alternative splicing
  # convert back to Granges
  ranges <- as(ranges, "GRanges")
  return(ranges)
}

# get GRanges representing all scaffolds
read_scaffold_gr <- function(scaf_len_file="Ehux_genome.fasta.len") {
  scaffolds <- read.table(scaf_len_file)
  colnames(scaffolds) <- c("chr", "len")
  genome_len <- sum(scaffolds$len)
  rownames(scaffolds) <- scaffolds$chr
  scaffolds.gr <- scaffolds[, 1:2]
  scaffolds.gr[, 2] <-1 
  scaffolds.gr[, 3] <- scaffolds$len
  colnames(scaffolds.gr) <- c("chr", "start", "end")
  scaffolds.gr <- as(scaffolds.gr, "GRanges")
  return(scaffolds.gr)
}
scaffold_gr = read_scaffold_gr("data/Ehux_genome.fasta.len")

# trim regions based the starts and ends of scaffolds
trimRegionByScaffolds <- function(gr) {
  # Get the promoter regions overlapped with the genome
  o = findOverlaps(gr, scaffolds.gr)
  grl1 = split(gr[queryHits(o)], 1:length(o)) # You can't mendoapply on a GRanges object
  grl2 = split(scaffolds.gr[subjectHits(o)], 1:length(o))
  # set the promoters min and max to the min and max of genome scaffold
  foo = function(x, y) {
    rv = x
    start(rv) = max(start(x), start(y))
    end(rv) = min(end(x), end(y))
    return(rv)
  }
  gr_trimmed = unlist(mendoapply(foo, grl1, y=grl2))
  return(gr_trimmed)
}

# set up regions of interest
library(genomation)
mapping <- read.table("./genbank_mapping.txt", row.names = 1)
colnames(mapping) <- c("gene", "ID")
gene.parts <- readTranscriptFeatures("./Ehux_genbank.bed", remove.unusual = FALSE, unique.prom = FALSE)

genes = readFeatureFlank("./Ehux_genbank.bed", feature.flank.name=c("Gene","UTR"))
genes$Gene$name <- gene.parts$TSSes$name
gene_regions <- as.data.frame(genes$Gene)
gene_regions$name = gene.parts$TSSes$name  
gene_regions$ID <- mapping[gene_regions$name, "ID"]

up2000 <- getPromoterRegions(up = 2000, down = 0, bed_file = "data/Ehux_genbank.bed")
up1000 <- getPromoterRegions(up = 1000, down = 0, bed_file = "data/Ehux_genbank.bed")
down2000 <- getPromoterRegions(up = 0, down = 2000, bed_file = "data/Ehux_genbank.bed")
promoters <- getPromoterRegions(up=1000, down=1000, bed_file = "data/Ehux_genbank.bed")
# promoters <- as.data.frame(gene.parts$promoters)
# promoters$ID = mapping[promoters$name, "ID"]
# promoters = promoters[!duplicated(promoters$ID), ]    # remove duplicates due to alternative splicing
# colnames(promoters)[1] <- "chr"
# promoters.gr = as(promoters, 'GRanges')



promoters <- trimRegionByScaffolds(promoters)
up2000 <- trimRegionByScaffolds(up2000)
down2000 <- trimRegionByScaffolds(down2000)
up1000 <- trimRegionByScaffolds(up1000)

# Methylation counts in promoter regions
promobj <- regionCounts(methobj, gene.parts$promoters)
#meth.prom <- methylKit::unite(filterByCoverage(promobj, lo.count = 100))
meth.prom <- getFeatureMethyl(methobj, gene.parts$promoters, sample.ids)

# Compute the mean and confidence inteval of methylation levels by regions
conf_interval_by_regions <- function(meth.data, region="promoters", conf.int=0.95) {
  meth.data$beta <- as.data.frame(meth.data$beta)
  meth.data$beta$EH1516_CpG <- 100 * (meth.data$beta$EH1516B + meth.data$beta$EH1516C) / 2
  meth.data$beta$EH217_CpG <- 100 * (meth.data$beta$EH217A + meth.data$beta$EH217B + 
                                       meth.data$beta$EH217C)/3
  d = data.frame("strain"="EH1516", "region" = region, "mean" = mean(meth.data$beta$EH1516_CpG), 
                 "lci" = t.test(meth.data$beta$EH1516_CpG, conf.level=conf.int)$conf.int[1], 
                 "uci" = t.test(meth.data$beta$EH1516_CpG, conf.level=conf.int)$conf.int[2]
                 )
  
  d1 = data.frame("strain"="EH217", "region" = region, "mean" = mean(meth.data$beta$EH217_CpG), 
                 "lci" = t.test(meth.data$beta$EH217_CpG, conf.level=conf.int)$conf.int[1], 
                 "uci" = t.test(meth.data$beta$EH217_CpG, conf.level=conf.int)$conf.int[2]
                 )
  d = rbind(d, d1)
  return(d)
}


CHG.prom <- getFeatureMethyl(bismark$methCHG, gene.parts$promoters, sample.ids)

# intron regions
intron_obj <- regionCounts(methobj, gene.parts$introns)
#meth.intron <- methylKit::unite(filterByCoverage(intron_obj, lo.count = 100))
meth.intron <- getFeatureMethyl(methobj, gene.parts$introns, sample.ids)
CHG.intron <- getFeatureMethyl(bismark$methCHG, gene.parts$introns, sample.ids)

#exon regions
exon_obj <- regionCounts(methobj, gene.parts$exons)
#meth.exon <- methylKit::unite(filterByCoverage(exon_obj, lo.count = 100))
meth.exon <- getFeatureMethyl(methobj, gene.parts$exons, sample.ids)
CHG.exon <- getFeatureMethyl(bismark$methCHG, gene.parts$exons, sample.ids)


# get methylation data in the upstream 2000 regions
#meth.up2000 <- selectByOverlap(meth, up2000)
meth.up2000 <- getFeatureMethyl(methobj, up2000, sample.ids)
CHG.up2000 <- getFeatureMethyl(bismark$methCHG, up2000, sample.ids)

meth.down2000 <- getFeatureMethyl(methobj, down2000, sample.ids)
CHG.down2000 <- getFeatureMethyl(bismark$methCHG, down2000, sample.ids)

meth_conf <- conf_interval_by_regions(meth.prom, region="promoters")
meth_conf <- rbind(meth_conf, conf_interval_by_regions(meth.exon, region="exons"))
meth_conf <- rbind(meth_conf, conf_interval_by_regions(meth.intron, region = "introns"))
meth_conf <- rbind(meth_conf, conf_interval_by_regions(meth.up2000, region = "up2000"))
meth_conf <- rbind(meth_conf, conf_interval_by_regions(meth.down2000, region = "down2000"))

meth_conf_0.01 <- conf_interval_by_regions(meth.prom, region="promoters", conf.int = 0.99)
meth_conf_0.01 <- rbind(meth_conf_0.01, conf_interval_by_regions(meth.exon, region="exons", conf.int = 0.99))
meth_conf_0.01 <- rbind(meth_conf_0.01, conf_interval_by_regions(meth.intron, region = "introns", conf.int = 0.99))
meth_conf_0.01 <- rbind(meth_conf_0.01, conf_interval_by_regions(meth.up2000, region = "up2000", conf.int = 0.99))
meth_conf_0.01 <- rbind(meth_conf_0.01, conf_interval_by_regions(meth.down2000, region = "down2000", conf.int = 0.99))

CHG_conf <- conf_interval_by_regions(CHG.prom, region="promoters")
CHG_conf <- rbind(CHG_conf, conf_interval_by_regions(CHG.exon, region="exons"))
CHG_conf <- rbind(CHG_conf, conf_interval_by_regions(CHG.intron, region = "introns"))
CHG_conf <- rbind(CHG_conf, conf_interval_by_regions(CHG.up2000, region = "up2000"))
CHG_conf <- rbind(CHG_conf, conf_interval_by_regions(CHG.down2000, region = "down2000"))

CHG_conf_0.01 <- conf_interval_by_regions(CHG.prom, region="promoters", conf.int = 0.99)
CHG_conf_0.01 <- rbind(CHG_conf_0.01, conf_interval_by_regions(CHG.exon, region="exons", conf.int = 0.99))
CHG_conf_0.01 <- rbind(CHG_conf_0.01, conf_interval_by_regions(CHG.intron, region = "introns", conf.int = 0.99))
CHG_conf_0.01 <- rbind(CHG_conf_0.01, conf_interval_by_regions(CHG.up2000, region = "up2000", conf.int = 0.99))
CHG_conf_0.01 <- rbind(CHG_conf_0.01, conf_interval_by_regions(CHG.down2000, region = "down2000", conf.int = 0.99))

# Plot mean and error bars for methylation by regions
library(ggrepel)
ggplot(meth_conf) + geom_point(aes(x=region, y=mean, color=strain), position = position_dodge(width=0.2)) +
  geom_errorbar(aes(x=region, ymin=lci, ymax=uci, color=strain), position = "dodge", width=0.2) +
  #theme_light() +
  ggtitle("CpG Methylation Level Distributions by Regions") +
  labs(x = "Region", y = "Methylation Level (%)") +
  theme(axis.title=element_text(size=14), plot.title=element_text(size=14, hjust=0.5)) 
  #theme(legend.position = c(0.45, 0.6), legend.direction = "horizontal", legend.text = element_text(size=14))
ggsave("meth_Region_Methyl.png", dpi=600)

ggplot(CHG_conf) + geom_point(aes(x=region, y=mean, color=strain), position = position_dodge(width=0.2)) +
  geom_errorbar(aes(x=region, ymin=lci, ymax=uci, color=strain), position = "dodge", width=0.2) +
  #theme_light() +
  ggtitle("CHG Methylation Level Distributions by Regions") +
  labs(x = "Region", y = "Methylation Level (%)") +
  theme(axis.title=element_text(size=14), plot.title=element_text(size=14, hjust=0.5)) 
#theme(legend.position = c(0.45, 0.6), legend.direction = "horizontal", legend.text = element_text(size=14))
ggsave("CHG_Region_Methyl.png", dpi=600)


#meth.up1000 <- selectByOverlap(meth, up1000)

m <- getFeatureMethyl(methobj, gene.parts$promoters, sample.ids)
m_up1000 <- getFeatureMethyl(methobj, gene.parts_up1000$promoters, sample.ids)
m_dn1000 <- getFeatureMethyl(methobj, gene.parts_dn1000$promoters, sample.ids)

m.EH1516.beta <- m$beta[, 1:2]
m.EH217.beta <- m$beta[, 3:5]

m.DMR.beta <- m$beta[which(abs(rowMedians(m.EH1516.beta) - rowMedians(m.EH217.beta)) >= 0.2), ]
#m10 <- getFeatureMethyl(filterByCoverage(methobj, lo.count = 10), gene.parts$promoters, sample.ids)
#m10_up1000 <- getFeatureMethyl(filterByCoverage(methobj, lo.count = 10), gene.parts_up1000$promoters, sample.ids)

m_exon <- getFeatureMethyl(methobj, gene.parts$exons, sample.ids)
m_intron <- getFeatureMethyl(methobj, gene.parts$introns, sample.ids)




getCorrelation(meth.prom, plot = FALSE)
meth.prom.DSS <- calculateDiffMeth(meth.prom)
meth.prom.DSS.20p <- getMethylDiff(meth.prom.DSS, difference = 20, qvalue = 0.01)

meth.prom.DSS.data <- merge(getData(meth.prom.DSS), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
meth.prom.DSS.20p <- merge(getData(meth.prom.DSS.20p), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
plot(meth.prom.DSS.20p$meth.diff, expr.coef[as.character(meth.prom.DSS.20p$ID), ]$LFC, pch=20, xlab="Promoter Methyl Diff", ylab="LFC", col="blue")
text(meth.prom.DSS.20p$meth.diff, expr.coef[as.character(meth.prom.DSS.20p$ID), ]$LFC+0.5, labels = as.character(meth.prom.DSS.20p$ID), cex=0.5)
abline(h=0)
abline(v=0)

