###############################################################################
# This script analyze the relation of differentially expressed genes (DEGs) with 
# DMCs in the their gene, promoter
# and upstream regions 
#
###############################################################################
library(rethinking)
library(tidyverse)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# set up gene_regions in ../methyl_by_regions.R
gene_regions

# remove duplicates in gene_regions
gene_regions_uq <- gene_regions %>% distinct(ID, .keep_all = TRUE)
# load gene expr data from Bayesian analysis
load("./RData/deg.RData")

# Calculate scaling factor using DESeq2
library(DESeq2)
samples <- data.frame(samples = colnames(deg_data[, 2:7]))
samples$strain <- c(rep("EH1516", 3), rep("EH217",3))
dds <- DESeqDataSetFromMatrix(countData = deg_data[, 2:7], colData=samples, design = ~1)
dds <- estimateSizeFactors(dds)
log(sizeFactors(dds))

#fix col names
colnames(deg_data)[8:9] <- c("expr_base", "std")
colnames(deg_data)[14:15] <- c("logfc", "std.1")
# merge gene_regions with DEG data
gene_data <- merge(gene_regions_uq, deg_data, by = "ID", no.dups=TRUE)

load("./RData/dmc.RData")
# Make a list of meth data, separate by chr 
#dmc_data_list <- list()
#for(chr in  unique(dmc_data_all$chr)) {
#  dmc_data_list[[chr]] <- dmc_data_all[dmc_data_all$chr == chr, ]
#}
# load the meth data list from file
load("./RData/dmc_list.RData")
# Get meth data in one gene region
meth_in_region <- function(gr) {
  nr <- nrow(gr)
  m_count = data.frame(ncpg = rep(0, nr), dmc_neg = rep(0, nr), dmc_pos = rep(0, nr),
                       b1 = rep(0, nr), b2 = rep(0, nr),
                       dmc_b1 = rep(0, nr), dmc_b2 = rep(0, nr))
  for (i in 1:nr) {
    chr <- as.character(gr[i, 1])
    start <- gr[i, 2]
    end <- gr[i, 3]
    #print(paste(chr, start, end))
    meth_d <- dmc_data_list[[chr]]
    if(!is.null(meth_d)) { 
      m_in_g <- meth_d[meth_d$start >= start & meth_d$end <= end, ]
      nc <- nrow(m_in_g)
      n_neg <- sum(m_in_g$dmc & ((m_in_g$m2-m_in_g$m1) < 0))
      n_pos <- sum(m_in_g$dmc & ((m_in_g$m2-m_in_g$m1) > 0))
      
      # calculate the mean meth levels for all CpGs in the region for two strains
      b1 <- 0
      b2 <- 0
      if(nc > 0) {
        b1 <- (sum(m_in_g$numCs1) + sum(m_in_g$numCs2)) / 
          ((sum(m_in_g$coverage1) + sum(m_in_g$coverage2)))
        b2 <- (sum(m_in_g$numCs3) + sum(m_in_g$numCs4) + sum(m_in_g$numCs5)) / 
          (sum(m_in_g$coverage3) + sum(m_in_g$coverage4) + sum(m_in_g$coverage5))
      }
      # calculate the mean meth levels for two strains for only DMCs
      dmc_b1 <- 0
      dmc_b2 <- 0
      if((n_pos + n_neg) > 0) {
        # DMCs in the region
        dmc_in_g <- m_in_g %>% filter (dmc)
        dmc_b1 <- (sum(dmc_in_g$numCs1) + sum(dmc_in_g$numCs2)) / 
             ((sum(dmc_in_g$coverage1) + sum(dmc_in_g$coverage2)))
        dmc_b2 <- (sum(dmc_in_g$numCs3) + sum(dmc_in_g$numCs4) + sum(dmc_in_g$numCs5)) / 
          (sum(dmc_in_g$coverage3) + sum(dmc_in_g$coverage4) + sum(dmc_in_g$coverage5))
      }
      m_count[i, ] <- c(nc, n_neg, n_pos, b1, b2, dmc_b1, dmc_b2)
    }
  }
  
  return(m_count)
}

get_DMC_in_regions <- function(regions) {
  library(methylKit)
  data <- merge(regions, deg_data, by="ID", no.dups=TRUE)
  # count # of methyl sites in gene regions
  x <- meth_in_region(data[, 2:4])
  data <- cbind(data, x)
  # selec regions with non-zero CpG sites
  data_dmc <- data[data$ncpg > 0, ]
  # Only keep expressed genes that have a least 3 nonzero counts in expression data 
  data_dmc <- data_dmc %>% filter(rowSums(data_dmc[, 8:13] > 0) >= 3)
  # dens((data_dmc$dmc_neg+data_dmc$dmc_pos) *1000 / data_dmc$width)
  # plot(data_dmc$dmc_neg/data_dmc$ncpg, data_dmc$dmc_pos/data_dmc$ncpg)
  # x <- (data_dmc$dmc_neg) / data_dmc$ncpg
  # which (x == max(x))
  # x2 <- log(gene_data_dmc$dmc_pos +1) 
  
  # get average methylation level per region by strains
  # load("RData/methobj.RData")
  # sample.ids = c("EH1516B", "EH1516C", "EH217A", "EH217B", "EH217C")
  # p <- getFeatureMethyl(methobj, as(data_dmc, "GRanges"), sample.ids, lo.count = 3)
  # m <- p$m[as.character(data_dmc$ID), ]      # reorder the data by gene_data_dmc ID's
  # m1 <- rowMeans(m[, 1:2])                   # average m values of EH1516
  # m2 <- rowMeans(m[, 3:5])
  # data_dmc$m1 <- m1
  # data_dmc$m2 <- m2
  # 
  # b <- p$beta[as.character(data_dmc$ID), ]
  # b1 <- rowMeans(b[, 1:2])
  # b2 <- rowMeans(b[, 3:5])
  # #s <- normalize(c(b1, b2))
  # data_dmc$b1 <- b1
  # data_dmc$b2 <- b2
  
  return(data_dmc)
}

# Get methylation data in gene regions
gene_data_dmc <- get_DMC_in_regions(gene_regions_uq)

# Test with a sample of gene_data_dmc
# colnames(gene_data_dmc)[1] <- "gid"
n <- sample(1:nrow(gene_data_dmc), 1000)
#d <- gene_data_dmc[n, ]
d <- gene_data_dmc

# reshape the data so that each row has one gene expression
# the first col assumed to be gid
# the reshaped data start from col2 and have nsamples
sf_v <- log(sizeFactors(dds))            # scaling factor vector of size 6
strain <- c(rep(0, 3), rep(1, 3))  # 1: EH1516, 2: EH217
reshape_data <- function(d, nsample, strain) {
  # the first sample
  d_reshaped <- d[, c(1, 2)]
  colnames(d_reshaped) <- c("gid", "expr")
  d_reshaped$sample <- 1
  d_reshaped$strain <- strain[1]
  d_reshaped$sf <- sf_v[1]
  for (i in 2:nsample) {
    d_i <- d[, c(1, i+1)]
    colnames(d_i) <- c("gid", "expr")
    d_i$sample <- i
    d_i$sf <- sf_v[i]
    d_i$strain <- strain[i]
    d_reshaped <- rbind(d_reshaped, d_i)
  }
  d_reshaped$id <- rep(1:nrow(d), nsample)
  return (d_reshaped)
}

# Test the relation between the count of pos and neg DMCs with the gene expression
test_Expr_DMC <- function(d) {
  d$id <- 1:nrow(d)
  d <- d[, c(1,8:13, 2:7, 26:30)]
  d_reshaped <- reshape_data(d, 6, strain)
  d_reshaped$ncpg <- rep(d$ncpg, 6)
  d_reshaped$dmc_pos <- rep(d$dmc_pos, 6)
  d_reshaped$dmc_neg <- rep(d$dmc_neg, 6)
  d_reshaped$DE <- rep(d$DE, 6)
  
  e_bar = log(median(d_reshaped$expr))
  # build the dat list for MCMC
  dat <- list (
    G = d_reshaped$id,              # gene id
    E = d_reshaped$expr,            # expression count
    S = d_reshaped$sample,          # sample ID
    T = d_reshaped$strain,          # treatment, i.e. strain = 0 for EH1516, 1 for EH217
    sf = d_reshaped$sf,             # scaling factor
    #W = standardize(log(rep(d$width, 6))),          # gene width
    e_bar = rep(e_bar, nrow(d_reshaped)),           # global mean expressionb
    XC = normalize(log(d_reshaped$ncpg)),           # normalized log number of CpG sites
    XDP = normalize(log(d_reshaped$dmc_pos + 1)),   # normalized log number of dmc_pos sites
    XDN = normalize(log(d_reshaped$dmc_neg + 1)),   # normalized log number of dmc_neg sites
    #XDP = normalize(d_reshaped$dmc_pos / d_reshaped$ncpg),
    #XDN = normalize(d_reshaped$dmc_neg / d_reshaped$ncpg),
    ng = max(d_reshaped$id)
  )
  
  # Model the relation between gene expression and DMCs
  m_Expr_DMC <- ulam(
    alist (
      E ~ dgampois(lambda, phi),
      # f[S]: log sample factor, e[G] log express of gene in treatment 1,
      log(lambda) <- e_bar + f[S] + e[G] + bDP * XDP * T + bDN * XDN *T,
      #log(lambda) <- e_bar + sf + e[G] + bDP * XDP * T + bDN * XDN *T,
      
      vector[6]: f ~ normal(0, 0.2),
      vector[ng]: e ~ normal(0, 3),
      bDP ~ normal(0, 1.5),
      bDN ~ normal(0, 1.5),
      phi ~ dexp(1)
    ), data = dat, chains=4, cores=4
  )
  return(m_Expr_DMC)
}

m_Expr_DMC <- test_Expr_DMC(gene_data_dmc[, 1:29])

m_Expr_DMC_gt_0 <- test_Expr_DMC(gene_data_dmc[, 1:29] %>% filter (dmc_pos > 0 | dmc_neg >0))
precis(m_Expr_DMC, depth=2, prob=0.95,
       pars = c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]", "f[6]"))
precis(m_Expr_DMC_gt_0)
precis(m_Expr_DMC_gt_0, depth=2, prob=0.95,
       pars = c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]", "f[6]"))

# Multi-level Model for gene expression and DMCs
# Only consider genes with non-zero DMCs
d <- gene_data_dmc %>% filter (dmc_pos > 0 | dmc_neg > 0)
d <- d[1:1000,]
d <- gene_data_dmc
d <- add_column(d, id = 1:nrow(d), .after = 29)
#d <- d[, c(1,8:13, 2:7, 26:30)]
d_reshaped <- reshape_data(d[, c(1,8:13, 2:7, 26:30)], 6, strain)
d_reshaped$ncpg <- rep(d$ncpg, 6)
d_reshaped$dmc_pos <- rep(d$dmc_pos, 6)
d_reshaped$dmc_neg <- rep(d$dmc_neg, 6)
d_reshaped$DE <- rep(d$DE, 6)

e_bar = log(median(d_reshaped$expr))
# build the dat list for MCMC
dat <- list (
  G = d_reshaped$id,
  E = d_reshaped$expr,
  S = d_reshaped$sample,
  T = d_reshaped$strain,
  sf = d_reshaped$sf,
  W = standardize(log(rep(d$width, 6))),
  e_bar = rep(e_bar, nrow(d_reshaped)),
  XC = normalize(log(d_reshaped$ncpg)),           # normalized log number of CpG sites
  XDP = normalize(log(d_reshaped$dmc_pos + 1)),   # normalized log number of dmc_pos sites
  XDN = normalize(log(d_reshaped$dmc_neg + 1)),   # normalized log number of dmc_neg sites
  #XDP = normalize(d_reshaped$dmc_pos / d_reshaped$ncpg),
  #XDN = normalize(d_reshaped$dmc_neg / d_reshaped$ncpg),
  ng = max(d_reshaped$id)
)

m_ML_Expr_DMC_new2 <- ulam(
  alist (
    E ~ dgampois(lambda, phi),
    
    #log(lambda) <- e_bar + f[S] + e[G] + bDP[G] * XDP * T + bDN[G] * XDN *T,
    #vector[6]: f ~ normal(0, 0.2),
    log(lambda) <- e_bar + sf + e[G] + bDP[G] * XDP * T + bDN[G] * XDN *T,
    vector[ng]: e ~ normal(0, 3),
    phi ~ exponential(1),
    #ebar ~ normal(0, 1.5),
    #sigma ~ exponential(1),
    vector[ng]: bDP ~ normal(bDP_bar, bDP_sigma),
    vector[ng]: bDN ~ normal(bDN_bar, bDN_sigma),
    c(bDP_bar, bDN_bar) ~ normal(0, 1.5),
    c(bDP_sigma, bDN_sigma) ~ exponential(1)
  ), data = dat, chains = 4, cores = 4
)
precis(m_ML_Expr_DMC)

precis(m_ML_Expr_DMC_new)

# Test GLM using difference of M values
test_Expr_Mval <- function(d) {
  d <- add_column(d, id = 1:nrow(d), .after = 29)
  #d <- d[, c(1,8:13, 2:7, 26:30)]
  d_reshaped <- reshape_data(d[, c(1,8:13, 2:7, 26:30)], 6, strain)
  d_reshaped$ncpg <- rep(d$ncpg, 6)
  d_reshaped$dmc_pos <- rep(d$dmc_pos, 6)
  d_reshaped$dmc_neg <- rep(d$dmc_neg, 6)
  d_reshaped$DE <- rep(d$DE, 6)
  
  e_bar = log(median(d_reshaped$expr))
  
  # Data list using average methylation level as input
  dat <- list (
    G = d_reshaped$id,
    E = d_reshaped$expr,
    S = d_reshaped$sample,
    T = d_reshaped$strain,
    e_bar = rep(e_bar, nrow(d_reshaped)),
    #m1 = rep(d$m1, 6),                                   # standardized methylation m values in EH1516
    #m2 = rep(d$m2, 6),                                   # standardized methylation m values in EH1516
    dm = c(rep(d$m1-d$m1, 3), rep(d$m2-d$m1, 3)),         # the first 3 samples from strain1 and last 3 from strain 2 
    ng = max(d_reshaped$id)
  )
  
  # Model the relation between gene expression and change of average methylation values
  m_Expr_M <- ulam(
    alist (
      E ~ dgampois(lambda, phi),
      # f[S]: log sample factor, e[G] log express of gene in treatment 1,
      log(lambda) <- e_bar + f[S] + e[G] + bM * dm,
      
      vector[6]: f ~ normal(0, 0.1),
      vector[ng]: e ~ normal(0, 3),
      bM ~ normal(0, 1.5),
      phi ~ dexp(1)
    ), data = dat, chains=4, cores=4
  )
  
  return(m_Expr_M)
}

# Test GLM using difference of B (meth) values
test_Expr_Bval <- function(d) {
  d <- add_column(d, id = 1:nrow(d), .after = 29)
  #d <- d[, c(1,8:13, 2:7, 26:30)]
  d_reshaped <- reshape_data(d[, c(1,8:13, 2:7, 26:30)], 6, strain)
  d_reshaped$ncpg <- rep(d$ncpg, 6)
  d_reshaped$dmc_pos <- rep(d$dmc_pos, 6)
  d_reshaped$dmc_neg <- rep(d$dmc_neg, 6)
  d_reshaped$DE <- rep(d$DE, 6)
  
  e_bar = log(median(d_reshaped$expr))
  
  # Data list using average methylation level as input
  dat <- list (
    G = d_reshaped$id,
    E = d_reshaped$expr,
    S = d_reshaped$sample,
    T = d_reshaped$strain,
    e_bar = rep(e_bar, nrow(d_reshaped)),
    #m1 = rep(d$m1, 6),                                   # standardized methylation m values in EH1516
    #m2 = rep(d$m2, 6),                                   # standardized methylation m values in EH1516
    db = c(rep(d$b1-d$b1, 3), rep(d$b2-d$b1, 3)),         # the first 3 samples from strain1 and last 3 from strain 2 
    ng = max(d_reshaped$id)
  )
  
  # Model the relation between gene expression and change of average methylation values
  m_Expr_B <- ulam(
    alist (
      E ~ dgampois(lambda, phi),
      # f[S]: log sample factor, e[G] log express of gene in treatment 1,
      log(lambda) <- e_bar + f[S] + e[G] + bB * db,
      
      vector[6]: f ~ normal(0, 0.1),
      vector[ng]: e ~ normal(0, 3),
      bB ~ normal(0, 1.5),
      phi ~ dexp(1)
    ), data = dat, chains=4, cores=4
  )
  
  return(m_Expr_B)
}

m_Expr_M_gt_0 <- test_Expr_Mval(gene_data_dmc %>% filter(dmc_neg > 0 | dmc_pos >0 ))

m_Expr_B_gt_0 <- test_Expr_Bval(gene_data_dmc %>% filter((dmc_neg+dmc_pos) > 0))


# Get DMCs in the promoter regions
# get GRanges representing all scaffolds
scaffolds.gr = read_scaffold_gr("data/Ehux_genome.fasta.len")

# Use functions in methyl_by_regions.R
promoters <- getPromoterRegions(up=1000, down=1000, bed_file = "data/Ehux_genbank.bed") 
promoters <- trimRegionByScaffolds(promoters)

prom_data_dmc <- get_DMC_in_regions(promoters)
# drop column 7
prom_data_dmc[, 7] = NULL

m_Expr_DMC_Prom <- test_Expr_DMC(prom_data_dmc[, 1:29])

m_Expr_DMC_Prom_gt_0 <- test_Expr_DMC(prom_data_dmc[, 1:29] %>% filter (dmc_neg > 0 | dmc_pos > 0))

m_Expr_M_Prom <- test_Expr_Mval(prom_data_dmc %>% filter (dmc_neg >0 | dmc_pos > 0))

# Get DMCs in the up2000 regions
up2000 <- getPromoterRegions(up = 2000, down = 0, bed_file = "data/Ehux_genbank.bed")
up2000 <- trimRegionByScaffolds(up2000)
up2000_data_dmc <- get_DMC_in_regions(up2000)
up2000_data_dmc[, 7] <- NULL

m_Expr_DMC_up2000 <- test_Expr_DMC(up2000_data_dmc[, 1:29])
m_Expr_DMC_up2000_gt_0 <- test_Expr_DMC(up2000_data_dmc[, 1:29] %>% filter (dmc_neg > 0 | dmc_pos > 0))

m_Expr_M_up2000 <- test_Expr_Mval(up2000_data_dmc %>% filter (dmc_neg >0 | dmc_pos > 0))

# Get DMCs in the up2000 regions
up1000 <- getPromoterRegions(up = 1000, down = 0, bed_file = "data/Ehux_genbank.bed")
up1000 <- trimRegionByScaffolds(up1000)
up1000_data_dmc <- get_DMC_in_regions(up1000)
up1000_data_dmc[, 7] <- NULL

m_Expr_DMC_up1000 <- test_Expr_DMC(up1000_data_dmc[, 1:29])
m_Expr_DMC_up1000_gt_0 <- test_Expr_DMC(up1000_data_dmc[, 1:29] %>% filter (dmc_neg > 0 | dmc_pos > 0))

m_Expr_M_up1000 <- test_Expr_Mval(up1000_data_dmc %>% filter (dmc_neg >0 | dmc_pos > 0))
