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
  m_count = data.frame(ncpg = rep(0, nr), dmc_neg = rep(0, nr), dmc_pos = rep(0, nr))
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
      m_count[i, ] <- c(nc, n_neg, n_pos)
    }
  }
  
  return(m_count)
}

# count # of methyl sites in gene regions
x <- meth_in_region(gene_data[, 2:4])
gene_data <- cbind(gene_data, x)
# select genes with non-zero Cpg sites
gene_data_dmc <- gene_data[gene_data$ncpg >0, ]
dens((gene_data_dmc$dmc_neg+gene_data_dmc$dmc_pos) *1000 / gene_data_dmc$width)
plot(gene_data_dmc$dmc_neg/gene_data_dmc$ncpg, gene_data_dmc$dmc_pos/gene_data_dmc$ncpg)
x <- (gene_data_dmc$dmc_neg) / gene_data_dmc$ncpg
which (x == max(x))
x2 <- log(gene_data_dmc$dmc_pos +1) 

# Only keep expressed genes that have a least 3 nonzero counts in expression data 
gene_data_dmc <- gene_data_dmc %>% filter(rowSums(gene_data_dmc[, 8:13] > 0) >= 3)

# get average methylation level per gene by strains
load("RData/methobj.RData")
sample.ids = c("EH1516B", "EH1516C", "EH217A", "EH217B", "EH217C")
p <- getFeatureMethyl(methobj, as(gene_data_dmc, "GRanges"), sample.ids, lo.count = 3)
m <- p$m[as.character(gene_data_dmc$ID), ] # reorder the data by gene_data_dmc ID's
m1 <- rowMeans(m[, 1:2])        # average m values of EH1516
m2 <- rowMeans(m[, 3:5])
#ng <- length(m1)
#s <- standardize(c(m1, m2))
#gene_data_dmc$m1 <- s[1:ng]
#gene_data_dmc$m2 <- s[(ng+1):(2*ng)]
gene_data_dmc$m1 <- m1
gene_data_dmc$m2 <- m2

b <- p$beta[as.character(gene_data_dmc$ID), ]
b1 <- rowMeans(b[, 1:2])
b2 <- rowMeans(b[, 3:5])
#s <- normalize(c(b1, b2))
#gene_data_dmc$b1 <- s[1:ng]
#gene_data_dmc$b2 <- s[(ng+1):(2*ng)]


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

# Test the relation between pos and neg DMCs with the gene expression
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
    W = standardize(log(rep(d$width, 6))),          # gene width
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
precis(m_Expr_DMC)

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
m_Expr_M2 <- ulam(
  alist (
    E ~ dgampois(lambda, phi),
    # f[S]: log sample factor, e[G] log express of gene in treatment 1,
    log(lambda) <- e_bar + f[S] + e[G] + bM * (m2 - m1),
    
    vector[6]: f ~ normal(0, 0.1),
    vector[ng]: e ~ normal(0, 3),
    bM ~ normal(0, 1.5),
    phi ~ dexp(1)
  ), data = dat, chains=4, cores=4
)


# Get DMCs in the promoter regions
#promoters <- as.data.frame(gene.parts$promoters)
#promoters$ID = mapping[promoters$name, "ID"]
#promoters = promoters[!duplicated(promoters$ID), ]    # remove duplicates due to alternative splicing
#colnames(promoters)[1] <- "chr"
prom_data <- merge(promoters, deg_data, by = "ID", no.dups=TRUE)
# drop column 7
prom_data[, 7] = NULL
x <- meth_in_region(prom_data[, 2:4])
prom_data <- cbind(prom_data, x)
prom_data_dmc <-prom_data[prom_data$ncpg >0, ]
dens((prom_data_dmc$dmc_neg+prom_data_dmc$dmc_pos)  / prom_data_dmc$width)
x <- (prom_data_dmc$dmc_neg) / prom_data_dmc$ncpg
which (x == max(x))

# filter out genes with at least 3 zero counts
prom_data_dmc <- prom_data_dmc %>% filter(rowSums(prom_data_dmc[, 8:13] > 0) >= 3)

m_Expr_DMC_Prom <- test_Expr_DMC(prom_data_dmc)

# Get DMCs in the up2000 regions
up2000_data <- merge(up2000, deg_data, by="ID", no.dups=TRUE)
up2000_data[, 7] = NULL
x <- meth_in_region(up2000_data[, 2:4])
up2000_data <- cbind(up2000_data, x)
up2000_data_dmc <- up2000_data[up2000_data$ncpg > 0, ]
dens((up2000_data_dmc$dmc_neg+up2000_data_dmc$dmc_pos)  / up2000_data_dmc$width)
# filter out genes with at least 3 zero counts
up2000_data_dmc <- up2000_data_dmc %>% filter(rowSums(up2000_data_dmc[, 8:13] > 0) >= 3)

m_Expr_DMC_up2000 <- test_Expr_DMC(up2000_data_dmc)

# Get DMCs in the up2000 regions
up1000_data <- merge(up1000, deg_data, by="ID", no.dups=TRUE)
up1000_data[, 7] = NULL
x <- meth_in_region(up1000_data[, 2:4])
up1000_data <- cbind(up1000_data, x)
up1000_data_dmc <- up1000_data[up1000_data$ncpg > 0, ]
dens((up1000_data_dmc$dmc_neg+up1000_data_dmc$dmc_pos)  / up1000_data_dmc$width)
# filter out genes with at least 3 zero counts
up1000_data_dmc <- up1000_data_dmc %>% filter(rowSums(up1000_data_dmc[, 8:13] > 0) >= 3)

m_Expr_DMC_up1000 <- test_Expr_DMC(up1000_data_dmc)
