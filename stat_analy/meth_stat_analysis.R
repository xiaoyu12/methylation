# Script to perform Bayesian stat analysis on E. hux methylation data
# Determine diff methylated loci

#For execution on a local, multicore CPU with excess RAM we recommend calling
#options(mc.cores = parallel::detectCores()).
#To avoid recompilation of unchanged Stan programs, we recommend calling
#rstan_options(auto_write = TRUE)

load("RData/meth.RData")
library(rethinking)
library(tidyverse)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

treatment <- slot(meth, "treatment")
chr1 <- meth[meth$chr == "chr1", ] 
# number the CpG loci in the data 
ncpg <- nrow(chr1)
chr1 <- cbind(chr1, 1:ncpg)
colnames(chr1)[20] <- "loc"

# reshape the data so that each row is the coverage and #C of one sample
reshape_data <- function(meth_data) {
  s1 <- meth_data[, c("chr", "start", "end", "coverage1", "numCs1", "loc")]
  colnames(s1)[c(4, 5)] <- c("coverage", "numCs")
  s1$treat <- treatment[1]
  s1$sample <- 1
  s2 <- meth_data[, c("chr", "start", "end", "coverage2", "numCs2", "loc")]
  colnames(s2)[c(4, 5)] <- c("coverage", "numCs")
  s2$treat <- treatment[2]
  s2$sample <- 2
  s3 <- meth_data[, c("chr", "start", "end", "coverage3", "numCs3", "loc")]
  colnames(s3)[c(4, 5)] <- c("coverage", "numCs")
  s3$treat <- treatment[3]
  s3$sample <- 3
  s4 <- meth_data[, c("chr", "start", "end", "coverage4", "numCs4", "loc")]
  colnames(s4)[c(4, 5)] <- c("coverage", "numCs")
  s4$treat <- treatment[4]
  s4$sample <- 4
  s5 <- meth_data[, c("chr", "start", "end", "coverage5", "numCs5", "loc")]
  colnames(s5)[c(4, 5)] <- c("coverage", "numCs")
  s5$treat <- treatment[5]
  s5$sample <- 5
  
  return(rbind(s1, s2, s3, s4, s5))
}



# Bayesian model for differentially methylated sites
m_DMC <- ulam(
  alist(
    C ~ dbinom(N, p),
    logit(p) <- a[L, T],
    
    matrix[L, T]: a ~ normal(0, 1.5)
  ), data=dat, chains=4, cores=4
)

# Try to model methylation at loc using a beta-binomial distribution.
# The power was lower because of the dispersion
m_DMC_beta <- ulam(alist(
  C ~ dbetabinom(N, pbar, theta),
  logit(pbar) <- a[L, T],
  matrix[L, T]: a ~ normal(0, 1.5),
  transpars> theta <<- phi + 10,
  phi ~ dexp(0.1)
), data=dat, chains=4, cores=4)

precis(m_DMC, depth=3)
m_DMC_code <- stancode(m_DMC)
m_DMC_model <- stan_model(model_code=m_DMC_code, verbose=FALSE)
res <- sampling(m_DMC_model, data=dat, cores=4, verbose=FALSE) # data must be the same dimension ?!
x <- precis(res, depth=3)

# sample posterior and determine differentially methylated C's (DMCs)
calc_DMCs <- function(model, level=0.99, diff=0.25) {
  post <- extract.samples(model)
  # count the samples with higher a values in treatment 1 vs. treatment 2
  x <- apply(post$a[, , 1]-post$a[, , 2], 2, function (x) sum(ifelse(x > 0, 1, 0)))
  # convert it into percentage
  x <- x / length(post$a[,1, 1])
  # check if the percentage is more significant than the confidence level
  y <- sapply(x, function(i) i > level | i < (1-level))
  # check if the diff in mean methyl level > 0.3
  m1 <- apply(post$a[, , 1], 2, mean)
  m1 <- sapply(m1, inv_logit)
  m2 <- apply(post$a[, , 2], 2, mean)
  m2 <- sapply(m2, inv_logit)
  z <- sapply(m1-m2, function(i) abs(i) > diff)
  y <- y & z
  return (data.frame(m1=m1, m2=m2, dmc=y))
}

# Process large arrays of methyl data by dividing them into blocks of size 1000
m_data <- getData(meth)
nr <- nrow(m_data)
m_data <- cbind(m_data, 1:nr)
colnames(m_data)[20] <- "loc"
dmc <- data.frame(matrix(ncol=3, nrow=0))
colnames(dmc) <- c("m1", "m2", "dmc")
batch = 10
for (i in 1:ceiling(nr/batch)) {
  print(paste0("processing block ", i))  
  d <- m_data[((i-1)*batch+1):(i*batch), ]
  print(dim(d))
  nval = min(batch, nr-(i-1)*batch);    #number of valid entries in the block
  if (nval < batch) {
    d[is.na(d)] <- 0        # convert NA into 0
  }
  # reshape the data
  d <- reshape_data(d)
  # set up list for stan sampling
  dat <- list(
    N = d$coverage,     # number of coverage
    C = d$numCs,        # number of Cs
    L = rep(1:batch, 5),          # loci number
    S = as.integer(d$sample),  # sample number, 1 - 5
    T = ifelse(d$treat == 0, 1, 2)  # treatment 1: EH1516, 2: E217
  )
  # run the sampling
  res <- sampling(m_DMC_model, data=dat, cores=4, verbose=FALSE) 
  # Determine DMCs
  x <- calc_DMCs(res)
  dmc <- rbind(dmc, x[1:nval, ])
}
dmc_data <- cbind(m_data, dmc)

# par(mfrow=c(1,1))
# d1 <- list(
#   N = rep(100, 1000),
#   L = 1:1000,
#   T = rep(1, 1000)
# )
# 
# sim_d1 <- sim(m_DMC, data=d1, n=2000)
# d1$T = rep(2, 1000)
# sim_d2 <- sim(m_DMC, data=d1, n=2000)
# k = 314
# dens(sim_d1[, k], col=4)
# dens(sim_d2[, k], col=2, add=TRUE)
# dens(sim_d1[, k]-sim_d2[, k])
# dens(inv_logit(post$a[, k, 2]), col=2, add=TRUE)
# dens(inv_logit(post$a[, k, 1]) - inv_logit(post$a[, k, 2]))
# dens(inv_logit(post$a[, k, 1]), col=4)
# dens(inv_logit(post$a[, k, 2]), col=2, add=TRUE)
# dens(inv_logit(post$a[, k, 1]) - inv_logit(post$a[, k, 2]))

# Merge the meth stat partitions 
load("./part_data/dmc_p_1.RData")
dmc_data_all <- dmc_data
for (i in 2:15) {
  file_name <- paste0("./part_data/dmc_p_",i,".RData")
  load(file_name)
  dmc_data_all <- rbind(dmc_data_all, dmc_data)
}
# set confidence interval for each loc
#dmc_data_all$pval <- rowMin(cbind(dmc_data_all$pval, 1-dmc_data_all$pval))

# split dmc_data_all into a named list index by chr
start_idx <- list(1)
cur_chr <- "chr1"
# find the start index of each chr
for(i in 2:nrow(dmc_data_all)) {
  chr_i <- dmc_data_all[i, "chr"]
  if(chr_i != cur_chr) {
    cur_chr <- chr_i
    start_idx <- append(start_idx, i)
  }
}
start_idx <- append(start_idx, nrow(dmc_data_all)+1)
dmc_data_list <- list()
for(i in 1:(length(start_idx)-1)) {
  start <- start_idx[[i]]
  end <- start_idx[[i+1]] - 1
  chr <- dmc_data_all[start, "chr"]
  dmc_data_list[[chr]] <- dmc_data_all[start:end, ]
}
save(dmc_data_list, file="dmc_list.RData")

# DMCS are selected as site withd diff > 0.25 and pval = 0.01
dmcs <- dmc_data_all %>% filter (abs(m1-m2) > 0.25 & pval < 0.01)
nrow(dmcs %>% filter (m1 < m2))
nrow(dmcs %>% filter (m1 > m2))


# Compare the DMCs with those of MOABS
moabs.DMC <- readMOABS_DMR_bedfile("./data/dmc_M2_EH1516.G.bed_vs_EH217.G.bed.bed")
# only consider sites with coverage >=3 in all samples
moabs.DMC <- as.data.frame(moabs.DMC)
colnames(moabs.DMC)[1] = "chr"
moabs.dmcs <- merge(dmc_data_all, moabs.DMC, by=c("chr", "start", "end"), sort=FALSE)
moabs.dmcs <- moabs.dmcs[, 1:24]
nrow(moabs.dmcs %>% filter (m1 < m2))
nrow(moabs.dmcs %>% filter (m1 > m2))
save(dmcs, moabs.dmcs, file="RData/moabs.dmcs.RData")
# count the number of overlapped DMCs
sum(countOverlaps(as(moabs.dmcs, "GRanges"), as(dmcs, "GRanges"))) #/ nrow(moabs.dmcs)

# The DMCs in MOABS but not in this analysis
moabs.Uniq <- moabs.dmcs[countOverlaps(as(moabs.dmcs, "GRanges"), as(dmcs, "GRanges")) == 0, ]
dmcs.Uniq <- dmcs[countOverlaps(as(dmcs, "GRanges"), as(moabs.dmcs, "GRanges")) == 0, ]

# Set the dmc flag to TRUE if it's in moabs.dmc
dmc_data_all$dmc <- countOverlaps(as(dmc_data_all, "GRanges"), as(moabs.dmcs, "GRanges")) ==1 
# recalculate m1 & m2 using coverage and Cs
dmc_data_all$m1 <- (dmc_data_all$numCs1+dmc_data_all$numCs2) / (dmc_data_all$coverage1+dmc_data_all$coverage2)
dmc_data_all$m2 <- (dmc_data_all$numCs3+dmc_data_all$numCs4+dmc_data_all$numCs5) / 
  (dmc_data_all$coverage3+dmc_data_all$coverage4+dmc_data_all$coverage5)
save(dmc_data_all, file="RData/dmc.RData")

# DMCs in methylKit analysis
load("RData/methtylKit.DMC.RData")
sum(countOverlaps(as(methylKit.DMC, "GRanges"), as(moabs.dmcs, "GRanges"))) / nrow(methylKit.DMC)

# The DMCs in DSS analysis
load("RData/dss.DMC.RData") 


