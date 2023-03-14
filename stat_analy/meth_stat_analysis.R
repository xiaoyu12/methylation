# Script to perform Bayesian stat analysis on E. hux methylation data
# Determine diff methylated loci

#For execution on a local, multicore CPU with excess RAM we recommend calling
#options(mc.cores = parallel::detectCores()).
#To avoid recompilation of unchanged Stan programs, we recommend calling
#rstan_options(auto_write = TRUE)

load("meth.RData")
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
precis(m_DMC, depth=3)
m_DMC_code <- stancode(m_DMC)
m_DMC_model <- stan_model(model_code=m_DMC_code, verbose=FALSE)
res <- sampling(m_DMC_model, data=dat, cores=4, verbose=FALSE) # data must be the same dimension ?!
x <- precis(res, depth=3)

# sample posterior and determine differentially methylated C's (DMCs)
calc_DMCs <- function(model, level=0.95) {
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
  z <- sapply(m1-m2, function(i) abs(i) > 0.3)
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
for (i in 1:ceiling(nr/1000)) {
  print(paste0("processing block ", i))  
  d <- m_data[(i*1000-999):(i*1000), ]
  print(dim(d))
  nval = min(1000, nr-(i-1)*1000);    #number of valid entries in the block
  if (nval < 1000) {
    d[is.na(d)] <- 0        # convert NA into 0
  }
  # reshape the data
  d <- reshape_data(d)
  # set up list for stan sampling
  dat <- list(
    N = d$coverage,     # number of coverage
    C = d$numCs,        # number of Cs
    L = rep(1:1000, 5),          # loci number
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
load("dmc_p_1.RData")
dmc_data_all <- dmc_data
for (i in 2:15) {
  file_name <- paste0("dmc_p_",i,".RData")
  load(file_name)
  dmc_data_all <- rbind(dmc_data_all, dmc_data)
}
save(dmc_data_all, file="dmc.RData")
# split dmc_data_all into a named list index by chr
start_idx <- list(1)
cur_chr <- "chr1"
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
