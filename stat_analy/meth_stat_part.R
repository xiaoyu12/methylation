library(rethinking)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  stop("Data file name is required")
} else if(length(args) == 1) {
  # output file
  args[2] = "out.RData"
}

# Load the input data d
load(args[1])

treatment = c(0, 0, 1, 1, 1)

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

# sample posterior and determine differentially methylated C's (DMCs)
calc_DMCs <- function(model, level=0.99) {
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
  return (data.frame(m1=m1, m2=m2, pval=x, dmc=y))
}

m_data <- d
nr <- nrow(m_data)
m_data <- cbind(m_data, 1:nr)
colnames(m_data)[20] <- "loc"

# Bayesian model for differentially methylated sites
# process the first 1000 CpG and build the model
batch = 1000    # batch size set to be 1000
d <- m_data[1:batch, ]
d <- reshape_data(d)
dat <- list(
  N = d$coverage,     # number of coverage
  C = d$numCs,        # number of Cs
  L = rep(1:batch, 5),          # loci number
  S = as.integer(d$sample),  # sample number, 1 - 5
  T = ifelse(d$treat == 0, 1, 2)  # treatment 1: EH1516, 2: E217
)
m_DMC <- ulam(
  alist(
    C ~ dbinom(N, p),
    logit(p) <- a[L, T],
    
    matrix[L, T]: a ~ normal(0, 1.5)
  ), data=dat, chains=4, cores=4
)
dmc <- calc_DMCs(m_DMC)
# precis(m_DMC, depth=3)
m_DMC_code <- stancode(m_DMC)
m_DMC_model <- stan_model(model_code=m_DMC_code, verbose=FALSE)

#dmc <- data.frame(matrix(ncol=3, nrow=0))
#colnames(dmc) <- c("m1", "m2", "dmc")
for (i in 2:ceiling(nr/batch)) {
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
save(dmc_data, file=args[2])