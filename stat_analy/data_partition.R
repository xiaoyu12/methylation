# Partition the large data into 10 (or more) partitions to run in parallel
library(methylKit)

load("RData/meth.RData")
meth_data <- getData(meth)
nrow <- nrow(meth_data)

npar <- 16
psize <- ceiling(nrow / npar)
par <- seq(0, nrow, psize)
par <- c(par, nrow)

for(i in 1:npar) {
  file_name <- paste0("part_data/meth_p_", i, ".RData")
  d <- meth_data[(par[i]+1):par[i+1], ]
  save(d, file=file_name)
  #write.csv(d, file=file_name, quote=FALSE)
}

# load and partition CHG data
load("RData/methObjs.RData")
meth.CHG.cov3 <- methylKit::unite(filterByCoverage(methCHG, lo.count = 3))
head(meth.CHG.cov3)
meth_CHG <- getData(meth.CHG.cov3)

nrow <- nrow(meth_CHG)

npar <- 16
psize <- ceiling(nrow / npar)
par <- seq(0, nrow, psize)
par <- c(par, nrow)

for(i in 1:npar) {
  file_name <- paste0("part_data/CHG_p_", i, ".RData")
  d <- meth_CHG[(par[i]+1):par[i+1], ]
  save(d, file=file_name)
  #write.csv(d, file=file_name, quote=FALSE)
}
