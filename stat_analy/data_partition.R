# Partition the large data into 10 (or more) partitions to run in parallel
load("meth.RData")
meth_data <- getData(meth)
nrow <- nrow(meth_data)

npar <- 15
psize <- ceiling(nrow / npar)
par <- seq(0, nrow, psize)
par <- c(par, nrow)

for(i in 1:npar) {
  file_name <- paste0("meth_p_", i, ".RData")
  d <- meth_data[(par[i]+1):par[i+1], ]
  save(d, file=file_name)
  #write.csv(d, file=file_name, quote=FALSE)
}
