# Combine the gene expression data from EHX1516 vs EH217 strains and EH217 strain
# w/o the spike of sodium carbonate to study the co-expression of genes
# We are particularly interested in the genes over-expressed in EH217 and with spike

# Data preparation
datExprA2 <- expr_data      # expression data of EH1516 and EH217 strains
datExprA1 <- read.table("data/Ehux_genbank_filtered_w_ids_r2.matrix", 
                        header=TRUE, row.names = 1, sep="\t")
# only consider 9 and 9S columns in datExprA1
datExprA1 <- datExprA1[, 5:8]

# Remove rows with all zeros @TODO: review those conditions
isexpr <- (rowSums(datExprA1 >=15) >=2) & (rowSums(datExprA2 >= 15) >= 3)
datExprA1 <- datExprA1[isexpr, ]
datExprA2 <- datExprA2[isexpr, ]

# reshape the data so that there is one gene & sample per line
d <- cbind(datExprA2, datExprA1)

nsample = ncol(d) - 1
strain <- c(rep(0, 3), rep(1, 7))
spike <- c(rep(0, 8), rep(1, 2))
# the first sample
d_reshaped <- d[, c(1, 2)]
colnames(d_reshaped) <- c("gid", "expr")
d_reshaped$sample <- 1
d_reshaped$strain <- strain[1]
d_reshaped$spike <- spike[1]
for (i in 2:nsample) {
  d_i <- d[, c(1, i+1)]
  colnames(d_i) <- c("gid", "expr")
  d_i$sample <- i
  d_i$strain <- strain[i]
  d_i$spike <- spike[i]
  d_reshaped <- rbind(d_reshaped, d_i)
}
d_reshaped$id <- rep(1:nrow(d), nsample)

# Set up data list for MCMC
e_bar <- log(median(d_reshaped$expr))
dat <- list (
  G = d_reshaped$id,           # gene id from 1:ng
  E = d_reshaped$expr,         # expression count
  S = d_reshaped$sample,       # sample id 
  T = d_reshaped$strain,        # strain 1&2
  SK = d_reshaped$spike,        # treatment of NaCO3 spike
  e_bar = rep(e_bar, nrow(d_reshaped)),
  ng = max(d_reshaped$id),
  ns = max(d_reshaped$sample)
)

m_co_expr <- ulam(
  alist(
    E ~ dgampois(lambda, phi),
    # e[G]: base experssion of gene (in EH1516)
    # b[G]: effect of strains on gene expression
    # c[G]: effect of spike on gene expression
    log(lambda) <- e_bar + f[S] + e[G] + b[G]*T + c[G]*SK,
    
    vector[ns]: f ~ normal(0, 0.1),
    vector[ng]: e ~ normal(0, 3),
    vector[ng]: b ~ normal(0, 3),
    vector[ng]: c ~ normal(0, 3),
    phi ~ dexp(1)
  ), data = dat, chains=4, cores=4
)
save(m_co_expr, file="m_co_expr.RData")
post <- extract.samples(m_co_expr)
params <- precis(m_co_expr, prob=0.95)
ng <- nrow(d)
p.expr <- params[11:(ng+10),]
p.strain <- params[(ng+11):(2*ng+10), ]
p.spike <- params[(2*ng+11):(3*ng+10) ,]
# genes over-expressed both in EH217 strain and with NaCO3 spike
both_up <- p.strain$`2.5%` > 0 & p.spike$`2.5%` > 0
both_up <- d[shared_up, ]

# A modified model with variable dispersion for genes
m_co_expr_vs <- ulam(
  alist(
    E ~ dgampois(lambda, phi),
    # e[G]: base experssion of gene (in EH1516)
    # b[G]: effect of strains on gene expression
    # c[G]: effect of spike on gene expression
    log(lambda) <- f[S] + e[G] + b[G]*T + c[G]*SK,
    phi <- dis[G],
    vector[ns]: f ~ normal(0, 0.1),
    vector[ng]: e ~ normal(0, 3),
    vector[ng]: b ~ normal(0, 3),
    vector[ng]: c ~ normal(0, 3),
    vector[ng]: dis ~ dexp(1)
  ), data = dat, chains=4, cores=4
)

# Compute the scale factor using median of ratios method used by DESeq2
scale_factor <- function(matrix) {
  #print(head(matrix))
  geo_mean <- apply(matrix, 1, function(x) exp(mean(log(x))))
  not_zero <- geo_mean != 0
  geo_mean <- geo_mean[not_zero]
  matrix <- matrix[not_zero, ]
  # divide each col by the geo_mean vector
  matrix <- sweep(matrix, 1, geo_mean, FUN="/")
  # get median of ratios as the scale factor
  sf <- apply(matrix, 2, median)
  return(sf)
} 

matrix <- d[2:11]
sf <- scale_factor(d[2:11])
