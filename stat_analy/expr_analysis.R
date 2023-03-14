# Analyze differentially expressted genes (DEGs) using Bayesian statistical models

# read gene expression data
expr_data <- read.table("data/EHX_1516_v_217_Jan2019.txt", sep="\t", header=TRUE )

expr_data = expr_data[rowSums(expr_data > 0) >= 3,]  # at least 3 columns Counts  > 0
#expr_data <- expr_data[1:100, ]
# reshape the gene expr data
expr_data$gid <- 1:nrow(expr_data)

strain <- c(rep(0, 3), rep(1, 3))  # 1: EH1516, 2: EH217
expr_reshape <- expr_data[, c(1, 2, 8)]  # the first sample
colnames(expr_reshape) <- c("gid", "expr", "id")
expr_reshape$sample <- 1
expr_reshape$strain <- strain[1]
for (i in 2:6) {  # sample 2-6 in the data
  d_i <- expr_data[, c(1, i+1, 8)]
  colnames(d_i) <- c("gid", "expr", "id")
  d_i$sample <- i
  d_i$strain <- strain[i]
  expr_reshape <- rbind(expr_reshape, d_i)
}

# set up data list for MCMC
dat <- list(
  G = expr_reshape$id,           # gene id from 1:ng
  E = expr_reshape$expr,         # expression count
  S = expr_reshape$sample,       # sample id 1:6
  T = expr_reshape$strain,        # treatment/strain 1&2
  e_bar = rep(log(median(expr_reshape$expr)), nrow(expr_reshape)),
  ng = max(expr_reshape$id)
)

# model gene expression 
# negative-binomial model
m_Expr <- ulam (
  alist (
    E ~ dgampois(lambda, phi),
    # f[S]: log sample factor, e[G] log express of gene in treatment 1, b[G] change in treatment 2
    log(lambda) <- e_bar + f[S] + e[G] + b[G]*T,      
  
    vector[6]: f ~ normal(0, 0.1),
    vector[ng]: e ~ normal(0, 3),
    vector[ng]: b ~ normal(0, 3),
    phi ~ dexp(1)
  ), data = dat, chains=4, cores=4
)

params <- precis(m_Expr, depth=2, prob=0.95)

# Determine DEGs
f <- params[1:6, ]
e <- params[7:(7+dat$ng-1), ]
b <- params[(7+dat$ng):(6+2*dat$ng), ]
DE <- b$`2.5%` > 0 | b$`97.5%` <0
deg_data <- cbind(expr_data, e)
deg_data <- cbind(deg_data, b)
deg_data <- cbind(deg_data, DE)
colnames(deg_data)[9] <- "expr_base"
colnames(deg_data)[15] <- "logfc"
deg_data <- deg_data[rowSums(deg_data[, 2:7] > 0) >= 3,] 
save(deg_data, file="deg.RData")
