library(tidyverse)

# Consider the genes enriched in GO:0030286 (dynein complex)
print(GO_dynein)

# Consider genes enriched in GO:0005929 (cilium)
GO_cilium <- enriched_genes_by_GO_term(GO_DE_up, "GO:0005929")

#read DMGs from .csv file
DMGs <- read.csv("../DMGs.csv", row.names = 1)

#read DMPs from .csv file
DMPs <- read.csv("../DMPs.csv", row.names = 1)

# overlap enriched genes with DMGs and DMPs
DMGs_GO_dynein = DMGs %>% filter(ID %in% GO_dynein)
DMPs_GO_dynein = DMPs %>% filter(ID %in% GO_dynein)

# DE shared up genes
DE_shared_up <- read.csv("../DE_shared_up.csv", row.names = 1)
DE_shared_up$ID <- rownames(DE_shared_up)
DE_shared_up_dynein <- DE_shared_up %>% filter(ID %in% GO_dynein)
DE_shared_up_cilium <- DE_shared_up %>% filter (ID %in% GO_cilium)

Expr_9_n_9S <- read.table("../Ehux_genbank_filtered_w_ids_r2.matrix", row.names = 1, header=1)
Expr_E1516_E217 <- read.table("../EHX_1516_v_217_Jan2019.txt", row.names = 1, header = 1)
Expr_9_n_9S[GO_dynein, 5:8]
Expr_E1516_E217[GO_dynein, ]

Expr_9_n_9S[GO_microtubule, 5:8]
Expr_E1516_E217[GO_microtubule, ]

Expr_9_n_9S[GO_cilium, 5:8]
Expr_E1516_E217[GO_cilium, ]
