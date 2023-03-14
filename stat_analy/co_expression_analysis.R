##############################################################################
# Co-expression analysis of E. hux EH217 vs EH1516 strains and 9S vs 9 
# growth conditions
##############################################################################
library('Rcpp')
library("WGCNA")
library("DESeq2")
library("limma")
library("flashClust")
library("tidyverse")

#library(doParallel)
#registerDoParallel(cores=16)

options(stringsAsFactors = FALSE);
#allowWGCNAThreads(nThreads = 16)
enableWGCNAThreads()

nSets = 2
setLabels = c("9S-vs-9", "EH217-vs-EH1516")
# Form multi-set expression data:
multiExpr = vector(mode = "list", length = nSets)

datExprA1 <- read.table("data/Ehux_genbank_filtered_w_ids_r2.matrix", header=TRUE, row.names = 1, sep="\t")
datExprA2 <- read.table("data/EHX_1516_v_217_Jan2019.txt", header=TRUE,  sep="\t")

# only consider 9 and 9S columns in datExprA1
datExprA1 <- datExprA1[, 5:8]
colnames(datExprA1) <- c("EH217_9.1", "EH217_9.2", "EH217_9S.1", "EH217_9S.2")

# Remove rows with all zeros
isexpr <- (rowSums(datExprA1 >=15) >=2) & (rowSums(datExprA2 >= 15) >= 3)
datExprA1 <- datExprA1[isexpr, ]
datExprA2 <- datExprA2[isexpr, ]

d <- data.frame(datExprA2, datExprA1)
matrix <- d[, 2:11]
#Normalization
samples <- data.frame(sample = colnames(matrix))
samples$strain <- as.factor(c(rep("EH1516", 3), rep("EH217", 7)))
samples$condition <- as.factor(c(rep("new", 6), rep("C9", 2), rep("C9S", 2)))
rownames(samples) <- samples$sample
dds <- DESeqDataSetFromMatrix(countData=matrix, colData=samples, ~strain + condition)  # Creating a DESeqDataSet object
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)
dds <- estimateSizeFactors(dds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)

dds <- DESeq(dds)
# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) #%>% t() # Transpose this data
normalized_counts <- log.norm.counts

keep = rank(-rowMeans(normalized_counts)) <= 5000
top_counts = normalized_counts[keep,]

sft <- pickSoftThreshold(t(top_counts),
                         dataIsExpr = TRUE,
                         powerVector = c(seq(1, 10, by=1), seq(10, 20, by=2), seq(22, 60, by=4)),
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

# Choose softpower to be 18 according to the plot 
softPower = 18

bwnet <- blockwiseModules(t(normalized_counts),
                          maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = softPower, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          minModuleSize = 100,
                          randomSeed = 12345, # there's some randomness associated with this calculation, so we should set a seed
                          verbose = 3,
)
# Write main WGCNA results object to file
readr::write_rds(bwnet,
                 file = file.path("results", "Ehux_wgcna_results.RDS")
)

# Another way to determine modules
adjacency = adjacency(t(normalized_counts),power=softPower,type="signed");
#diag(adjacency)=0
TOM = TOMsimilarity(adjacency, TOMType="signed")
dissTOM = 1 - TOM
geneTree = flashClust(as.dist(dissTOM), method="ward")

moduleLabels = cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = 2,                             
                             #cutHeight = 30,
                             minClusterSize = 100,
                             pamRespectsDendro = FALSE);
moduleColors = labels2colors(moduleLabels)
table(moduleColors)
modules <- cbind(d, moduleLabels)
modules <- cbind(modules, moduleColors)
modules[, 1:10] <- NULL
modules$gene <- rownames(modules)

# Get module eigengenes
#module_eigengenes <- bwnet$MEs
module_eigengenes <- moduleEigengenes(t(normalized_counts), colors=moduleColors)
module_eigengenes <- module_eigengenes$eigengenes

# Print out a preview
head(module_eigengenes)

#samples$condition = as.factor(samples$condition)
mod <- model.matrix(~samples$condition)
fit <- lmFit(t(module_eigengenes), design=mod)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")
head(stats_df)

#ME5 and ME16 seems to be the most interesting
# "101141" in ME1 doesn't make sense
module_df <- module_eigengenes %>%
  tibble::rownames_to_column("accession_code") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(samples %>%
                      dplyr::select(sample, condition),
                    by = c("accession_code" = "sample")
  )

ggplot(
  module_df,
  aes(
    x = condition,
    y = MEturquoise,
    color = condition
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()

# What genes are a part of module 16?
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key %>%
  dplyr::filter(module == "ME16")

# Plot mean profile of modules
labeledHeatmap(Matrix = t(module_eigengenes), xLabels = rownames(module_eigengenes), 
               yLabels = colnames(module_eigengenes),
               colors= blueWhiteRed(50),
               cex.text = 0.5)

# Read shared up genes
shared_up <- read.table("data/DE_shared_up.txt")
#shared_up <- gene_module_key %>% filter(gene %in% shared_up$V1)
shared_up <- modules %>% filter(gene %in% shared_up$V1)
table(shared_up$moduleColors)

