# Calculate the CHG DMC effects on gene expressions
source ("../meth_analysis_funcs.R")

# first load CHG DMCs
moabs.HG.DMC <- readMOABS_DMR_bedfile("../dmc_M2_EH1516.HG.bed_vs_EH217.HG.bed.bed")
dmc_CHG <- as.data.frame(moabs.HG.DMC)
dmc_CHG_P <- dmc_CHG %>% filter (meth.diff > 0)  # hyper-meth in EH217
dmc_CHG_N <- dmc_CHG %>% filter (meth.diff < 0)  # hypo-meth in EH217

# count DMCs in one gene region

count_dmcs_in_region <- function(gr, dmc_P, dmc_N) {
  nr <- length(gr)
  dmc_P_gr <- as(dmc_P, "GRanges")
  dmc_N_gr <- as(dmc_N, "GRanges")
  m_count = data.frame(dmcs = rep(0, nr), dmc_neg = rep(0, nr), dmc_pos = rep(0, nr))
  m_count$dmc_pos <- countOverlaps(gr, dmc_P_gr, ignore.strand=TRUE)
  m_count$dmc_neg <- countOverlaps(gr, dmc_N_gr, ignore.strand=TRUE)
  m_count$dmcs <- m_count$dmc_pos + m_count$dmc_neg
  
  return(m_count)
}


dmc_CHG_count <- count_dmcs_in_region(gene_gr, dmc_CHG_P, dmc_CHG_N)































