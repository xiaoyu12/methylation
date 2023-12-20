library("methylKit")

# read Bismark CpG and CHG files
#base_dir = "I:/data/methylKit/"
read_Bismark_CpG_CHG <- function(base_dir="E:/xiaoyu/Documents/methylKit/") { 
  # Those files are output the processBismarkAln() function of the methylKit package
  list.CpG <- list(paste0(base_dir, "EH1516C_CpG.txt"),
                 paste0(base_dir, "EH1516D_CpG.txt"),
                 paste0(base_dir, "EH217A_CpG.txt"),
                 paste0(base_dir, "EH217B_CpG.txt"),
                 paste0(base_dir, "EH217C_CpG.txt"))
  list.CHG <- list(paste0(base_dir, "EH1516C_CHG.txt"),
                 paste0(base_dir, "EH1516D_CHG.txt"),
                 paste0(base_dir, "EH217A_CHG.txt"),
                 paste0(base_dir, "EH217B_CHG.txt"),
                 paste0(base_dir, "EH217C_CHG.txt"))

  sample.ids = list("EH1516B", "EH1516C", "EH217A", "EH217B", "EH217C")
  methCpG <- methRead(list.CpG, sample.id = sample.ids, assembly = "ehux", treatment = c(0, 0, 1, 1, 1),
                    context="CpG", pipeline = "bismark", mincov = 3)

  methCHG <- methRead(list.CHG, sample.id = sample.ids, assembly = "ehux", treatment = c(0, 0, 1, 1, 1),
                    context="CHG", pipeline = "bismark", mincov = 3)

  #save(methCpG, methCHG, file="methObjs.RData")
  meth_list = list("methCpG" = methCpG, "methCHG" = methCHG)
  #meth_list[[1]] = methCpG
  #meth_list[[2]] = methCHG
  return (meth_list)
}

# read Bismark coverage files
read_Bismark_coverage <- function(base_dir = "./") {
  # List of CpG coverage files
  file.list <- list(paste0(base_dir, "EH1516C.merged_CpG_evidence.cov"),
                    paste0(base_dir, "EH1516D.merged_CpG_evidence.cov"),
                    paste0(base_dir, "EH217A.merged_CpG_evidence.cov"),
                    paste0(base_dir, "EH217B.merged_CpG_evidence.cov"),
                    paste0(base_dir, "EH217C.merged_CpG_evidence.cov"))
  
  sample.ids = c("EH1516B", "EH1516C", "EH217A", "EH217B", "EH217C")
  # Read methylation data
  methobj <- methRead(file.list, sample.id = as.list(sample.ids),
                      assembly = "ehux", pipeline = "bismarkCoverage", treatment = c(0, 0, 1, 1, 1),
                      context = "CpG", mincov = 3, header = FALSE)
  return(methobj)
}