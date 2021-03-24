library("methylKit")

base_dir = "E:/xiaoyu/Documents/methylKit/"
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

sample.ids = list("EH1516C", "EH1516D", "EH217A", "EH217B", "EH217C")
methCpG <- methRead(list.CpG, sample.id = sample.ids, assembly = "ehux", treatment = c(0, 0, 1, 1, 1),
                    context="CpG", pipeline = "bismark", mincov = 3)

methCHG <- methRead(list.CHG, sample.id = sample.ids, assembly = "ehux", treatment = c(0, 0, 1, 1, 1),
                    context="CHG", pipeline = "bismark", mincov = 3)

save(methCpG, methCHG, file="methObjs.RData")
