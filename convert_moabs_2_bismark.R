######################################################
# This script converts MOABS bed file into bismark
# merged CpG format
######################################################
args = commandArgs(trailingOnly = TRUE)
print(paste(args, collapse = " "))


# Read input bed file
print(paste("Reading file", args[1]))
bed <- read.table(args[1], sep = "\t", header=TRUE, comment.char = "", strip.white = TRUE, blank.lines.skip = TRUE)
bed$start <- bed$start + 1
bed$ratio <- bed$ratio * 100.0
bed$nonMethC <- bed$totalC - bed$methC
bed <- bed[, c(1, 2, 3, 4, 6, 9)]

# Write to output file
print(paste("Writing to file", args[2]))
write.table(format(bed, scientific=FALSE), file = args[2], sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
