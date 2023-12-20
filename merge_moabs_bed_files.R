##############################################################
# This script merges two MOABS mcall bed files into a 
# combined bed file, which may be used for mcomp DMR analysis
##############################################################
args = commandArgs(trailingOnly = TRUE)
print(paste(args, collapse = " "))

# Input bed file names
file1 = args[1]
file2 = args[2]

# Read input bed files
print(paste("Read file", args[1]))
bed1 <- read.table(file1, sep="\t", header = TRUE, comment.char = "", strip.white = TRUE, blank.lines.skip = TRUE)
print(paste("Read file", args[2]))
bed2 <- read.table(file2, sep="\t", header = TRUE, comment.char = "", strip.white = TRUE, blank.lines.skip = TRUE)

bed <- merge(bed1, bed2, by = colnames(bed1)[1:3], all.x= FALSE, all.y=FALSE, sort=FALSE)
bed$totalC.x <- bed$totalC.x + bed$totalC.y
bed$methC.x <- bed$methC.x + bed$methC.y
bed$ratio.x <- bed$methC.x / bed$totalC.x
bed$totalC.1.x <- bed$totalC.1.x + bed$totalC.1.y
bed$methC.1.x <- bed$methC.1.x + bed$methC.1.y
bed$totalC.2.x <- bed$totalC.2.x + bed$totalC.2.y
bed$methC.2.x <- bed$methC.2.x + bed$methC.2.y
bed <- bed[, 1:15]
colnames(bed) <- colnames(bed1)

print(paste("Write to file", args[3]))
write.table(bed, file=args[3], sep = '\t', col.names = TRUE, row.names = FALSE, quote=FALSE)