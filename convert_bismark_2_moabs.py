'''
This program converts Bismark CpG_report file to moabs (mcall) .bam.G.bed format,
in order to use mcomp for DMR analysis.

The reason for this conversion is that moabs default alignment program bsmap seems to
produce inconsistent results.

copyright: 2019 Xiaoyu Zhang
'''
from optparse import OptionParser
import pandas as pd

def make_a_bed_row(chr, start, end, methC_p, totalC_p, methC_n, totalC_n, context):
    totalC = totalC_p + totalC_n
    assert totalC > 0
    methC = methC_p + methC_n
    ratio = float(methC) / totalC
    strand = '-'
    if totalC_p > 0 and totalC_n > 0:
        strand = 'B'
    elif totalC_p > 0:
        strand = '+'
    bedrow = [chr, start, end, ratio, totalC, methC, strand, 'G', '+', totalC_p, methC_p, '-',
              totalC_n, methC_n, context]
    return bedrow


parser = OptionParser(description="This program converts Bismark coverage file to moabs (mcall) .bam.G.bed format")
parser.add_option("-i", "--input", dest="bismark_file", help="Bismark coverage2cytosine output CpG_report file")
parser.add_option("-o", "--output", dest="bed_file", help="Output Bed file in moabs (mcall) format")

(options, args) = parser.parse_args()

# Read input file
print("Reading input file", options.bismark_file)
cpg_table = pd.read_table(options.bismark_file, header=None)

#print(bed_table)

# Iterate through the cpg table
nrows = cpg_table.shape[0]
i = 0
bedrows = []
while i < nrows:
    row = cpg_table.iloc[i]
    totalC_p = totalC_n = 0
    methC_p = methC_n = 0
    context = row[6]
    if row[2] == '+':
        totalC_p = row[3] + row[4]
        methC_p = row[3]
        start = row[1]
        end = row[1]+2
        if i+1 < nrows:
            next = cpg_table.iloc[i+1]
            # the two rows form a CpG pair
            if next[2] == '-' and next[0] == row[0] and next[1] == row[1]+1:
                totalC_n = next[3] + next[4]
                methC_n = next[3]
                i = i+1
    else:  # Only - strand has data
        totalC_n = row[3] + row[4]
        methC_n = row[3]
        start = row[1]-1
        end = row[1]+1

    if totalC_p + totalC_n > 0:
        # bedrow = {"chrom" : row[0],  "start" : row[1]-1,   "end" : next[1],
        #          "ratio" : ratio, "totalC" : totalC, "methC" : methC,
        #          "strand": strand, "next" : 'G', "Plus" : '+', "totalC" : totalC_p,
        #          "methC" : row[3], "Minus" : '-', "totalC" : totalC_n, "methC" : next[3], "localSeq" : row[6]}
        bedrow = make_a_bed_row(row[0], start, end, methC_p, totalC_p, methC_n, totalC_n, context)
        bedrows.append(bedrow)
            # print(bedrow)
    i = i+1     # move to next row
    if (i % 100001) == 0:
        print("Processing row {} out of {} rows".format(i, nrows))

# Moabs bed file has 15 columns
bed_table = pd.DataFrame(bedrows, columns=["chrom",  "start",   "end",  "ratio", "totalC", "methC", "strand", "next", "Plus", "totalC", "methC", "Minus", "totalC", "methC", "localSeq"])

bed_table.to_csv(options.bed_file, sep='\t', header=True, index=False)
