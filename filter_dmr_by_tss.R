# Filter moabs.DMR to keep only regions whose seqnames match with gene.parts$TSSes

# Get unique seqnames from gene.parts$TSSes
tss_seqnames <- unique(seqnames(gene.parts$TSSes))

# Filter moabs.DMR to keep only those whose seqnames are in tss_seqnames
filtered_dmr <- moabs.DMR[seqnames(moabs.DMR) %in% tss_seqnames]

# Print summary of filtering
cat("Original number of DMRs:", length(moabs.DMR), "\n")
cat("Filtered number of DMRs:", length(filtered_dmr), "\n")

# Assign filtered results back to moabs.DMR if needed
# moabs.DMR <- filtered_dmr
