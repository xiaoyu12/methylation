library("stringr")

# Read a KOBAS output file
kobas_DE_up <- read.csv("kobas_DE_E1516_vs_E217_up.txt", comment.char = "#", sep = "\t")

# Read GO enrichment file
GO_DE_up <- read.csv("DE_E1516_vs_E217_up_ids.txt.GOseq.enriched", sep="\t", header = TRUE, row.names = 1)

enriched_genes_by_GO_term <- function(enrich_res, GO_term) {
  # Get the string of enriched genes, seperated by ", "
  enriched_genes <- enrich_res[GO_term, "gene_ids"]
  # split the string into a list
  enriched_genes <- strsplit(enriched_genes, ", ")
  return(enriched_genes[[1]])
}

# Get the genes enriched in GO:0030286 (dynein complex)
GO_dynein <- enriched_genes_by_GO_term(GO_DE_up, "GO:0030286")

# Get the genes in GO:0007017 (microtubule-based process)
GO_microtubule <- enriched_genes_by_GO_term(GO_DE_up, "GO:0007017")

# Loop through pathways in KOBAS output and find the overlapps with GO_dynein
for (i in 1:nrow(kobas_DE_up)) {
  ids = kobas_DE_up[i, "Input"]
  ids <- strsplit(ids, split="|", fixed=TRUE)
  ids <- ids[[1]]
  # find intersection
  common <- intersect(ids, GO_microtubule)
  if(length(common) > 0) {
    print(kobas_DE_up[i, "Term"])
    print(common)
  }
}
