library(ggplot2)
library(tidyverse)

plot_enrichment_bar_chart <- function(enrich_file, out_img = "enriched_barchart.png",
                                      p_threshold = 0.01) {
  # Read GO enrichment file
  enriched_go <- read.csv(enrich_file, sep="\t", header = TRUE, row.names = 1)
  # Filter by pvalues
  enriched_go <- enriched_go %>% filter(over_represented_pvalue < p_threshold) %>% as_tibble()
  # Sort by ontology and pvalues
  enriched_go <- enriched_go[order(enriched_go$ontology, enriched_go$over_represented_pvalue, 
                                   decreasing = TRUE),]
  # lock in factor level order
  enriched_go$term <- factor(enriched_go$term, levels = enriched_go$term)
  
  # Plot the barchart
  ggplot(enriched_go) + geom_bar(aes(x=term, y=numDEInCat, fill=ontology), stat="identity" ) + 
    geom_text(aes(x=term, y=numDEInCat, label=numDEInCat), hjust=-0.3, size=2.5) +
    coord_flip() +
    scale_fill_manual(breaks = c("BP", "CC", "MF"), values = c("blue3", "green3", "red3")) +
    theme_light() +
    scale_y_log10(limit=c(1, 2500), breaks=c(10, 100, 1000)) +
    theme(legend.direction = "vertical", legend.text = element_text(size=8),
          legend.key.width = unit(0.1, 'cm')) +
    labs(y="Num of Genes") +
    theme(axis.title.y = element_blank(),
          text=element_text(size=8))
  
  ggsave(out_img, dpi=1200)
}


plot_enrichment_bar_chart("DMGs.txt.GOseq.enriched", out_img = "DMGs.GOseq.enriched.png")

# Plot enriched GO terms DE between strains E217 and E1516
plot_enrichment_bar_chart("DE_E1516_vs_E217_ids.txt.GOseq.enriched", 
                          p_threshold = 0.001,
                          out_img = "DE_E1516_vs_E217.enriched.png")

# Plot enriched GO terms up-regulated in strains E217 vs. E1516
plot_enrichment_bar_chart("DE_E1516_vs_E217_up_ids.txt.GOseq.enriched", 
                          p_threshold = 0.001,
                          out_img = "DE_E1516_vs_E217.up.enriched.png")

# Plot enriched GO terms down-regulated in strains E217 vs. E1516
plot_enrichment_bar_chart("DE_E1516_vs_E217_down_ids.txt.GOseq.enriched", 
                          p_threshold = 0.001,
                          out_img = "DE_E1516_vs_E217.down.enriched.png")

# Plot enriched GO terms for upregulated genes shared in E217 vs. E1516 and 9S vs. 9
plot_enrichment_bar_chart("DE_shared_up.txt.GOseq.enriched", 
                          p_threshold = 0.001,
                          out_img = "DE_shared_up.enriched.png")

# Plot DMGs enrichment analysis file
plot_enrichment_bar_chart("DMGs.txt.GOseq.enriched",
                          p_threshold = 0.01,
                          out_img = "DMGs.GOseq.enriched.png")

# Plot DMPs enrichment analysis file
plot_enrichment_bar_chart("DMPs.txt.GOseq.enriched",
                          p_threshold = 0.01,
                          out_img = "DMPs.GOseq.enriched.png")


# Select top enriched GO terms for revigo for GO reduction and visualization
select_top_for_revigo <- function(enrich_file, p_threshold = 0.001, 
                                  cat_size = 200, out_file = "revigo_input.txt") {
  # Read GO enrichment file
  enriched_go <- read.csv(enrich_file, header = TRUE, sep="\t")
  # Filter by pvalues
  enriched_go <- enriched_go %>% filter(over_represented_pvalue < p_threshold)
  # Filter by category size. Only keep small and more specific categories
  enriched_go <- enriched_go %>% filter(numInCat < cat_size)
  
  write.table(enriched_go[, 1:2], file=out_file, 
              sep="\t", row.names = FALSE, quote=FALSE, col.names =FALSE)  
}

select_top_for_revigo("DE_shared_up.txt.GOseq.enriched")
select_top_for_revigo("DE_E1516_vs_E217_up_ids.txt.GOseq.enriched", 
                      out_file = "DE_E1516_vs_E217_up_revigo.txt")
select_top_for_revigo("DE_E1516_vs_E217_down_ids.txt.GOseq.enriched",
                      out_file = "DE_E1516_vs_E217_down_revigo.txt")
select_top_for_revigo("DMGs.txt.GOseq.enriched",
                      p_threshold = 0.01,
                      out_file = "DMGs_revigo.txt")
select_top_for_revigo("DMPs.txt.GOseq.enriched",
                      p_threshold = 0.01,
                      out_file = "DMPs_revigo.txt")
