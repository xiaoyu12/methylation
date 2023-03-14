library(ggplot2)
library(tydyverse)

plot_kobas_bubble <- function(KOBAS_file, p_threshold=0.05, out_img="kobas_plot.png",
                              width=NA, height=NA, units="in") {
  # Read KOBAS output file
  enriched <- read.delim(KOBAS_file,comment.char = "#", blank.lines.skip = TRUE)
  # Filter out pathway with pvalue > threshold
  top <- enriched %>% filter(P.Value <= p_threshold)
  
  colnames(top)[7] = "qvalue"
  colnames(top)[4] = "Gene_number"
  # Sort by pvalues
  top <- top[order(top$P.Value, decreasing = TRUE),]
  # lock in factor level order
  top$Term <- factor(top$Term, levels = top$Term)
  
  top$Rich.factor = top$Gene_number / top$Background.number
  sp <- ggplot(top, aes(x=Rich.factor, y=Term, size=Gene_number, color=P.Value)) + geom_point(alpha=0.7)
  sp <- sp + theme_bw() #+ scale_color_gradient(low = "blue", high = "red", limits=c(0,1)) 
  sp <- sp +scale_color_gradientn(colours=rev(rainbow(20)), limits=c(0, p_threshold)) 
  sp <- sp + ggtitle("Statistics of Pathway Enrichment")+ theme(plot.title = element_text(hjust = 0.5))
  sp + xlab("Rich factor") + ylab("")
  
  ggsave(out_img, dpi=1200, width=width, height = height, units=units)
  
  return (top)
}

kobas <- plot_kobas_bubble("kobas_ehux_DE_shared_up.txt", 
                           p_threshold = 0.08, 
                           out_img = "bubble_kegg_DE_shared_up.png", height=5)

kobas <- plot_kobas_bubble("kobas_DMGs.txt", out_img = "bubble_kegg_DMGs.png", height=5)
kobas <- plot_kobas_bubble("kobas_DMPs.txt", out_img = "bubble_kegg_DMPs.png", 
                           p_threshold = 0.08, height=5)

kobas <- plot_kobas_bubble("kobas_DE_E1516_vs_E217_up.txt", out_img = "bubble_kegg_E217_up.png",
                           p_threshold = 0.05, height=5)

kobas <- plot_kobas_bubble("kobas_DE_E1516_vs_E217_down.txt", out_img = "bubble_kegg_E217_down.png",
                           p_threshold = 0.05, height=5)
