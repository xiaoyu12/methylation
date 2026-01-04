library(ggplot2)
library(tidyverse)

prepare_kobas_data <- function(KOBAS_file, p_threshold = 0.05, direction_label = NULL) {
  enriched <- read.delim(KOBAS_file, comment.char = "#", blank.lines.skip = TRUE)
  top <- enriched %>% filter(P.Value <= p_threshold)

  colnames(top)[7] <- "qvalue"
  colnames(top)[4] <- "Gene_number"

  top <- top[order(top$P.Value, decreasing = TRUE), ]
  top$Term <- factor(top$Term, levels = top$Term)
  top$Rich.factor <- top$Gene_number / top$Background.number

  if (!is.null(direction_label)) {
    top$Direction <- direction_label
    top$Term_for_plot <- paste(direction_label, top$Term, sep = "__")
  }

  top
}

plot_kobas_bubble <- function(KOBAS_file, p_threshold = 0.05, out_img = "kobas_plot.png",
                              width = NA, height = NA, units = "in") {
  top <- prepare_kobas_data(KOBAS_file, p_threshold)

  sp <- ggplot(top, aes(x = Rich.factor, y = Term, size = Gene_number, color = P.Value)) +
    geom_point(alpha = 0.7) +
    theme_bw() +
    scale_color_gradientn(colours = rev(rainbow(20)), limits = c(0, p_threshold)) +
    ggtitle("Statistics of Pathway Enrichment") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Proportion of genes") +
    ylab("")

  ggsave(out_img, sp, dpi = 1200, width = width, height = height, units = units)

  return(top)
}

kobas <- plot_kobas_bubble("kobas_ehux_DE_shared_up.txt",
                           p_threshold = 0.08,
                           out_img = "bubble_kegg_DE_shared_up.png", height = 5)

kobas <- plot_kobas_bubble("kobas_DMGs.txt", out_img = "bubble_kegg_DMGs.png", height = 5)
kobas <- plot_kobas_bubble("kobas_DMPs.txt", out_img = "bubble_kegg_DMPs.png",
                           p_threshold = 0.08, height = 5)

kobas <- plot_kobas_bubble("kobas_DE_E1516_vs_E217_up.txt", out_img = "bubble_kegg_E217_up.png",
                           p_threshold = 0.05, height = 5)

kobas <- plot_kobas_bubble("kobas_DE_E1516_vs_E217_down.txt", out_img = "bubble_kegg_E217_down.png",
                           p_threshold = 0.05, height = 5)

# Combine up and down results into a single figure (up on top, down on bottom) with different shapes
up_data <- prepare_kobas_data("kobas_DE_E1516_vs_E217_up.txt", p_threshold = 0.05, direction_label = "Up")
down_data <- prepare_kobas_data("kobas_DE_E1516_vs_E217_down.txt", p_threshold = 0.05, direction_label = "Down")

# Preserve individual ordering within each half
up_levels <- up_data$Term_for_plot
down_levels <- down_data$Term_for_plot
combined <- bind_rows(up_data, down_data) %>%
  mutate(Term_for_plot = factor(Term_for_plot, levels = c(up_levels, down_levels)),
         Direction = factor(Direction, levels = c("Up", "Down")))

combined_plot <- ggplot(combined,
                        aes(x = Rich.factor, y = Term_for_plot, size = Gene_number,
                            color = P.Value, shape = Direction)) +
  geom_point(alpha = 0.7) +
  facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0.15, NA)) +
  scale_y_discrete(labels = function(x) sub("^[^_]+__", "", x)) +
  scale_color_gradientn(colours = rev(rainbow(20)), limits = c(0, 0.05)) +
  scale_shape_manual(values = c("Up" = 19, "Down" = 19)) +  # Up = square (15), Down = circle (19)
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Statistics of Pathway Enrichment (M217 Up vs Down)",
       x = "Proportion of genes", y = "")

ggsave("bubble_kegg_E217_combined.png", combined_plot, dpi = 1200, width = 7, height = 5, units = "in")
