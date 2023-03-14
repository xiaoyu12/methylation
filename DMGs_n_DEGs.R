# This file studies the overlap of DMGs and DEGs
DEGs = gene_regions %>% filter(ID %in% rownames(annot.DE))

DMG_n_DEGs = DEGs %>% filter(ID %in% DMGs$ID)
DMG_n_DEGs_up = DMG_n_DEGs %>% filter(ID %in% rownames(annot.DE.up))
DMG_n_DEGs_down = DMG_n_DEGs %>% filter(ID %in% rownames(annot.DE.down))

# write DMGs_n_DEG_up and down IDs to a file for enrichment analysis
write.table(DMG_n_DEGs_up$ID, file="DMG_n_DEGs_up.txt", quote=FALSE, row.names = FALSE)
cilium_up = c(106349, 117244, 122188, 194332, 198139, 198190, 207330, 211505, 212369, 
              216525, 219951, 223628, 229544, 230030, 243453, 252292, 438199, 63754, 66449, 77679)
cilium_up = as.character(cilium_up)
cilium_up = annot.DE [cilium_up, ]

