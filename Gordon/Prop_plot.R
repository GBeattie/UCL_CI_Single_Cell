library(ggplot2)
library(dplyr)

################
## Prop plot (snippet)
## Bulky but quick count/proportion plots with split UMAPs.
## Will plot cluster by replicate (modify "seurat_clusters" and "replicate", alongside surat object name, to suit needs).

prop.plot <- plot_grid(
  ggplot(data.frame(table(seu.int@meta.data$replicate, seu.int$seurat_clusters)) %>% filter(Freq != 0), aes(fill = Var1, y = Freq, x = Var2)) +
    geom_col(position=position_dodge(width=0.75), aes(group = Var1)) +
    labs(x = "Cluster", y = "Cell Freq.", fill = "Sample") +
    theme_classic() + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    scale_fill_brewer(palette="Set1") +
    scale_colour_brewer(palette="Set2") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 4800)) 
  ,
  ggplot(data.frame(prop.table(table(seu.int@meta.data$replicate, seu.int$seurat_clusters), c(1))), aes(fill = Var1, y = Freq, x = Var2)) +
    geom_col(position=position_dodge(width=0.75), aes(group = Var1)) +
    labs(x = "Cluster", y = "Proportion", fill = "Sample") +
    theme_classic() + 
    theme(axis.text = element_text(size = 10, angle = 45, hjust = 1)) +
    scale_fill_brewer(palette="Set1") +
    scale_colour_brewer(palette="Set2") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.24)) 
  , ncol = 1)

plot_grid(DimPlot(seu.int, split.by = "orig.ident", pt.size = NULL, group.by = "seurat_clusters", ncol = 1, reduction = "umap"), prop.plot, ncol = 2, rel_widths = c(0.35, 1))