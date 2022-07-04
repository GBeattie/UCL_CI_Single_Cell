#########################
## Custom umap

## Packages
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

## Script
umap.df <- FetchData(seu, vars = c("UMAP_1", "UMAP_2", "cluster.id")) ## change to match desired reduction and grouping (throughout script)
umap.cols <- DiscretePalette(length(unique(umap.df$cluster.id)), palette = "polychrome") 
ggrepel.df <- umap.df %>% group_by(cluster.id) %>% summarise_at(c("UMAP_1", "UMAP_2"), mean) 
ggrepel.df <- rbind(ggrepel.df, umap.df %>% mutate(cluster.id = "")) 
ggplot(umap.df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(fill = cluster.id), size = 1.5, colour = "black", pch = 21, stroke = 0.1) + 
  geom_text_repel(ggrepel.df,
                  mapping = aes(x = UMAP_1, y = UMAP_2, label = cluster.id), 
                  max.overlaps = Inf, box.padding = 0.1, force = 10, max.time = 1, force_pull = 2) +
  scale_fill_manual(values = umap.cols) +
  xlim(-9, 9) + ylim(-9, 9) + ## alter if too much/little room around points
  theme_void()