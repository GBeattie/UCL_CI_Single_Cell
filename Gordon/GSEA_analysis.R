library(clusterProfiler)
library(org.Hs.eg.db)
# library(org.Mm.eg.db) ## for mouse
library(dplyr)

#########
## GSEA (snippet)
## quick GSEA analyses, bulk.de.list is a list of DEs generated through FindMarkers. Script will rank genes
## by log2FC and run through GSEA analysis, followed by simple dotplots.

## bulk DE
gs.bulk.list <- list()
for (de in names(bulk.de.list)) {
  gs.bulk.list[[de]] <- gseGO(gene = rev(sort(setNames(bulk.de.list[[de]]$avg_log2FC, bulk.de.list[[de]]$Gene))),
                              OrgDb    = org.Hs.eg.db,
                              ont      = "ALL",
                              keyType = "SYMBOL",
                              minGSSize    = 100,
                              maxGSSize    = 500,
                              pvalueCutoff = 0.05
  )
}

## Plotting 
gs.bulk.plot.list <- list()
for (gs in names(gs.bulk.list)) {
  df <- as.data.frame(gs.bulk.list[[gs]]) %>% shonarrr::slice_min_max(NES, n = 20) %>% arrange(NES) %>% mutate(Description = factor(Description, levels = unique(Description)))
  gs.bulk.plot.list[[gs]] <- ggplot(df, aes(x = NES, y = Description, size = setSize, colour = -log10(p.adjust))) +
    geom_point() +
    scale_colour_gradient(low = "blue", high = "red") +
    geom_vline(xintercept = 0) +
    ggtitle(gs) +
    scale_y_discrete(label = function(x) stringr::str_wrap(stringr::str_trunc(x, 100), width = 50))
}

plot_grid(plotlist = gs.bulk.plot.list, nrow = 1, align = "hv", axis = "tblr")