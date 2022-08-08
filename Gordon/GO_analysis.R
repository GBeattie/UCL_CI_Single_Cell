library(clusterProfiler)
library(org.Hs.eg.db)
# library(org.Mm.eg.db) ## for mouse
library(dplyr)

#########
## GO (snippet)
## quick GO analyses, bulk.de.list is a list of DEs generated through FindMarkers. Script will take top +ve and
## -ve log2FC genes and run through GO analysis, followed by simple dotplots.
## NOTE: GO filter can be changed (see note below)

go.list <- list()
for (i in c("+", "-")) {
  for (de in names(bulk.de.list)) {
    name <- paste0(de, " ", i)
    if(i == "+") { genes <- filter(bulk.de.list[[de]], avg_log2FC > 0.25)$Gene }
    else { genes <- filter(bulk.de.list[[de]], avg_log2FC < -0.25)$Gene }
    go.list[[name]] <- enrichGO(gene = genes,
                                     OrgDb    = org.Hs.eg.db, ## change for mouse
                                     ont      = "BP",
                                     keyType = "SYMBOL", pool = T,
                                     qvalueCutoff = 0.05
    )
  }
}

## Plotting. Note: filtering (i.e. pooling hierarchically) GO terms on the fly with gofilter command,
## level 1 = highest level, level 5+ (or remove filter altogether) = lowest level (i.e. most terms)
go.bulk.plot.list <- list()
for (go in names(go.list)) {
  go.bulk.plot.list[[go]] <- dotplot(gofilter(go.list[[go]], level = 3), showCategory=15, font.size = 8.5) + ggtitle(go) +
    scale_y_discrete(label = function(x) stringr::str_wrap(stringr::str_trunc(x, 100), width = 50))
}

plot_grid(plotlist = go.bulk.plot.list, nrow = 2)