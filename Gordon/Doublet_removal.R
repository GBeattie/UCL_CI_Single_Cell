library(DoubletFinder)

####################
## Doublet removal (snippet), using DoubletFinder. Input is a list of seurat objects, follows
## authors recommended pipeline generating optimal pK value by sample. Normally run if any
## samples in analysis exceed 10,000 cells. Takes long to run so normally I save downstream metadata
## (post-integration, UMAP and clustering) which has doublets removed as a proxy for this step.

for (samp in names(seu.list)) {
  ## Subset and find sample specific pK
  seu.list[[samp]] <- seu.list[[samp]] %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  sweep <- paramSweep_v3(seu.list[[samp]], PCs = 1:10, sct = F)
  sweep.summ <- summarizeSweep(sweep, GT = FALSE)
  pK <- as.numeric(as.character((find.pK(sweep.summ) %>% slice_max(BCmetric, n = 1))$pK))
  ## Find sample specifid nExp
  homotypic.prop <- modelHomotypic(seu.list[[samp]]@meta.data$seurat_clusters)
  multiplet.rate <- 0.01*(orig.count[[samp]]/1000)
  nExp_poi <- round(multiplet.rate*nrow(seu.list[[samp]]@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Detect doublets and store
  seu.list[[samp]] <- doubletFinder_v3(seu.list[[samp]], PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, sct = F)
  doublet.meta <- FetchData(seu.list[[samp]] , vars = grep("DF.class", colnames(seu.list[[samp]]@meta.data), value = T))
  seu.list[[samp]] <- seu.list[[samp]][, which(doublet.meta == "Singlet")]
}
clipr::write_clip(as.data.frame(unlist(lapply(seu.list, function(x) ncol(x))))) ## writes cell counts by sample post doublet removal to clipboard