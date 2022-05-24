library(Seurat)
library(tidyr)
library(sctransform)

####################
## Integration (snippet)
## This is a basic pipeline in which three integration (or batch effect) methods can be tried and checked to see if
## any of them are suitable for the data. Occasionally modified to experiment with any known batches, different batch
## correction algorithms, or removal of troublesome genes from expression matrix or reductions (not shown here). The
## goal is to achieve good overlap between samples with the least ammount of correction, hence why no batch correction
## is tried first, followed by Harmony (good for weak/medium effects) and CCA (good for strong effects).
## seu.list = list of seurat objects.

## No Batch Correction
seu.list <- lapply(seu.list, function(x) SCTransform(x))
features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)
seu.int <- merge(seu.list[[1]], y=seu.list[2:length(seu.list)])
VariableFeatures(seu.int) <- features
seu.int <- RunPCA(seu.int, verbose = FALSE)
seu.int <- RunUMAP(seu.int, dims = 1:30, verbose = FALSE)
seu.int <- FindNeighbors(seu.int, dims = 1:30, verbose = FALSE)
for (res in c(0.8, 0.5, 0.3)) {
  seu.int <- FindClusters(seu.int, verbose = FALSE, resolution = res)}
DimPlot(seu.int, group.by = c("orig.ident", grep("SCT_snn_", colnames(seu.int@meta.data), value = T)), label = T)

## Harmony
seu.int <- RunHarmony(seu.int, group.by.vars = "orig.ident", assay.use = "SCT")
seu.int <- RunUMAP(seu.int, reduction = "harmony", dims = 1:30)
seu.int <- FindNeighbors(seu.int, dims = 1:30, reduction = "harmony", verbose = FALSE)
for (res in c(0.8, 0.5, 0.3)) {
seu.int <- FindClusters(seu.int, verbose = FALSE, resolution = res)}
DimPlot(seu.int, group.by =  c("orig.ident", grep("SCT_snn_", colnames(seu.int@meta.data), value = T)), label = T)

## CCA
seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", anchor.features = features)
seu.int <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
seu.int <- RunPCA(seu.int, verbose = TRUE)
seu.int <- RunUMAP(seu.int, dims = 1:30, verbose = FALSE)
seu.int <- FindNeighbors(seu.int, dims = 1:30, verbose = FALSE)
for (res in c(0.8, 0.5, 0.3)) {
seu.int <- FindClusters(seu.int, verbose = FALSE, resolution = res)
}
DimPlot(seu.int, group.by =  c("orig.ident", grep("integrated_snn_", colnames(seu.int@meta.data), value = T)), label = T)
