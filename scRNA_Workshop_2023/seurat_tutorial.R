## R script to perform analysis of mouse heart and neuron cells

## Load packages
library(Seurat)
library(SeuratObject)
library(harmony) ## Only if correcting for batch
library(ggplot2)
library(dplyr)
library(EnhancedVolcano) ## Follow install instructions here: https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#installation
library(pheatmap)

## Set working directory
setwd("~/Documents/Tutorial")

## Load in matrices
mtx.heart <- Read10X_h5(filename = "./GEX/heart/filtered_feature_bc_matrix.h5")
mtx.neuron <- Read10X_h5(filename = "./GEX/neuron/filtered_feature_bc_matrix.h5")

## Create Seurat Objects
seu.heart <- CreateSeuratObject(mtx.heart, project = "Heart")
seu.neuron <- CreateSeuratObject(mtx.neuron, project = "Neuron")

## Calculate mitochondiral %
seu.heart <- PercentageFeatureSet(seu.heart, pattern = "^mt-", col.name = "percent.mito")
seu.neuron <- PercentageFeatureSet(seu.neuron, pattern = "^mt-", col.name = "percent.mito")

## Plot feature counts and mitochondrial %
VlnPlot(seu.heart, features = c("nFeature_RNA", "percent.mito"))
VlnPlot(seu.neuron, features = c("nFeature_RNA", "percent.mito"))

## Filters
seu.heart <- subset(seu.heart, nFeature_RNA > 400 & percent.mito < 10)
seu.neuron <- subset(seu.neuron, nFeature_RNA > 400 & percent.mito < 10)

##########################################
## UMAP and clustering

seu.int <- merge(seu.heart, seu.neuron)
seu.int <- NormalizeData(seu.int)
seu.int <- FindVariableFeatures(seu.int)
seu.int <- ScaleData(seu.int)
seu.int <- RunPCA(seu.int)

## Keep if you don't want to batch correct
seu.int <- RunUMAP(seu.int, dims = 1:30)
seu.int <- FindNeighbors(seu.int)

## Keep if you want to batch correct
# seu.int <- RunHarmony(seu.int, group.by.vars = "orig.ident")
# seu.int <- RunUMAP(seu.int, dims = 1:30, reduction = "harmony")
# seu.int <- FindNeighbors(seu.int, reduction = "harmony")

seu.int <- FindClusters(seu.int)

## Visualise UMAP
DimPlot(seu.int)
DimPlot(seu.int, group.by = "orig.ident")

#########################################
## Save Seurat object, and load back in
## When closing RStudio (make sure to save script!) and reopening, start running script until
## you've loaded packaged and set working directory, then skip steps until readRDS command below

saveRDS(seu.int, "seurat_object.rds")
seu.int <- readRDS("seurat_object.rds")

#########################################
## Proportion plots
## Count
ggplot(seu.int@meta.data, aes(x = seurat_clusters, fill = orig.ident)) + geom_bar()
## proportion
ggplot(seu.int@meta.data, aes(x = seurat_clusters, fill = orig.ident)) + geom_bar(position = "fill")

#########################################
## Differential expression

## Cluster markers
cluster.markers <- FindAllMarkers(seu.int, logfc.threshold = 1, only.pos = T)
top3.markers <- group_by(cluster.markers, cluster) %>% slice_max(avg_log2FC, n = 3) %>% pull(gene)
DoHeatmap(seu.int, features = top3.markers)

## DE between samples
de.sample <- FindMarkers(seu.int, group.by = "orig.ident", ident.1 = "Heart", ident.2 = "Neuron", logfc.threshold = 1)
EnhancedVolcano(de.sample, x = "avg_log2FC", y = "p_val_adj", lab = rownames(de.sample), pCutoff = 0.05, FCcutoff = 3)

## DE between samples within a cluster
de.sample.clust3 <- FindMarkers(subset(seu.int, subset = seurat_clusters == "3"), group.by = "orig.ident", ident.1 = "Heart", ident.2 = "Neuron", logfc.threshold = 1)
EnhancedVolcano(de.sample.clust3, x = "avg_log2FC", y = "p_val_adj", lab = rownames(de.sample.clust3), pCutoff = 0.05, FCcutoff = 3)

#########################################
## Gene expression visualisations

## UMAPs
FeaturePlot(seu.int, features = "Xist")
FeaturePlot(seu.int, features = "Xist", split.by = "orig.ident")

FeaturePlot(seu.int, features = c("Xist", "Sox11"))
FeaturePlot(seu.int, features = c("Xist", "Sox11"), split.by = "orig.ident")

## Violin plots
VlnPlot(seu.int, features = "Xist")
VlnPlot(seu.int, features = "Xist", split.by = "orig.ident")

VlnPlot(seu.int, features = c("Xist", "Sox11"))
VlnPlot(seu.int, features = c("Xist", "Sox11"), split.by = "orig.ident")
VlnPlot(seu.int, features = c("Xist", "Sox11"), split.by = "orig.ident", stack = T)

## Heatmap
heatmap.genes <- c("Cd74", "Krt18", "Dcn", "H2-Ab1", "Mgp", "Col3a1", "Fabp7", "Krt8", "Col1a2", "Col1a1", "Krt19", "Ctsl")

DoHeatmap(seu.int, features = heatmap.genes)

seu.int.av <- AverageExpression(seu.int, return.seurat = T, features = heatmap.genes)
pheatmap(seu.int.av@assays$RNA@scale.data)

## Dotplot
DotPlot(seu.int, features = heatmap.genes)

#########################################
## Gene signature

gene.sig <- c("Myl7",
              "Crip1",
              "Myl3",
              "Eln",
              "Rac2",
              "Krt8")
seu.int <- AddModuleScore(seu.int, features = list(gene.sig), name = "heart_gene_signature")
FeaturePlot(seu.int, features = "heart_gene_signature1", split.by = "orig.ident")

#########################################
#########################################
#########################################
## Additional analyses

#########################################
## Interacting with the Seurat object
## To view the contents of a Seurat object you can click on the object in the viewer

## Metadata: a table of information for each cell, rows = cells, columns = data
View(seu.int@meta.data) ## open metadata table in viewer
colnames(seu.int@meta.data) ## show column names of metadata
seu.int$orig.ident ## print metadata column in console
table(seu.int$orig.ident) ## print frequency table
table(seu.int$orig.ident, seu.int$seurat_clusters) ## print 2-way frequency table

#########################################
## subclustering
seu.int.clust3 <- subset(seu.int, subset = seurat_clusters == "3") 
## Follow this with UMAP and clustering steps AFTER NormalizeData() function,
## all other steps can be carries out as is (just remember to use the new object,
## seu.int.clust3, instead of the original, seu.int)

#########################################
## GO term analysis
## There are programmatic ways of doing this, but for a quick and easy answer do this:
## Say you want GO terms for the cluster markers of cluster 3 (you may need to install clipr)
## the below command will extract those marker genes and put them on your clipboard,
## you can then paste these into http://geneontology.org/
clipr::write_clip(cluster.markers %>% filter(cluster == "3") %>% pull(gene))



