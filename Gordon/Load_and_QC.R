library(Seurat)
library(tibble)
library(dplyr)
library(cowplot)
library(ggplot2)
setwd("~/Documents/Projects/...") ## SET

##############
## Script to summarize, read in, and filter

## Assumes working directory contains a directory called GEX, which contains directories
## for each sample (with appropriate names), each of these subdirectories should include
## the metrics_summary.csv file and filtered_feature_bc_matrix.h5 file

##############
## Summarizer
## Simply collates each metrics_summary.csv and outputs to clipboard for quick summary of
## all runs which can be pasted into excel/word/email for sharing or reporting.

dirs <- list.dirs("./GEX", recursive = F)
summ.list <- list()
for (dir in dirs) {
  summ.list[[gsub("./GEX/", "", dir)]] <- read.csv(paste0(dir, "/metrics_summary.csv"))
}
summaries <- do.call(rbind, summ.list)
clipr::write_clip(summaries)

##############
## Read in data

## Create seurat objects for each sample, placing each in a list, followed by quick checks
## of feature and %mito reads

## Create seu
seu.list <- list()
for (dir in list.dirs("./GEX")[-1]) {
  name <- gsub("./GEX/", "", dir)
  seu.list[[name]] <- Read10X_h5(paste0(dir, "/filtered_feature_bc_matrix.h5"))
  seu.list[[name]] <- CreateSeuratObject(seu.list[[name]], project = name, min.cells = 3)
  seu.list[[name]] <- RenameCells(seu.list[[name]], add.cell.id = name)
}

clipr::write_clip(as.data.frame(unlist(lapply(seu.list, function(x) ncol(x))))) ## writes cell counts by sample pre filtering to clipboard

## Mito + feature check ## NOTE: change pattern to "^MT-" for human
seu.list <- lapply(seu.list, function(x) PercentageFeatureSet(x, pattern = "^mt-", col.name = "percent.mito"))
plot_grid(plotlist = lapply(seu.list, function(x) VlnPlot(x, features = "percent.mito", y.max = 100, pt.size = 0.1) + geom_hline(yintercept = 10, colour = "red") & NoLegend()), nrow = 1, align = "h", axis = "tblr")
plot_grid(plotlist = lapply(seu.list, function(x) VlnPlot(x, features = "nFeature_RNA", pt.size = 0.1, y.max = 8000) + geom_hline(yintercept = 400, colour = "red") & NoLegend()), nrow = 1, align = "hv", axis = "tblr") 

## Filters
orig.count <- lapply(seu.list, function(x) ncol(x)) ## For multiplet.rate calc in DoubletFinder step (see doublet removal script).
seu.list <- lapply(seu.list, function(x) subset(x, subset = percent.mito < 10))
seu.list <- lapply(seu.list, function(x) subset(x, subset = nFeature_RNA > 400))

clipr::write_clip(as.data.frame(unlist(lapply(seu.list, function(x) ncol(x))))) ## writes cell counts by sample post filtering to clipboard