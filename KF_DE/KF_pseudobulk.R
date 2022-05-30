library(zellkonverter) #to read scanpy .h5ad to sce
library(SingleCellExperiment)
library(tidyverse)
library(stringr)
library(scuttle)
library(edgeR)
library(limma)


###############
## Pseuodbulk DE
## This script was inspired by https://doi.org/10.1038/s41467-021-25960-2
## limma() is used below, but any other bulk RNAseq approach (DESeq2, edgeR) would work


## Import
## Data from https://doi.org/10.1038/s41590-021-00931-3
## Only one cell type (activated CD8+Tem), CD38-ve sort only
sce = readH5AD('Documents/ResearchProjects/_scRNAseqWorkflows/Increate_v4/UCL_CI_SC/CD8TemPD1_only.h5ad')
## scanpy.X = sce@assay
sce@assays@data@listData[['counts']] <- sce@assays@data@listData[['X']] 
## covar: disease = healthy donor (HD) or mutiple myeloma (MM)
table(sce$donor, sce$disease)
covar='disease'
## clustering used was leiden_1.5

## Pseudobulk expression
## Pseudobulk to sample-cluster ("seq_name:leiden_1.5")
ids = paste0(sce$seq_name, ':', as.character(sce$leiden_1.5))
pb = scuttle::aggregateAcrossCells(sce,ids=ids, statistics='sum', use.assay.type = "counts")
pb$total_counts<-colSums(counts(pb))
## Remove null columns
colData(pb)<-colData(pb)[ , colSums(is.na(colData(pb))) == 0]

## DE with limma: data prep
pb.test = pb
## Min cells filter
minCells=10
pb.test<-pb.test[,pb.test$ncells > minCells]
## Removing genes that are lowly expressed
keep.exprs = edgeR::filterByExpr(pb.test, group=pb.test$seq_name)
pb.test = pb.test[keep.exprs,]
## Normalising gene expression distributions
pb.test.DGElist <- edgeR::calcNormFactors(pb.test, method = "TMM")

## DE with limma: inspecting data
## Inspect normalised expression by PCA
pca = prcomp(t(pb.test.DGElist$counts))$x
as_tibble(cbind(colData(pb),pca)) -> pca.pl
## Inspect covars of interest
## condition=disease
ggplot(pca.pl,aes(PC1,PC2,color=disease,label=seq_name))+geom_point(size=4)+geom_text(color='black',size=3)
## n cells (can see contributes to outliers)
ggplot(pca.pl,aes(PC1,PC2,color=ncells,label=seq_name))+geom_point(size=4)+geom_text(color='black',size=3)
## total_counts (similarly contributes)
ggplot(pca.pl,aes(PC1,PC2,color=total_counts,label=seq_name))+geom_point(size=4)+geom_text(color='black',size=3)
## can inspect i.e. stress genes
stress_exprssn = c(colMeans(pb.test.DGElist$counts[c('JUNB','JUN','FOS','FOSB'),]))
ggplot(pca.pl,aes(PC1,PC2,color=stress_exprssn,label=seq_name))+geom_point(size=4)+geom_text(color='black',size=3)

## DE with limma: performing tests
## Creating a design matrix and contrast matrix
## https://doi.org/10.12688%2Ff1000research.27893.1
condition<-colData(pb.test)[,covar]
design <- model.matrix(~0+condition)
colnames(design) <- gsub("group", "", str_replace_all(colnames(design),'condition',''))
makeContrasts(contrast = A-B, levels=c('A','B')) -> contr.matrix
colnames(contr.matrix)<-covar
rownames(contr.matrix)<-colnames(design)
## Removing heteroscedascity from count data
v <- voom(pb.test.DGElist$counts, design, plot=T) #also converts to EList
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

## DE with limma: inspecting results
summary(decideTests(efit)) 
## treat method calculates p-values from empirical Bayes moderated t-statistics with a 
## minimum log-FC requirement
tfit <- treat(vfit, fc=1.2)
dt <- decideTests(tfit)
print(summary(dt))
## Examining individual DE genes
res <- topTreat(tfit, coef=1, n=Inf) #coef selects which contrast's genes
res$gene<-rownames(res)
res = as_tibble(res)
## quick volcano
ggplot(res, aes(logFC, -log10(adj.P.Val)))+
  geom_point(size=1,alpha=0.5, aes(color=adj.P.Val<0.01))+
  geom_text(aes(label=ifelse(adj.P.Val<0.01,gene,'')),size=2)



