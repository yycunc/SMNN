## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE---------------
library(knitr)
opts_chunk$set(echo = TRUE)

## ----Set up the library, message=TRUE--------------------------------------
library("SMNN")

## ----set up for input expression data--------------------------------------
data("data_SMNN")
dim(data_SMNN$batch1.mat)
data_SMNN$batch1.mat[1:5, 1:5]
dim(data_SMNN$batch2.mat)

## ----define the marker genes for cluster matching, warning=FALSE------------
# Maker genes
markers <- c("Col1a1", "Pdgfra", "Ptprc", "Pecam1")
# Corresponding cell type labels for each marker gene
cluster.info <- c(1, 1, 2, 3)

## ----results='hide', fig.show="hide", message=FALSE-------------------------
matched_clusters <- unifiedClusterLabelling(data_SMNN$batch1.mat, data_SMNN$batch2.mat, features.use = markers, cluster.labels = cluster.info, min.perc = 0.3)

## ----set python version, results='hide'-------------------------------------
library(reticulate)
# please change the path below to your local path for python
use_python("/nas/longleaf/apps/python/3.5.1/bin/python3")

## ----perform batch effect correction using SMNNcorrect----------------------
corrected.results <- SMNNcorrect(data_SMNN$batch1.mat, data_SMNN$batch2.mat, batch1.cluster.labels = matched_clusters[[1]], batch2.cluster.labels = matched_clusters[[2]], num.defined.clusters = 3, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE)

## ----output from SMNNcorrect------------------------------------------------
# Output after correction for batch 
# Output (1): the batch-corrected expression matrix
corrected.results$corrected[[2]][1:10,1:10]
# Output (2): mutual nearest neighbor information
corrected.results$pairs[[2]]
