#!/usr/bin/env Rscript

#' A list of two expression matrices for two batches. The first batch contains 400 cells of
#' three cell types, fibroblasts, macrophages and endothelial cells. And the second batches
#' has 500 cells of the same three cell types.
#'
#' @docType data
#' @usage data("data_SMNN")
#' @examples
#' # Load the example data data_SMNN
#' data("data_SMNN")
#' 
#' # Provide the marker genes for cluster matching
#' markers <- c("Col1a1", "Pdgfra", "Ptprc", "Pecam1")
#' 
#' # Specify the cluster labels for each marker gene
#' cluster.info <- c(1, 1, 2, 3)
#' 
#' # Call function unifiedClusterLabelling to identify the corresponding clusters between two batches
#' matched_clusters <- unifiedClusterLabelling(data_SMNN$batch1.mat, data_SMNN$batch2.mat, features.use = markers, cluster.labels = cluster.info, min.perc = 0.3)
#' 
#' # Set python version to be compatible with SMNNcorrect implementation
#' library(reticulate)
#' use_python("/nas/longleaf/apps/python/3.5.1/bin/python3")
#'
#' # Perform batch effect correction using SMNNcorrect
#' corrected.results <- SMNNcorrect(batches = list(data_SMNN$batch1.mat, data_SMNN$batch2.mat), batch.cluster.labels = matched_clusters, matched.clusters = c(1,2,3), k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE)
"data_SMNN"
