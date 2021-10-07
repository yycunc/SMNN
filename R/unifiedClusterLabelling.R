# !/usr/bin/env Rscript
# 04/26/2019
#
#' @title SMNN
#' 
#' @description This function unifiedClusterLabelling is used to match the clusters/cell types across multiple scRNA-seq batches. 
#' It takes as input raw expression matrices from two or more batches, a list of marker genes and their corresponding cluster labels.
#' It outputs cluster corresponding labels for the cells in each batch.
#' @usage unifiedClusterLabelling(..., features.use, cluster.labels, datatype="count", ident_list=NULL, cluster.use=NULL, cluster.names=NULL, min.exp.thresh=0, min.perc=0.6, min.average.expr=0)
#' @param batches is a list of two or more expression matrices where each row corresponds to a gene and each column corresponds to a single cell. All matriecs should contain the same number of rows (i.e., all batches should have the exact same gene set and in the same order).
#' @param features.use is a vector of marker genes used to match clusters across batches.
#' @param cluster.labels specifies the corresponding cluster label for each marker gene.
#' @param datatype defines the type of data, which can be "count", "CPM", "RPKM" and "FPKM". Default is "count".
#' @param ident_list is a list specifying cluster identities of the cells from each batch.
#' @param cluster.use defines the clusters used for matching.
#' @param cluster.names specifies the labels of clusters.
#' @param min.exp.thresh sets the minimum expression value for each marker genes in each cell. Default is 0.
#' @param min.perc sets the minimum percentage of cells whose expression level of marker gene \code{>=min.exp.thresh} during cluster definition. Default is 0.6.
#' @param min.average.expr sets the minimum average expression value of each marker gene during cluster definition. Default is 0.
#' @return unifiedClusterLabelling returns a list of unified cluster labels for the cells in each batch.
#' @author Yuchen Yang <yyuchen@email.unc.edu>, Gang Li <franklee@live.unc.edu>, Huijun Qian <hjqian@live.unc.edu>, Yun Li <yunli@med.unc.edu>
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
#' @import Matrix
#' @importFrom reshape melt
#' @export
unifiedClusterLabelling <- function(batches, features.use, cluster.labels, datatype="count", ident_list=NULL, cluster.use=NULL, cluster.names=NULL, min.exp.thresh=0, min.perc=0.6, min.average.expr=0) {
  group.used_all <- NULL
  nbatches <- length(batches)
  if (nbatches < 2L) { 
    stop("at least two batches must be specified") 
  }
  if (is.null(features.use) == 'TRUE'){
    stop("Error : a list of marker genes is required for cluster matching.")
  }
  if (is.null(cluster.labels) == 'TRUE'){
    stop("Error : a list of cluster labels is required for defining markers for certain clusters.")
  }
  if (!is.null(ident_list)){
    if (length(ident_list) != nbatches) { 
      stop("Error : ident_list must be of the same length as the number of batches.") 
    }
  } else {
    for (b in 1:nbatches){
      ident <- seurat_Smnn(batches[[b]], datatype = datatype, SEED = 1)
      ident_list <- c(ident_list, list(ident))
    }
  }
  
  for (b in 1:nbatches){
    data <- batches[[b]]
    ident <- ident_list[[b]]
    features.use <- features.use[features.use %in% rownames(data)]
    if (is.null(cluster.use)){
      cluster.use_b <- levels(ident)
    } else {
      cluster.use_b <- cluster.use
    }
    if (is.null(cluster.names)){
      cluster.names_b <- cluster.use_b
    } else{
      cluster.names_b <- cluster.names
    }
    if (length(cluster.names_b) != length(cluster.use_b)){
      print("Error : cluster.names must be of the same length as the clusters.use/ number of clusters. Using cluster numbers as labels ...")
      cluster.names_b <- cluster.use_b
    }
    
    #Initialize matrix of percent expressing cells
    PercMat <- matrix(0, nrow = length(features.use), ncol = 0)
    rownames(PercMat) <- features.use; 
    
    #Initialize matrix of average transcript levels
    ExpMat <- PercMat;
    
    #Count mat
    Count.mat <- data[features.use, colnames(data)]
    
    if (length(features.use) > 1){
      for (i in cluster.use_b){
        cells.in.cluster <- names(ident)[which(ident == i)]
        vec.exp <- apply(data[features.use, cells.in.cluster], 1, function(x) sum(x > min.exp.thresh)/length(x)) 
        PercMat <- cbind(PercMat,vec.exp)
        
        vec.exp <- apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x > min.exp.thresh) > 1){ mean(x[x > min.exp.thresh]) } else {sum(x)})
        ExpMat <- cbind(ExpMat, vec.exp)
      }
    } else if (length(features.use) == 1){
      for (i in cluster.use_b){
        cells.in.cluster <- names(ident)[which(ident == i)]
        vec.exp <- sum(data[features.use, cells.in.cluster] > min.exp.thresh)/length(data[features.use, cells.in.cluster]) 
        PercMat <- cbind(PercMat,vec.exp)
        
        if (sum(Count.mat[cells.in.cluster] > 0) > 1){
          vec.exp <- mean(Count.mat[cells.in.cluster][Count.mat[cells.in.cluster] > 0])
        } else {
          vec.exp <- sum(Count.mat[cells.in.cluster])
        }
        ExpMat <- cbind(ExpMat, vec.exp)
      }
    }
    colnames(ExpMat) <- cluster.names_b
    colnames(PercMat) <- cluster.names_b
    
    #if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
    #if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
    
    
    ExpVal <- melt(ExpMat)
    PercVal <- melt(PercMat)
    colnames(ExpVal) <- c("gene","cluster","nTrans")
    ExpVal$percExp <- PercVal$value*100
    
    ### define clusters according to corresponding markers of each cell type
    group.index <- matrix(0, nrow = length(cluster.names_b), ncol = 0)
    rownames(group.index) <- cluster.names_b
    for (i in 1:length(cluster.labels)){
      index <- as.numeric(ExpVal$percExp[which(ExpVal$gene == features.use[i])] >= min.perc*100 & ExpVal$nTrans[which(ExpVal$gene == features.use[i])] >= min.average.expr)
      group.index <- cbind(group.index, index)
    }
    colnames(group.index) <- cluster.labels
    
    group.assigned <- rep(0, length(cluster.names_b))
    for (i in 1:length(cluster.names_b)){
      if (sum(group.index[i,]) == 0){
        group.assigned[i] <- 0
      } else if (sum(group.index[i,]) >= 1){
        cluster.multi_assigned <- unique(cluster.labels[which(group.index[i,] != 0)])
        if (length(cluster.multi_assigned) == 1){
          group.assigned[i] <- cluster.multi_assigned
        }
      }
    }
    
    dist.calc <- function(x1, x2){
      sum <- 0
      for (i in x1){
        for (j in x2){
          sum <- sum + sqrt(sum((i - j) ^ 2))
        }
      }
      dist <- sum/(length(x1) * length(x2))
      
      return(dist)
    } 
    
    ### refine cell type labels for those clusters expressing markers of more than one cell types
    for (i in 1:length(cluster.names_b)){
      if (sum(group.index[i,]) > 1){
        cluster.multi_assigned <- cluster.labels[which(group.index[i,] != 0)]
        if (length(unique(cluster.multi_assigned)) > 1){
          dist.i <- NULL
          for (j in 1:length(cluster.multi_assigned)){
            if (sum(group.assigned == cluster.multi_assigned[j]) == 0){
              dist <- 0
            } else {
              cell.object <- which(ident == cluster.names_b[i])
              cell.query <- NULL
              for (c in cluster.names_b[which(group.assigned == cluster.multi_assigned[j])]){
                cell.query <- c(cell.query, which(ident == c))
              }
              dist <- dist.calc(data[features.use[j],cell.query], data[features.use[j], cell.object])
            }
            dist.i <- c(dist.i, dist)
          }
          group.assigned[i] <- cluster.multi_assigned[which(dist.i == min(dist.i))]
        }
      }
    }
    
    ### assign cell type label to each cell
    group.used <- rep(0, length(ident))
    for (c in 1:length(group.assigned)){
      group.used[which(ident == (c-1))] <- group.assigned[c]
    }
    
    group.used_all <- c(group.used_all, list(group.used))
  }
  
  return(group.used_all)
}


#' @import Seurat
#' @import dplyr
#' @import Matrix
seurat_Smnn <- function(inputTags, datatype = "count", resolution = 0.9, seurat_min_cell=200, resolution_min = 1.2, SEED = 1){
  seuratOUTPUT <- NULL
  
  # Initialize the Seurat object with the raw data (non-normalized data)
  # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
  seuratOUTPUT <- CreateSeuratObject(inputTags, min.cells = 0, min.features = 0, project = "single-cell clustering")
  
  # Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
  if (datatype == "count"){
    seuratOUTPUT = NormalizeData(object = seuratOUTPUT, normalization.method = "LogNormalize", scale.factor = 10000)
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM"){
    raw.data <- GetAssayData(object = seuratOUTPUT, slot = "raw.data")
    normalized.data <- log(raw.data+1)
    colnames(x = normalized.data) <- colnames(x = raw.data)
    rownames(x = normalized.data) <- rownames(x = raw.data)
    seuratOUTPUT <- SetAssayData(object = seuratOUTPUT, assay.type = "RNA",slot = "data", new.data = normalized.data)
  }
  
  # Detection of variable genes across the single cells
  seuratOUTPUT <- FindVariableFeatures(object = seuratOUTPUT, selection.method = "vst")
  
  # Regress out unwanted sources of variation
  seuratOUTPUT <- ScaleData(object = seuratOUTPUT, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
  
  ### Perform linear dimensional reduction
  seuratOUTPUT <- RunPCA(object = seuratOUTPUT, features = VariableFeatures(object = seuratOUTPUT), verbose = FALSE)
  
  seuratOUTPUT <- FindNeighbors(object = seuratOUTPUT, dims = 1:20)
  if (length(inputTags[1,]) >= seurat_min_cell){
    ### Clustering the cells by Seurat
    seuratOUTPUT <- FindClusters(object = seuratOUTPUT, algorithm = 3, resolution = resolution, verbose = FALSE, random.seed = SEED)
  } else {
    resolution <- resolution_min
    seuratOUTPUT <- FindClusters(object = seuratOUTPUT, algorithm = 3, resolution = resolution_min, verbose = FALSE, random.seed = SEED)
  }
  
  return(seuratOUTPUT@active.ident)
}
