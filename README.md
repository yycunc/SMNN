# SMNN
### SMNN: Batch Effect Correction for Single-cell RNA-seq data via Supervised Mutual Nearest Neighbor Detection

Batch effect correction has been recognized to be indispensable when integrating scRNA-seq data from multiple batches. A recent study proposed an effective batch effect correction method based on mutual nearest neighbors (MNN) across batches (Haghverdi *et al.*, 2018). However, the original MNN method is unsupervised in that it ignores cluster label information of single cells, which can further improve effectiveness of batch effect correction, particularly under realistic scenarios where true biological differences are not orthogonal to batch effect. 

Thus, we propose SMNN for batch effect correction of scRNA-seq data via supervised mutual nearest neighbor detection. SMNN either takes cluster/cell-type label information as input or infers cell types using scRNA-seq clustering. It then detects mutual nearest neighbors within matched cell types and corrects batch effect accordingly. Compared to MNN, SMNN provides improved merging within the corresponding cell types across batches, and retains more cell type specific features after correction, especially under realistic scenarios where different batches differ in many aspects including samples used, single cell capture technology employed, or library preparation approach adopted.

SMNN is maintained by Yuchen Yang [yyuchen@email.unc.edu] and Gang Li [franklee@live.unc.edu].

## News and Updates
May 23, 2019
* Version 0.99.0 released
  + First offical release
  + Now it can only work on Mac and Linux platform


## Brief introduction

The current implementation of SMNN encompasses two major steps: one optional clustering step and the other batch effect correction step. In the first step, SMNN takes the expression matrix as input, and performs clustering using Seurat v. 3.0 (Butler *et al.*, 2018). Corresponding clusters/cell types are then matched across batches based on marker genes specified by the user. This entire clustering step can be by-passed by feeding SMNN cell cluster labels. With cell cluster label information, SMNN searches mutual nearest neighbors within each cell type, and performs batch effect correction using the *SMNNcorrect* function.

In this tutorial, we will perform batch effect correction using SMNN in a toy example containing two batches. The first batch contains 400 cells from three cell types, namely fibroblasts, macrophages and endothelial cells. And the second batches has 500 cells from the same three cell types. Both two batches contain 3000 genes.


## Installation

SMNN package can be directly installed from GitHub with:
```{r installation}
install.packages("devtools")

devtools::install_github("yycunc/SMNN")
```


## Set up the library
```{r init, message=TRUE}
library("SMNN")
```


## Load the input expression matrix

Once installed, one can use the following command lines to load the data used in our toy example using the following command lines: 
```{r set up for input expression data}
data("data_SMNN")
dim(data_SMNN$batch1.mat)
data_SMNN$batch1.mat[1:5, 1:5]
dim(data_SMNN$batch2.mat)
```


## The optional clustering step
### Provide cell-type specific marker gene information

Two pieces of information are needed:
- Marker genes
- Corresponding cell type labels for each marker gene

```{r define the marker genes for cluster matching, warning=FALSE}
# Maker genes
markers <- c("Col1a1", "Pdgfra", "Ptprc", "Pecam1")
# Corresponding cell type labels for each marker gene
cluster.info <- c(1, 1, 2, 3)
```

### Harmonize cluster labels across batches

The function *unifiedClusterLabelling* is used to match/harmonize the clusters/cell type labels across multiple scRNA-seq batches. It takes as input raw expression matrices from two or more batches, a list of marker genes and their corresponding cluster labels, and outputs harmonized cluster label for every single cells across all batches.

```{r, results='hide', fig.show="hide", message=FALSE}
matched_clusters <- unifiedClusterLabelling(data_SMNN$batch1.mat, data_SMNN$batch2.mat, features.use = markers, cluster.labels = cluster.info, min.perc = 0.3)
```

### The batch effect correction step using SMNNcorrect function
#### Set python version to be compatible with SMNNcorrect implementation
In the *SMNNcorrection* function, we call several python functions from *mnnpy* package at [chriscainx/mnnpy](https://github.com/chriscainx/mnnpy) to speed up the computation.

Package *mnnpy* can be installed with pip install mnnpy,

or

git clone https://github.com/chriscainx/mnnpy.git <br>
cd mnnpy <br>
pip install .

To ensure successful execution of these python functions, we need to first set the python version (we recommand python3 for SMNN version 0.99).

```{r set python version, results='hide'}
library(reticulate)

# please change the path below to your local path for python
use_python("/nas/longleaf/apps/python/3.5.1/bin/python3")
```

## Perform batch effect correction
With harmonized cluster label information for single cells across batches, we perform batch effect correction. Specifically, we apply cosine normalization on both input and output data and set the number of mutual nearest neighbors at 20.

```{r perform batch effect correction using SMNNcorrect}
corrected.results <- SMNNcorrect(data_SMNN$batch1.mat, data_SMNN$batch2.mat, batch.cluster.labels = matched_clusters, num.defined.clusters = 3, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE)
```

*SMNNcorrect* function will output the following information: (1) the batch-corrected expression matrix for each batch; and (2) information regarding mutual nearest neighbors.

In the example below, we treat batch 1 as the reference batch and batch 2 as the batch to be corrected (such that batch 2 will be corrected towards the reference batch 1). Note that the reference batch (i.e., batch 1 in our example) will only applied cosine normalization.

```{r output from SMNNcorrect}
# Output after correction for batch
## Output (1): the batch-corrected expression matrix
corrected.results$corrected[[2]][1:10,1:10]
## Output (2): mutual nearest neighbor information
corrected.results$pairs[[2]]
```
