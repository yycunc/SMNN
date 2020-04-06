# !/usr/bin/env Rscript
# 04/26/2019
#
#' @title SMNN
#' 
#' @description This function SMNNcorrect is designed to perform supervised batch effect correction for scRNA-seq data by first identifying nearest neighbors (NNs) within corresponding clusters (or cell types) and then leveraging information from these NNs.
#' It takes as input raw expression matrices from two or more batches and a list of the unified cluster labels (output from unifiedClusterLabelling).
#' It outputs batch-corrected expression matrix for each batch.
#' @usage SMNNcorrect(batches, batch.cluster.labels, matched.labels=c(1,2,3), correct.others=FALSE, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE, subset.genes=NULL, order=NULL, n.jobs=NULL)
#' @param batches is a list of two or more expression matrices each corresponding to one batch, where each row corresponds to a gene, and each colname correspond to a cell. 
#' The number and order of rows should be identical across all maxtices (i.e., all batches should have the exact same gene set and in the same order).
#' @param batch.cluster.labels is a list of vectors specifying the cluster labels of each cell from each batch. Cells not belonging to any clusters should be set to 0.
#' SMNN performs batch effect correction without any prior knowledge on cell clusters if {batch.cluster.labels = NULL}.
#' @param matched.clusters specifies the cell clusters matched between two or more batches.
#' @param correct.others is a Boolean variable that defines whether to search nearest neighbors among the cells not belonging to any clusters. Default is FALSE, that is, cells not belonging to any clusters will not be considered as candidate nearest neighbors. 
#' @param k defines the maximum number of nearest neighbors to be identified. Default is 20.
#' @param sigma defines the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell. Default is 1.
#' @param cos.norm.in is a boolean variable that defines whether to do cosine normalization on input data before computing distances between cells.
#' Default is "TRUE".
#' @param cos.norm.out is a boolean variable that defines whether to do cosine normalization on output data before computing corrected expression results.
#' Default is "TRUE".
#' @param var.adj is a Boolean variable that indicates whether to do variance adjustment on the correction vectors. Default is "TRUE".
#' @param subset.genes is a vector specifying the gene set that is used for computing correction vectors.
#' Default is {subset.genes = NULL}, which means to use all the genes to compute conrrection vectors. 
#' @param order is a vector defining the reference batch and the order of the other batches to be corrected.
#' @param n.jobs specifies the number of parallel jobs. It would be set to the number of cores when \code{n.jobs = NULL}.
#' @return SMNNcorrect returns the following:
#'     \item corrected expression matrix for each batch
#'     \item information regarding NNs between the current batch under correction and the reference batch
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
#' 
#' # Set python version to be compatible with SMNNcorrect implementation
#' library(reticulate)
#' use_python("/nas/longleaf/apps/python/3.5.1/bin/python3")
#'
#' # Perform batch effect correction using SMNNcorrect
#' corrected.results <- SMNNcorrect(batches = list(data_SMNN$batch1.mat, data_SMNN$batch2.mat), batch.cluster.labels = matched_clusters, matched.clusters = c(1,2,3), k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE)
#' @import reticulate
#' @importFrom S4Vectors DataFrame Rle
#' @export
SMNNcorrect <- function(batches, batch.cluster.labels, matched.clusters, correct.others=FALSE, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE, subset.genes=NULL, order=NULL, n.jobs=NULL){
    #batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }

    if(!is.null(batch.cluster.labels)){
        ncluster.labels <- length(batch.cluster.labels)
        if(ncluster.labels != nbatches) { 
           stop("The number of cluster labels specified must be equal to the number of batches") 
        }
    }
    
    batches.t <- list()
    for (i in 1:nbatches){
        batches.t[[i]] <- t(batches[[i]])
    }

    if (is.null(n.jobs)){
       multiprocessing = import("multiprocessing")
       n.jobs <- multiprocessing$cpu_count()
    }
    
    subset.index <- match(subset.genes, rownames(batches[[1]]), nomatch = 0L)
    subset.index <- subset.index[subset.index != 0]
    if (length(subset.index) == 0){
       stop("At least one of the genes selected for correction vectors should exist in the input matrix")
    } else{
       subset.index <- as.character(subset.index - 1)
    }
    
    mnnpy <- import("mnnpy")
    print("Data preparation ...")
    prep.out <- mnnpy$utils$transform_input_data(datas=batches.t, cos_norm_in=cos.norm.in, cos_norm_out=cos.norm.out, var_index=as.character(c(0:ncol(batches.t[[1]]))), var_subset=subset.index, n_jobs=n.jobs)
    in.batches <- prep.out[[1]]
    out.batches <- prep.out[[2]]
    subset.genes <- prep.out[[3]]
    same.set <- prep.out[[4]]

    # Setting up the order.
    if (is.null(order)) {
        order <- seq_len(nbatches)
    } else {
        order <- as.integer(order)
        if (!identical(seq_len(nbatches), sort(order))) { 
            stop(sprintf("'order' should contain values in 1:%i", nbatches))
        }
    }
   
    # Setting up the variables.
    ref <- order[1]
    ref.batch.in <- in.batches[[ref]]
    if (!same.set) { 
        ref.batch.out <- out.batches[[ref]]
    }
    
    ### output for supervised batch correction
    output <- vector("list", nbatches)
    output[[ref]] <- t(out.batches[[ref]])
    mnn.list <- vector("list", nbatches)
    mnn.list[[ref]] <- DataFrame(current.cell=integer(0), other.cell=integer(0), other.batch=integer(0))
    original.batch <- rep(ref, nrow(ref.batch.in)) 

    # Looping through all the other batches.
    for (b in 2:nbatches) { 
        target <- order[b]
        other.batch.in <- in.batches[[target]]
        if (!same.set) { 
            other.batch.out <- out.batches[[target]]
        }
        
        if(!is.null(batch.cluster.labels)){
          # Supervised finding pairs of mutual nearest neighbours.
          s1_supervised <- integer()
          s2_supervised <- integer()
        
          for(c in matched.clusters){
            sets1 <- NULL          
            ref.batch.rank <- which(batch.cluster.labels[[1]] == c)
            other.batch.rank <- which(batch.cluster.labels[[target]] == c)

            ### MNNs for each subset. The index is python format (starting from 0 ..)
            print(paste0("Finding the mutual nearest neighbors for cell type ", c, " ..."))
            sets1 <- mnnpy$utils$find_mutual_nn(data1=ref.batch.in[which(batch.cluster.labels[[1]] == c),], data2=other.batch.in[which(batch.cluster.labels[[target]] == c),], k1=k, k2=k, n_jobs=n.jobs)
          
            s1.trace_back <- integer()
            for(i in 1:length(unlist(sets1[[1]]))){
              s1.trace_back <- c(s1.trace_back, ref.batch.rank[unlist(sets1[[1]])[i]+1])
            }
            s1_supervised <- c(s1_supervised, s1.trace_back)
          
            s2.trace_back <- integer()
            for(j in 1:length(unlist(sets1[[2]]))){
              s2.trace_back <- c(s2.trace_back, other.batch.rank[unlist(sets1[[2]])[j]+1])
            }
            s2_supervised <- c(s2_supervised, s2.trace_back)
          }
        
          if(correct.others == TRUE){
            sets1 <- NULL
          
            ref.batch.rank <- which(batch.cluster.labels[[1]] == 0)
            other.batch.rank <- which(batch.cluster.labels[[target]] == 0)
            sets1 <- mnnpy$utils$find_mutual_nn(data1=ref.batch.in[which(batch.cluster.labels[[1]] == 0),], data2=other.batch.in[which(batch.cluster.labels[[target]] == 0),], k1=k, k2=k, n_jobs=n.jobs)
  
            s1.trace_back <- integer()
            for(i in 1:length(unlist(sets1[[1]]))){
              s1.trace_back <- c(s1.trace_back, ref.batch.rank[unlist(sets1[[1]])[i] + 1])
            }
            s1_supervised <- c(s1_supervised, s1.trace_back)
          
            s2.trace_back <- integer()
            for(j in 1:length(sets1[[2]])){
              s2.trace_back <- c(s2.trace_back, other.batch.rank[unlist(sets1[[2]])[j] + 1])
            }
            s2_supervised <- c(s2_supervised, s2.trace_back)
          }
        
          # Change index from R format (starting from 1 ..) to python format (starting from 0 ..)
          s1 <- as.integer(s1_supervised)
          s2 <- as.integer(s2_supervised)
          sets_supervised <- list(first = s1, second = s2)
        } 
        else if(is.null(batch.cluster.labels)){
          batch.cluster.labels_ref = rep(0, dim(ref.batch.in)[1])
          batch.cluster.labels_other = rep(0, dim(other.batch.in)[1])
           
          # Unsupervised finding pairs of mutual nearest neighbours.
          s1_unsupervised <- integer()
          s2_unsupervised <- integer()
        
          sets1 <- NULL
          ref.batch.rank <- which(batch.cluster.labels_ref == 0)
          other.batch.rank <- which(batch.cluster.labels_other == 0)
          ### MNNs for each subset. The index is python format (starting from 0 ..)
          print("Unsupervised searching the mutual nearest neighbors...")
          sets1 <- mnnpy$utils$find_mutual_nn(data1=ref.batch.in[which(batch.cluster.labels_ref == 0),], data2=other.batch.in[which(batch.cluster.labels_other == 0),], k1=k, k2=k, n_jobs=n.jobs)
          
          for(i in 1:length(unlist(sets1[[1]]))){
            s1_unsupervised <- c(s1_unsupervised, ref.batch.rank[unlist(sets1[[1]])[i]+1])
          }
          
          for(j in 1:length(unlist(sets1[[2]]))){
            s2_unsupervised <- c(s2_unsupervised, other.batch.rank[unlist(sets1[[2]])[j]+1])
          }
        
          # Change index from R format (starting from 1 ..) to python format (starting from 0 ..)
          s1 <- as.integer(s1_unsupervised)
          s2 <- as.integer(s2_unsupervised)
          sets_unsupervised <- list(first = s1, second = s2)
        }

        # Computing the correction vector
        print("Computing the correction vector ...")
        correction.in <- compute_correction(ref.batch.in, other.batch.in, s1, s2, other.batch.in, sigma)

        if (!same.set) {
            correction.out <- compute_correction(ref.batch.out, other.batch.out, s1, s2, other.batch.in, sigma)
        }
       
        # Adjusting the shift variance; done after any SVD so that variance along the correction vector is purely technical.
        if (var.adj) { 
            print("Adjusting the shift variance ...")
            multiprocessing <- import("multiprocessing")
            builtins <- import_builtins()
            np <- import("numpy")
            p_n <- multiprocessing$Pool(n.jobs)
            #x= mnnpy$utils$adjust_s_variance(ref.batch.in, other.batch.in, other.batch.in[1,], vect[1,], 0.1)
            vect <- correction.in
            scaling <- c()
            for (i in 1:nrow(vect)){
                scaling <- c(scaling, adjust_s_variance(ref.batch.in, other.batch.in, other.batch.in[i,], vect[i,], sigma))
            }
            scaling.max <- np$fmax(scaling, 1)
            scaling <- np$ndarray$astype(scaling.max, np$float32)
            for (i in 1:nrow(correction.in)){
                correction.in[i,] <- correction.in[i,] * scaling[i]
            }

            #correction.in <- mnnpy$utils$adjust_shift_variance(ref.batch.in, other.batch.in, correction.in, sigma, n_jobs=1L)
            if (!same.set) {
            #x= mnnpy$utils$adjust_s_variance(ref.batch.in, other.batch.in, other.batch.in[1,], vect[1,], 0.1)
                vect <- correction.out
                scaling <- c()
                for (i in 1:nrow(vect)){
                    scaling <- c(scaling, adjust_s_variance(ref.batch.out, other.batch.out, other.batch.out[i,], vect[i,], sigma))
                }
#            scaling = p_n$starmap(mnnpy$utils$adjust_v_worker(ref.batch.in, other.batch.in, 0.1), builtins$zip(other.batch.in, vect), chunksize=floor(nrow(other.batch.in)/n_jobs)+1)
                scaling.max <- np$fmax(scaling, 1)
                scaling <- np$ndarray$astype(scaling.max, np$float32)
                for (i in 1:nrow(correction.in)){
                    correction.out[i,] <- correction.out[i,] * scaling[i]
                }
                #correction.out <- mnnpy$utils$adjust_shift_variance(ref.batch.out, other.batch.out, correction.out, sigma=sigma, n_jobs=1L, subset.row) 
            }
        }

        # Applying the correction and expanding the reference batch. 
        other.batch.in <- other.batch.in + correction.in
        ref.batch.in <- rbind(ref.batch.in, other.batch.in)
        if (same.set) {
            output[[target]] <- t(other.batch.in)
        } else {
            other.batch.out <- other.batch.out + correction.out
            ref.batch.out <- rbind(ref.batch.out, other.batch.out)
            output[[target]] <- t(other.batch.out)
        }

        # Storing the identities of the MNN pairs (using RLEs for compression of runs).
        mnn.list[[target]] <- DataFrame(current.cell=s2, other.cell=Rle(s1), other.batch=Rle(original.batch[s1]))
        original.batch <- c(original.batch, rep(target, nrow(other.batch.in)))
    }    
    
    # Formatting output to be consistent with input.
    names(output) <- names(batches)
    names(mnn.list) <- names(batches)
    final <- list(corrected=output, pairs=mnn.list)
    
    print("Done!")
    return(final)
}

compute_correction <- function(data1, data2, mnn1, mnn2, data2_or_raw2, sigma){
    np <- import("numpy")
    vect <- data1[mnn1,] - data2[mnn2,]
    mnn2_py <- as.integer(mnn2-1)
    mnn_output <- np$unique(mnn2_py, return_counts = TRUE)
    mnn_index <- mnn_output[[1]]
    mnn_count <- mnn_output[[2]]
    shape <- as.integer(c(nrow(data2),ncol(vect)))
    vect_reduced <- np$zeros(shape, dtype = np$float32)
    for (i in 1:length(mnn2)){
        vect_reduced[mnn2[i],] <- vect_reduced[mnn2[i],] + vect[i,]
    }
    vect_avg <- np$divide(vect_reduced[mnn_index + 1,], as.matrix(np$ndarray$astype(mnn_count, np$float32), ncol = 1))
    vect_avg <- as.matrix(vect_avg, ncol = 1)

    dist <- t(apply(data2, 1, FUN = function(x) x%*%t(data2[mnn_index + 1,])))
    exp_distance <- np$exp(-dist/sigma)
    density <- np$sum(exp_distance[mnn_index + 1,], axis = 0L)
    mult <- np$divide(exp_distance, density)
    total_prob <- np$sum(mult, axis = 1L, keepdims = TRUE)
    output_mult <- np$dot(mult, vect_avg)
    correction.in <- np$divide(output_mult, total_prob)
    
    return(correction.in)
}

adjust_s_variance <- function(data1, data2, curcell, curvect, sigma){
    np <- import("numpy")
    shape <- as.integer(c(nrow(data1), 2))
    distance1 <- np$zeros(shape, dtype = np$float32)
    l2_norm <- np$linalg$norm(curvect)
    grad <- np$divide(curvect, l2_norm)
    curproj <- np$dot(grad, curcell)
    prob2 <- 0
    totalprob2 <- 0
    for (samecell in 1:nrow(data2)){
        sameproj <- np$dot(grad, data2[samecell,])
        working <- curcell - data2[samecell,]
        scale <- np$dot(working, grad)
        working <- working - grad * scale
        samedist <- np$dot(working, working)
    #   samedist = mnnpy$utils$sq_dist_to_line(curcell, grad, samecell)
        sameprob <- np$exp(-samedist / sigma)
        if (sameproj <= curproj){
           prob2 <- prob2 + sameprob
        }
        totalprob2 <- totalprob2 + sameprob
    }
    prob2 <- prob2/totalprob2
    totalprob1 <- 0
    for (other in 1:nrow(data1)){
        othercell <- data1[other,]
        distance1[other, 1] <- np$dot(grad, othercell)
        working <- curcell - othercell
        scale <- np$dot(working, grad)
        working <- working - grad * scale
        otherdist <- np$dot(working, working)
        #otherdist = sq_dist_to_line(curcell, grad, othercell)
        weight <- np$exp(-otherdist/sigma)
        distance1[other, 2] <- weight
        totalprob1 <- totalprob1 + weight
    }
    index <- np$argsort(distance1[,1])
    distance1 <- distance1[index + 1,]
    target <- prob2 * totalprob1
    cumulative <- 0
    ref_quan <- distance1[nrow(distance1), 1]
    for (i in 1:nrow(distance1)){
        cumulative <- cumulative + distance1[i,][2]
        if (cumulative > target){
           ref_quan <- distance1[i,][1]
           break
        }
    }
    scaling <- (ref_quan - curproj)/l2_norm

    return(scaling)
}

