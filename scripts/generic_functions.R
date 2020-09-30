#' Creation of counting sparse matrix
#' @description Get a seurat object and return a sparse matrix with the counst
#' @name create_count_matrix
#' @param seurat_object ExpressionSet object for bulk samples
#' @return dgCMatrix object with the count matrix of the seurat object passed as a parameter
#' @export
create_count_matrix <- function(seurat_object){
  counts <-  as(seurat_object@assays$RNA@counts, "dgCMatrix")
  
  # Make the column names as the cell IDs and the row names as the gene IDs
  rownames(counts) <- rownames(seurat_object@assays$RNA@counts)
  colnames(counts) <- seurat_object@meta.data$index_cell_id
  
  counts
}

#' Creation of ESET object
#' @description Base on some parameters the function creates the ESET functions
#' @name getESET
#' @param exprs counting matrix
#' @param fdata Feature data, for genes, usually it's gene name
#' @param pdata Pheno data, for samples, usually it's the characteristics for each single cell/bulk sample, including name, gender, age, cluster, disease,...
#' @return ESET object with counting matrix, samples and additional attributes related with the RNA-seq experiment
#' @export
getESET <- function (exprs, fdata, pdata) 
{
  pdata <- as.data.frame(pdata)
  fdata <- as.data.frame(fdata)
  exprs <- Matrix(exprs, sparse = TRUE)
  rownames(pdata) <- colnames(exprs)
  rownames(fdata) <- rownames(exprs)
  eset <- ExpressionSet(exprs, AnnotatedDataFrame(pdata), 
                        AnnotatedDataFrame(fdata))
}

#' Creation of a sparse ESET object base on a counting matrix and a seurat object
#' @description Get a seurat object and return a sparse matrix with the counst
#' @name create_sparse_eset_object
#' @param count Counting matrix
#' @param seurat_object Seurat object with all attributes related with the RNA-seq experiment: sample, clustering and so on.
#' @return ESET object with counting matrix, samples and additional attributes related with the RNA-seq experiment
#' @export
create_sparse_eset_object <- function(count, seurat_object){
  sparse_matrix <- as(count, "dgCMatrix")
  fdata_split_dataset <- rownames(sparse_matrix)
  pdata_split_dataset <- cbind(cellname = seurat_object@meta.data$index_cell_id, 
                               stage = as.character(seurat_object@meta.data$stage),
                               age = as.character(seurat_object@meta.data$age),
                               sample_id = as.character(seurat_object@meta.data$sample.ID), 
                               cluster = as.character(seurat_object@meta.data$cell_types_name),
                               cluster_normalized = as.character(seurat_object@meta.data$cluster_normalized))
  eset.sc.sparse <- getESET(sparse_matrix, fdata = fdata_split_dataset, pdata = pdata_split_dataset)
  
  eset.sc.sparse
}

#' Addition of the subclustering attribute to the ESET object 
#' @description Get a seurat object and return a sparse matrix with the counst
#' @name add_metacluster_pseudotime
#' @param sc.eset ESET object with clustering
#' @return ESET object with a especific subclustering added.
#' @export
add_metacluster_pseudotime <- function (sc.eset){
  
  sc.eset$metacluster[sc.eset$cluster_normalized %in%  paste0('cluster_',c(8, 9))] <-  "Subcluster 1"
  sc.eset$metacluster[sc.eset$cluster_normalized %in%  paste0('cluster_',c(10, 11, 28))] <-  "Subcluster 2"
  sc.eset$metacluster[sc.eset$cluster_normalized %in%  paste0('cluster_',c(3, 4))] <-  "Subcluster 3"
  
  sc.eset
}

#' Calculation of proportion for the thymus single-cell data.
#' @description Calculate the proportions of cell types for the three time points 4w, 8w, 24w
#' @name calculate_proportions_4w_8w_24w
#' @param sc.eset ESET object with clustering
#' @param last_cluster_number Last number in the cluster normalized
#' @return Proportion for the three time points for each cell type
#' @export
calculate_proportions_4w_8w_24w <- function(sc.eset, last_cluster_number ){
  
  #get cluster number for each cell
  cluster_total_numbers <- as.matrix(sc.eset$cluster_normalized)
  cluster_total_numbers <- t(sapply(cluster_total_numbers, trimws))
  colnames(cluster_total_numbers) <- NULL
  cluster_total_numbers <- as.integer(t(paste0(trimws(str_remove(cluster_total_numbers, "cluster_")))))
  length(cluster_total_numbers)
  summary(cluster_total_numbers)
  
  
  #call the function
  thymus_4w <- calculate_cluster_proportion('4W', sc.eset, cluster_total_numbers, last_cluster_number)
  thymus_8w <- calculate_cluster_proportion('8W', sc.eset, cluster_total_numbers, last_cluster_number)
  thymus_24w <- calculate_cluster_proportion('24W', sc.eset, cluster_total_numbers, last_cluster_number)
  
  join <- rbind('4w_thymus' = thymus_4w[3,], '8w_thymus'=thymus_8w[3,], '24w_thymus'=thymus_24w[3,])
  sum(join)
  join
}

#' Calculation of proportion for a specific cluster for the thymus single-cell data.
#' @description Calculate the proportions of cell types for the three time points 4w, 8w, 24w
#' @name calculate_cluster_proportion
#' @param filter_string Filter to be used to select the set of cells for the cell type
#' @param sc.eset ESET object with clustering
#' @param cluster_total_numbers Number of clusters in the analyisis
#' @param last_cluster_number Last number in the cluster normalized
#' @param verbose For more information of the execution of the function. For default FALSE.
#' @return Proportion for the specific cluster
#' @export
calculate_cluster_proportion <- function(filter_string, sc.eset, cluster_total_numbers, last_cluster_number, verbose = FALSE){
  
  number_cluster_dataset <- last_cluster_number
  cluster_proportion <- data.frame()
  
  #iterate into the number of clusters in order to sum up.
  for(counter in 1:number_cluster_dataset){
    
    if(verbose)
      print(paste0('Cluster No:', counter))
    
    #filter for each cell type in the clustering
    columns_to_sum <- ((trimws(t(sc.eset$age))==filter_string) & t(cluster_total_numbers) == counter)
    summary(t(columns_to_sum))
    
    #matrix with all cell count for each cluster
    counts <- as.matrix(sc.eset@assayData$exprs[,which(columns_to_sum)])
    nrow(columns_to_sum)
    ncol(columns_to_sum)
    
    #sum for the matrix
    total_sum <- sum(counts)
    total_sum
    
    cluster_proportion[1,counter] <- total_sum
  }
  
  sum_cluster <- sum(cluster_proportion[1,])
  
  #Calculate the percentage and the normalized
  for(row_counter in 1:number_cluster_dataset){
    cluster_proportion[2, row_counter] <- ((cluster_proportion[1, row_counter]*100)/sum_cluster)
    cluster_proportion[3, row_counter] <- (cluster_proportion[1, row_counter]/sum_cluster)
  }
  
  rownames(cluster_proportion) <- c('count','percentage','normalized')
  colnames(cluster_proportion) <- paste0('cluster_',seq(1:number_cluster_dataset))
  cluster_proportion
}

#' Sorting of columns by cluster number name (cluster_##)
#' @description Sorting of columns by cluster number name (cluster_##)
#' @name order_result_generic
#' @param matrix Matrix with the data
#' @param start Number in the cluster that is the begining of the counting
#' @param end Last number in the cluster sequence.
#' @param verbose For more information of the execution of the function. For default FALSE.
#' @return Matrix with the sorted columns
#' @export
order_result_generic <- function(matrix, start = 1, end = 29, verbose = FALSE){
  
  matrix_return <- NULL
  names_list <- NULL
  index <- NULL
  to <- end
  
  for(counter in 1:to){
    index <- counter+start-1
    
    if(verbose)
      print(paste0(counter, ", ", index))
    
    if(paste0('cluster_',index) %in% colnames(matrix)){
      matrix_return <- cbind(matrix_return, matrix[,paste0('cluster_',index)])
      names_list <- c(names_list, paste0('cluster_',index))
    }
  }
  colnames(matrix_return) <- names_list
  
  matrix_return
}

#' Applying the proportions to the summary scRNA-seq object
#' @description Applying the proportions to the summary scRNA-seq object
#' @name apply_proportion_celltype_summary
#' @param pseudotime.summary.by.celltype Object with the summary of the sc-RNA-seq summarized by cell type
#' @param sample.proportions Vector with the sample proportions of each cell type
#' @param factor Allows to increment or decrement the weight of the simulated bulk data, 1 means the same weight that the original sc-RNA-seq data. The default value is 1. 
#'               Note: The factor says how scale the counting due to the fact that we are going to round each individual count, therefore some of then are going to be zero because are lesser than zero, 
#'               however the size of the replica or sample is importante, compared with the actual replica, it is about 15M.
#' @param verbose For more information of the execution of the function. For default FALSE.
#' @return Bulk data simulated with the proportions in the parameter sample.proportions
#' @export
apply_proportion_celltype_summary <- function(pseudotime.summary.by.celltype, sample.proportions, factor = 1, verbose = FALSE){
  
  #I have to calculate the proportion of the single cell and then compare that with the proportion that I will try to get, 
  # after I can calculate how should I modify the counting genes in order to have the expected proportion. If I don't do this
  # the final bulk data is going to have always the same proportion that the original sc-RNA-seq data.
  proportions.summary.real.sc <- apply(pseudotime.summary.by.celltype, 2, sum)/sum(pseudotime.summary.by.celltype)
  multiplicator.by.cluster <- sample.proportions / proportions.summary.real.sc
  
  pseudotime.proportional.celltype <- NULL
  for(counter in 1:length(sample.proportions)){
    
    celltype.cluster.name <- names(sample.proportions[counter])
    
    if(verbose)
      print(paste0("celltype name: ", celltype.cluster.name))
    
    #I get the counting of the specific celltype multiplied by the proportion of the sample.
    proportinal.celltype.counting <- pseudotime.summary.by.celltype[, celltype.cluster.name]  * multiplicator.by.cluster[celltype.cluster.name] * factor
    
    pseudotime.proportional.celltype <- cbind(pseudotime.proportional.celltype, as.integer(proportinal.celltype.counting))
  }
  colnames(pseudotime.proportional.celltype) <- names(sample.proportions)
  rownames(pseudotime.proportional.celltype) <- rownames(pseudotime.summary.by.celltype)
  
  if(verbose){
    print(paste0("Real proportions: ", sum(sample.proportions)))
    print(sample.proportions)
  }
  
  total_proportions <- apply(pseudotime.proportional.celltype, 2, sum) / sum(pseudotime.proportional.celltype)
  
  if(verbose){
    print(paste0("bulk proportions (simulations): ", sum(total_proportions)))
    print(total_proportions)
    print(paste0('Correlation: ', cor(sample.proportions, total_proportions)))
  }
  
  #summing up (rows=genes) of the whole set of celltypes
  bulk.data.simulated.sample <- apply(pseudotime.proportional.celltype, 1, sum)
  
  if(verbose)
    print(sum(bulk.data.simulated.sample))
  
  bulk.data.simulated.sample
}

#' Creating a summary object by cell type based on a scRNA-seq object
#' @description Creating a summary object by cell type based on a scRNA-seq object
#' @name create_summary_single_cell
#' @param pseudotime.celltypes List with the cell types of interest (Normalized)
#' @param shared.genes.bulk.single_cell Shared genes between the bulk and single cell data.
#' @param seurat_object.4w_8w_24w Seurat object with the single cell data.
#' @param sc.data.count Sparse matrix of the single cell data with the counting information
#' @param verbose For more information of the execution of the function. For default FALSE
#' @return Matrix with the counting for each cell type as columns and genes as rows.
#' @export
create_summary_single_cell <- function(pseudotime.celltypes, shared.genes.bulk.single_cell, seurat_object, sc.data.count, verbose = FALSE){
  pseudotime.summary.by.celltype = NULL
  for(counter in 1:length(pseudotime.celltypes)){
    
    if(verbose)
      print(paste0('Celltype: ', pseudotime.celltypes[counter]))
    
    celltype.counting <- sc.data.count[shared.genes.bulk.single_cell, seurat_object$cluster_normalized %in% pseudotime.celltypes[counter]]
    
    #sum up by row. We are losing the number of cells here
    celltype.counting.by.gene <- apply(celltype.counting, 1, sum)
    
    pseudotime.summary.by.celltype <- cbind(pseudotime.summary.by.celltype, as.integer(celltype.counting.by.gene))
  }
  colnames(pseudotime.summary.by.celltype) <- pseudotime.celltypes
  rownames(pseudotime.summary.by.celltype) <- shared.genes.bulk.single_cell
  
  pseudotime.summary.by.celltype.proportions <- apply(pseudotime.summary.by.celltype, 2, sum)/sum(pseudotime.summary.by.celltype)
  
  if(verbose){
    print('Proportions:')
    print(pseudotime.summary.by.celltype.proportions)    
  }
  
  pseudotime.summary.by.celltype
}

#' Creation of bulk data applying proportions to the complete set of replicas and samples in the original bulk data.
#' @description Creation of bulk data applying proportions to the complete set of replicas and samples in the original bulk data.
#' @name apply_proportions_complete_single_cell_dataset
#' @param bulk.data.count Matrix with the count expression of the bulk data.
#' @param deconvolution.real_data Result of the deconvolution with the calculated proportions.
#' @param pseudotime.summary.by.celltype Object with the summary of counting for each cell type in the single cell data.
#' @param verbose For more information of the execution of the function. For default FALSE
#' @return Matrix with the counting summarized for each sample simulating a bulk data.
#' @export
apply_proportions_complete_single_cell_dataset <- function(bulk.data.count, deconvolution.real_data, pseudotime.summary.by.celltype, verbose = FALSE){
  sample.names.list <- colnames(bulk.data.count)
  bulk.simulated.thymus.pseudotime <- NULL
  for(counter in 1:length(sample.names.list)){
    sample.name <- sample.names.list[counter]
    
    if(verbose)
      print(paste0("Creating bulk data for sample: ", sample.name))
    
    sample.proportions <- deconvolution.real_data$prop.est[sample.name,]
    
    #by using the function, I get each bulk data for each sample/replica
    bulk.data.simulated.replica <- apply_proportion_celltype_summary(pseudotime.summary.by.celltype, sample.proportions, 100)
    bulk.simulated.thymus.pseudotime <- cbind(bulk.simulated.thymus.pseudotime, bulk.data.simulated.replica)
  }
  colnames(bulk.simulated.thymus.pseudotime) <- paste0("simu_", sample.names.list)
  
  bulk.simulated.thymus.pseudotime
}

#
#' @title Compute TPM for a read count matrix. Source: https://www.biostars.org/p/335187/
#' @param dfr A numeric data.frame of read counts with samples (columns) and genes (rows).
#' @param len A vector of gene cds length equal to number of rows of dfr.
#' 
r_tpm <- function(dfr,len)
{
  dfr1 <- sweep(dfr,MARGIN=1,(len/10^4),`/`)
  scf <- colSums(dfr1)/(10^6)
  return(sweep(dfr1,2,scf,`/`))
}

#' Creation of a differential expression analysis usign DESeq2 for all simulated samples vs a real sample
#' @description Creation of a differential expression analysis usign DESeq2 for all simulated samples vs a real sample
#' @name create_differential_expression_set
#' @param reference.sample.name Reference name of the differential expression analysis: Thy_WT, Thy_KO, Thy_ALKpos_WT, Thy_ALKpos_KO
#' @param selected.genes.TPM Genes selected as top 10%
#' @param bulk.data.simulated.plus.sample Bulk data with all simulated samples and the reference real sample
#' @return Object with the different combination of samples based on the experiment configuration
#' @export
create_differential_expression_set <- function(reference.sample.name, selected.genes.TPM, bulk.data.simulated.plus.sample){
  
  countdata <- bulk.data.simulated.plus.sample[rownames(bulk.data.simulated.plus.sample) %in% selected.genes.TPM, ]
  
  #I put numbers on it because other wise I wouldn't get  ko vs wt.
  type <- c("simulation","simulation","simulation",
            "simulation","simulation","simulation",
            "simulation","simulation","simulation",
            "simulation","simulation","simulation",
            "real_tumor","real_tumor","real_tumor")
  
  if(reference.sample.name=='Thy_ALKpos_WT'){
    treatment <- c("simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO",
                   "simu_Thy_KO","simu_Thy_KO","simu_Thy_KO",
                   "simu_Thy_WT","simu_Thy_WT","simu_Thy_WT",
                   "simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT",
                   "Thy_ALKpos_WT","Thy_ALKpos_WT","Thy_ALKpos_WT")
  }else if(reference.sample.name=='Thy_WT'){
    treatment <- c("simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO",
                   "simu_Thy_KO","simu_Thy_KO","simu_Thy_KO",
                   "simu_Thy_WT","simu_Thy_WT","simu_Thy_WT",
                   "simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT",
                   "Thy_WT","Thy_WT","Thy_WT")
    
  }else if(reference.sample.name=='Thy_KO'){
    treatment <- c("simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO",
                   "simu_Thy_KO","simu_Thy_KO","simu_Thy_KO",
                   "simu_Thy_WT","simu_Thy_WT","simu_Thy_WT",
                   "simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT",
                   "Thy_KO","Thy_KO","Thy_KO")
    
  }else if(reference.sample.name=='Thy_ALKpos_KO'){
    treatment <- c("simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO","simu_Thy_ALKpos_KO",
                   "simu_Thy_KO","simu_Thy_KO","simu_Thy_KO",
                   "simu_Thy_WT","simu_Thy_WT","simu_Thy_WT",
                   "simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT","simu_Thy_ALKpos_WT",
                   "Thy_ALKpos_KO","Thy_ALKpos_KO","Thy_ALKpos_KO")
    
  }
  
  time <- c("18w","18w","18w","18w","18w","18w",
            "18w","18w","18w","18w","18w","18w",
            "18w","18w","18w")
  
  replicates  <- c("1","2","3","1","2","3",
                   "1","2","3","1","2","3",
                   "1","2","3")
  
  #secondary information
  sample <- paste0("count_", 1:15)
  
  coldata <- as.data.frame(cbind(replicates, treatment, time, sample, type))
  rownames(coldata) <- sample
  
  #building the object DESeqDataSet
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = countdata,
    colData = coldata,
    design = ~treatment)
  
  if(reference.sample.name=='Thy_ALKpos_WT'){
    #giving the order and setting the control sample for the analyisis, which is - in this case - the real Tumor_ALKpos_WT.
    ddsFullCountTable$treatment <- factor(ddsFullCountTable$treatment, levels = c("Thy_ALKpos_WT","simu_Thy_WT", "simu_Thy_KO", "simu_Thy_ALKpos_WT", "simu_Thy_ALKpos_KO"))
    ddsFullCountTable$treatment <- relevel(ddsFullCountTable$treatment, ref = "Thy_ALKpos_WT")
  }else if(reference.sample.name=='Thy_WT'){
    ddsFullCountTable$treatment <- factor(ddsFullCountTable$treatment, levels = c("Thy_WT","simu_Thy_WT", "simu_Thy_KO", "simu_Thy_ALKpos_WT", "simu_Thy_ALKpos_KO"))
    ddsFullCountTable$treatment <- relevel(ddsFullCountTable$treatment, ref = "Thy_WT")
  }else if(reference.sample.name=='Thy_KO'){
    ddsFullCountTable$treatment <- factor(ddsFullCountTable$treatment, levels = c("Thy_KO","simu_Thy_WT", "simu_Thy_KO", "simu_Thy_ALKpos_WT", "simu_Thy_ALKpos_KO"))
    ddsFullCountTable$treatment <- relevel(ddsFullCountTable$treatment, ref = "Thy_KO")
  }else if(reference.sample.name=='Thy_ALKpos_KO'){
    ddsFullCountTable$treatment <- factor(ddsFullCountTable$treatment, levels = c("Thy_ALKpos_KO","simu_Thy_WT", "simu_Thy_KO", "simu_Thy_ALKpos_WT", "simu_Thy_ALKpos_KO"))
    ddsFullCountTable$treatment <- relevel(ddsFullCountTable$treatment, ref = "Thy_ALKpos_KO")  
  }
  
  dds <- ddsFullCountTable
  
  ## Pre-filter
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  
  ## The rlog and variance stabilizing transformations => homoskedastic data
  rld <- rlog(dds, blind = FALSE)
  vsd <- vst(dds, blind = FALSE)
  dds <- estimateSizeFactors(dds)
  normalized.wt.bulk <- counts(dds, normalized=TRUE)
  
  dds <- DESeq(dds)
  
  if(reference.sample.name=='Thy_ALKpos_WT'){
    RES.thymus.lfc <- list(treatment_simu_Thy_WT_vs_Thy_ALKpos_WT = results(dds, name = "treatment_simu_Thy_WT_vs_Thy_ALKpos_WT"), 
                           treatment_simu_Thy_KO_vs_Thy_ALKpos_WT = results(dds, name = "treatment_simu_Thy_KO_vs_Thy_ALKpos_WT"), 
                           treatment_simu_Thy_ALKpos_WT_vs_Thy_ALKpos_WT = results(dds, name = "treatment_simu_Thy_ALKpos_WT_vs_Thy_ALKpos_WT"), 
                           treatment_simu_Thy_ALKpos_KO_vs_Thy_ALKpos_WT = results(dds, name = "treatment_simu_Thy_ALKpos_KO_vs_Thy_ALKpos_WT"))
  }else if(reference.sample.name=='Thy_WT'){
    RES.thymus.lfc <- list(treatment_simu_Thy_WT_vs_Thy_WT = results(dds, name = "treatment_simu_Thy_WT_vs_Thy_WT"), 
                           treatment_simu_Thy_KO_vs_Thy_WT = results(dds, name = "treatment_simu_Thy_KO_vs_Thy_WT"), 
                           treatment_simu_Thy_ALKpos_WT_vs_Thy_WT = results(dds, name = "treatment_simu_Thy_ALKpos_WT_vs_Thy_WT"), 
                           treatment_simu_Thy_ALKpos_KO_vs_Thy_WT = results(dds, name = "treatment_simu_Thy_ALKpos_KO_vs_Thy_WT"))
  }else if(reference.sample.name=='Thy_KO'){
    RES.thymus.lfc <- list(treatment_simu_Thy_WT_vs_Thy_KO = results(dds, name = "treatment_simu_Thy_WT_vs_Thy_KO"), 
                           treatment_simu_Thy_KO_vs_Thy_KO = results(dds, name = "treatment_simu_Thy_KO_vs_Thy_KO"), 
                           treatment_simu_Thy_ALKpos_WT_vs_Thy_KO = results(dds, name = "treatment_simu_Thy_ALKpos_WT_vs_Thy_KO"), 
                           treatment_simu_Thy_ALKpos_KO_vs_Thy_KO = results(dds, name = "treatment_simu_Thy_ALKpos_KO_vs_Thy_KO"))
  }else if(reference.sample.name=='Thy_ALKpos_KO'){
    
    RES.thymus.lfc <- list(treatment_simu_Thy_WT_vs_Thy_ALKpos_KO = results(dds, name = "treatment_simu_Thy_WT_vs_Thy_ALKpos_KO"), 
                           treatment_simu_Thy_KO_vs_Thy_ALKpos_KO = results(dds, name = "treatment_simu_Thy_KO_vs_Thy_ALKpos_KO"), 
                           treatment_simu_Thy_ALKpos_WT_vs_Thy_ALKpos_KO = results(dds, name = "treatment_simu_Thy_ALKpos_WT_vs_Thy_ALKpos_KO"), 
                           treatment_simu_Thy_ALKpos_KO_vs_Thy_ALKpos_KO = results(dds, name = "treatment_simu_Thy_ALKpos_KO_vs_Thy_ALKpos_KO"))
  }
  
  RES.thymus.lfc
}
