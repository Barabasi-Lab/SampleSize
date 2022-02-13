################################################################################
# This script contains the functions needed to make the analyses with
# co-expression networks.

# Author: Joaquim Aguirre Plans
# Date: October 2021
################################################################################

require(WGCNA)
require(wTO)
require(dplyr)
require(data.table)
require(Hmisc)
require(igraph)
require(minet)


################################################################################

#'  Calculate disparity filter algorithm from MA Sierra et al.
#'
#'  Method to calculate the disparity filter to remove spurious correlations.
#'  
#'  Function extracted from https://github.com/menchelab/MultiOme/blob/main/functions/sim_mat_functions.R
#'  
#'  @param weighted_adj_mat   Weighted adjacency matrix of correlation values.
#'
backbone_cal = function(weighted_adj_mat){
  pval_mat = weighted_adj_mat
  #====
  message("Calculate weight and degree of all nodes")
  W = colSums(weighted_adj_mat, na.rm = T)
  k = apply(weighted_adj_mat, 2, function(x) sum(!is.na(x)))
  #====
  message("Calculate p-value from each connection")
  for(i in 1:ncol(pval_mat)){
    pval_mat[,i] = (1-(weighted_adj_mat[,i]/W[i]))^(k[i]-1)
  }
  return(pval_mat)
}


################################################################################

#'  Write correlation networks to file.
#'
#'  Method to calculate correlation between genes from gene expression data.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Spearman ("s", "spearman") or Pearson ("p", "pearson") correlation. Default is "pearson". 
#'
calculate_correlation <- function(rnaseq, out_name, cor_method="pearson", disparity_filter=TRUE, corr_threshold=NA, pval_threshold=NA){

  correlation_result = CorrelationOverlap(Data = rnaseq, Overlap = row.names(rnaseq), method = "p") %>% as.data.frame() 
  
  if (!(is.na(corr_threshold))){
    correlation_result[abs(correlation_result) <= corr_threshold] = NA
  }

  if (isTRUE(disparity_filter)){
    correlation_pval_mat = backbone_cal(correlation_result)
    correlation_result = correlation_result %>% wTO.in.line() %>% rename(!!cor_method := "wTO")
    correlation_pval_mat = correlation_pval_mat %>% wTO.in.line() %>% rename("pval" := "wTO")
    correlation_result = correlation_result %>% inner_join(correlation_pval_mat)
    rm(correlation_pval_mat)
    if (!(is.na(pval_threshold))){
      correlation_result = correlation_result %>% filter(pval < pval_threshold)
    }
  }
  
  fwrite(correlation_result, out_name)
  
}


################################################################################

#'  Write correlation networks to file and includes calculation of p-value.
#'
#'  Method to calculate correlation and p-value between genes from gene expression data.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Spearman ("s", "spearman") or Pearson ("p", "pearson") correlation. Default is "pearson". 
#'  @param correction_method Method to correct p-value for multiple testing. Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr". Default is "bonferroni". 
#'
calculate_correlation_and_pvalue <- function(rnaseq, out_name, cor_method="spearman", correction_method="bonferroni", disparity_filter=FALSE, disparity_filter_pval_threshold=0.05){
  
  correlation_result <- rcorr(as.matrix(rnaseq), type=cor_method)
  correlation_result_mat = correlation_result$r
  correlation_pval_mat = correlation_result$P %>% wTO.in.line() %>% rename("pval" := "wTO")
  correlation_result = correlation_result$r %>% wTO.in.line() %>% rename("score" := "wTO")
  correlation_result = correlation_result %>% inner_join(correlation_pval_mat)
  rm(correlation_pval_mat)
  correlation_result$pval.adj = p.adjust(correlation_result$pval, method=correction_method)
  
  if (isTRUE(disparity_filter)){
    # Calculate disparity p-value without filtering
    correlation_disp_mat = backbone_cal(correlation_result_mat)
    correlation_disp_mat = correlation_disp_mat %>% wTO.in.line() %>% rename("disp.pval" := "wTO")
    correlation_result = correlation_result %>% left_join(correlation_disp_mat)
    rm(correlation_disp_mat)
    rm(correlation_result_mat)
    correlation_result_filt = correlation_result
    correlation_result_filt[correlation_result_filt$pval.adj >= disparity_filter_pval_threshold,]$score = NA
    correlation_result_filt = correlation_result_filt %>% rename("weight" := "score") %>% dplyr::select(Node.1, Node.2, weight)
    correlation_result_mat = graph_from_data_frame(correlation_result_filt)
    correlation_result_mat = as.data.frame(as.matrix(as_adjacency_matrix(correlation_result_mat, attr="weight")))
    correlation_disp_mat = backbone_cal(correlation_result_mat)
    rm(correlation_result_mat)
    correlation_disp_mat = correlation_disp_mat %>% wTO.in.line() %>% rename(!!paste("disp.pval.after.filt.",disparity_filter_pval_threshold,sep="") := "wTO")
    correlation_result = correlation_result %>% left_join(correlation_disp_mat)
    rm(correlation_disp_mat)
  }
  #correlation_result = correlation_result %>% rename(!!cor_method := "score")
  
  fwrite(correlation_result, out_name)
  
}


################################################################################

#'  Write mutual information matrix to file.
#'
#'  Method to calculate mutual information between genes from gene expression data and write a mutual information matrix in the output file.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are genes, rows are samples).
#'  @param out_name The name of the file to write the correlations to.
#'  @param estimator Name of the entropy estimator to be used. The package can use the four mutual information estimators implemented in the package "infotheo": "mi.empirical", "mi.mm", "mi.shrink", "mi.sg" and three estimators based on correlation: "pearson","spearman","kendall"(default:"spearman") - see details.
#'  "mi.empirical" : This estimator computes the entropy of the empirical probability distribution.
#'  "mi.mm" : This is the Miller-Madow asymptotic bias corrected empirical estimator.
#'  "mi.shrink" : This is a shrinkage estimate of the entropy of a Dirichlet probability distribution.
#'  "mi.sg" : This is the Schurmann-Grassberger estimate of the entropy of a Dirichlet probability distribution.
#'  "pearson" : This computes mutual information for normally distributed variable.
#'  "spearman" : This computes mutual information for normally distributed variable using Spearman's correlation instead of Pearson's correlation.
#'  "kendall" : This computes mutual information for normally distributed variable using Kendall's correlation instead of Pearson's correlation.
#'
calculate_mutual_information <- function(rnaseq, out_name, estimator="pearson"){
  
  mim <- build.mim(rnaseq, estimator = estimator)
  fwrite(mim, out_name)
  
}

################################################################################

#'  Create gene co-expression networks using WGCNA.
#'
#'  Method to create create co-expression networks using WGCNA.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are genes, rows are samples).
#'  @param out_name The name of the file to write the network to.
#'  @param type network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid", "distance".
#'  @param power soft thresholding power.
#'
calculate_network_WGCNA <- function(rnaseq, out_name, type="signed", power=6){
  # Information about differences between Adjacency matrix and TOM:
  # https://www.researchgate.net/post/What_do_adjacency_matrix_and_Topology_Overlap_Matrix_from_WGCNA_package_tell_about_the_data
  # WGCNA constructs two matrices, first it defines a correlation matrix up to a power beta so the degree distribution will fit a small-word network.
  # This matrix given only information about the expression correlation between genes.
  # WGCNA thinks that co-expression is not enough and the similarity between genes should be reflected at the expression and the network topology level.
  # This is why it defines the TOM matrix which uses the co-expression Adjacency matrix and build another adjacency matrix that considers topological similarity.
  # Normally the TOM matrix is the final result of WGCNA.
  
  adjacency.res = adjacency(rnaseq, type=type, power = power)
  TOM = TOMsimilarity(adjacency.res, TOMType = type)
  rownames(TOM) = rownames(adjacency.res)
  colnames(TOM) = colnames(adjacency.res)
  #TOM = TOM %>% wTO.in.line() %>% rename(WGCNA=wTO)
  fwrite(TOM, out_name)
  
}

################################################################################

#'  Create gene co-expression networks using ARACNE
#'
#'  Method to create co-expression networks using ARACNE.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are genes, rows are samples).
#'  @param out_name The name of the file to write the network to.
#'  @param estimator Name of the entropy estimator to be used. The package can use the four mutual information estimators implemented in the package "infotheo": "mi.empirical", "mi.mm", "mi.shrink", "mi.sg" and three estimators based on correlation: "pearson","spearman","kendall"(default:"spearman") - see details.
#'  @param eps Numeric value indicating the threshold used when removing an edge : for each triplet of nodes (i,j,k), the weakest edge, say (ij), is removed if its weight is below min{(ik),(jk)}-eps - see references.
#'
calculate_network_ARACNE <- function(rnaseq, out_name, estimator="pearson", eps=0){
  # Information about how to calculate a network using ARACNE with the minet package:
  # https://www.biostars.org/p/123084/
  # First, we use the function build.mim to calculate the mutual information matrix:
  # mim <- build.mim(dataset, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset)))
  # 'dataset' is a matrix where each row corresponds to a sample and each column corresponds to gene.
  # 'estimator' is the name of the entropy estimator to be used. The package can use the four mutual information estimators implemented in the package "infotheo": "mi.empirical", "mi.mm", "mi.shrink", "mi.sg" and three estimators based on correlation: "pearson","spearman","kendall"(default:"spearman") - see details.
  #
  # Then, we use the function aracne, which takes the mutual information matrix as input in order to return the infered network according to the Aracne algorithm. 
  # This algorithm applies the data processing inequality to all triplets of nodes in order to remove the least significant edge in each triplet.
  # net <- aracne(mim,eps=0)
  # 'eps' is a numeric value indicating the threshold used when removing an edge : for each triplet of nodes (i,j,k), the weakest edge, say (ij), is removed if its weight is below min{(ik),(jk)}-eps - see references.
  mim <- build.mim(rnaseq, estimator = estimator)
  net <- aracne(mim,eps=eps)
  fwrite(net, out_name)
  
}

