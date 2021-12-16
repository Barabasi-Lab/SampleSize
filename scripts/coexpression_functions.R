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

################################################################################

#'  Write correlation networks to file (to avoid memory errors).
#'
#'  Method to calculate correlation between genes from gene expression data.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Spearman ("s", "spearman") or Pearson ("p", "pearson") correlation. Default is "pearson". 
#'
calculate_correlation <- function(rnaseq, out_name, cor_method="pearson"){
  
  correlation_result = CorrelationOverlap(Data = rnaseq, Overlap = row.names(rnaseq), method = "p") %>% as.data.frame() 
  #correlation_result = correlation_result %>% wTO.in.line() %>% rename(!!cor_method := "wTO")
  fwrite(correlation_result, out_name)
  
}


################################################################################

#'  Write correlation networks to file (to avoid memory errors).
#'
#'  Method to calculate correlation between genes from gene expression data.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Default is "pearson".
#'
#calculate_correlation_parallelized <- function(rnaseq, out_name, cor_method="pearson"){
#  write("to\tfrom\tcor\tpval", file=out_name, append=FALSE)
#  res <- foreach(i = seq_len(ncol(rnaseq.t)),
#                 .combine = rbind,
#                 .multicombine = TRUE,
#                 .inorder = FALSE,
#                 .packages = c('data.table', 'doParallel')) %dopar% {
#                   cor(mat[,i], mat, method = 'pearson')
#                 }
#}


################################################################################

#'  Write correlation networks from WGCNA to file.
#'
#'  Method to calculate correlation between genes from gene expression data
#'  using WGCNA.
#'  
#'  @param rnaseq   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the network to.
#'  @param type network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid", "distance".
#'  @param power soft thresholding power.
#'
calculate_correlation_WGCNA <- function(rnaseq, out_name, type="signed", power=6){
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

