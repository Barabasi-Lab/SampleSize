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
  correlation_result = correlation_result %>% wTO.in.line() %>% rename(!!cor_method := "wTO")
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
  
  adjacency.res = adjacency(rnaseq, type=type, power = power)
  adjacency.res = adjacency.res %>% wTO.in.line() %>% rename(WGCNA=wTO)
  fwrite(adjacency.res, out_name)
  
}

