################################################################################
# This script contains the functions needed to make the analyses with
# co-expression networks.

# Author: Joaquim Aguirre Plans
# Date: October 2021
################################################################################

library(WGCNA)

################################################################################

#'  Write correlation networks to file (to avoid memory errors).
#'
#'  Method to calculate correlation between genes from gene expression data.
#'  
#'  Based on the function calculate_correlation from K. Ovens in: 
#'  https://github.com/klovens/compare/blob/master/coexpression_comparison.R
#'  
#'  @param rnaseq_df   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Default is "pearson".
#'
calculate_correlation <- function(rnaseq_df, out_name, cor_method="pearson"){
  write("to\tfrom\tcor\tpval", file=out_name, append=FALSE)
  for(i in 1:nrow(rnaseq_df)){
    for(j in i:nrow(rnaseq_df)){
      # do not require edge between the same gene (uninformative)
      if (i == j){
        c <- 0
        pval <- 1
      }
      else{
        c.res <- cor.test(x=as.numeric(rnaseq_df[i,]), y=as.numeric(rnaseq_df[j,]),method=cor_method, exact=FALSE)
        c <- as.numeric(c.res$estimate)
        pval <- c.res$p.value
        #c <- cor(x=as.numeric(rnaseq_df[i,]), y=as.numeric(rnaseq_df[j,]),method=cor_method,  use="complete.obs")
        if(is.na(c)){
          c <- 1
          pval <- 0
        }
        write(paste(rownames(rnaseq_df)[i],rownames(rnaseq_df)[j],c, pval,sep="\t"), file=out_name, append=TRUE)
      }
    }
  }
}


################################################################################

#'  Write correlation networks to file (to avoid memory errors).
#'
#'  Method to calculate correlation between genes from gene expression data.
#'  
#'  Based on the function calculate_correlation from K. Ovens in: 
#'  https://github.com/klovens/compare/blob/master/coexpression_comparison.R
#'  
#'  @param rnaseq_df   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Default is "pearson".
#'
#calculate_correlation_parallelized <- function(rnaseq_df, out_name, cor_method="pearson"){
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
#'  Based on the function calculate_correlation from K. Ovens in: 
#'  https://github.com/klovens/compare/blob/master/coexpression_comparison.R
#'  
#'  @param rnaseq_df   RNAseq expression data frame
#'                     (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Default is "pearson".
#'
calculate_correlation_WGCNA <- function(rnaseq_df, out_name, type="signed", power=12){
  
  adjacency.res = adjacency(rnaseq_df, type=type, power = power)
  
  write("to\tfrom\tcor", file=out_name, append=FALSE)
  for(i in 1:ncol(adjacency.res)){
    for(j in i:nrow(adjacency.res)){
      # do not require edge between the same gene (uninformative)
      if (i == j){
        c <- 0
        #pval <- 1
      }
      else{
        c <- adjacency.res[i,j]
        #pval <- c.res$p.value
        if(is.na(c)){
          c <- 1
          #pval <- 0
        }
        write(paste(colnames(adjacency.res)[i],rownames(adjacency.res)[j],c,sep="\t"), file=out_name, append=TRUE)
      }
    }
  }
}
