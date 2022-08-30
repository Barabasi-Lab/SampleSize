#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(wTO)
library(CoDiNA)
set.seed(1510)



option_list = list(
  make_option(c("-d", "--coexpression_network_file_D"), action="store", type="character", 
              help="Co-expression network file of condition 1 (e.g., disease)", metavar="character"),
  make_option(c("-n", "--coexpression_network_file_N"), action="store", type="character", 
              help="Co-expression network file of condition 2  (e.g., normal)", metavar="character"),
  make_option(c("-l", "--output_edges_file"), action="store", type="character", default="./results_diff_net.txt", 
              help="Output results of the differentially co-expressed edges [default= %default]", metavar="character"),
  make_option(c("-v", "--output_nodes_file"), action="store", type="character", default="./results_diff_nodes.txt", 
              help="Output results of the differentially co-expressed nodes [default= %default]", metavar="character"),
  make_option(c("-t", "--threshold"), action="store", type="double", default = 0.05,
              help="P-value threshold", metavar="double"),
  make_option(c("-s", "--stretch_normalization"), action="store_true",
              help="If TRUE, uses stretch parameter to normalize networks before running CODINA"),
  make_option(c("-f", "--filter_by_common_nodes"), action="store_true",
              help="If TRUE, filters networks by keeping nodes that are common in both networks")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_differentially_coexpressed_genes.R -d /work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_tcga/tumor/TCGA-BRCA/pearson_tcga_TCGA-BRCA_size_440_rep_5.net -n /work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_gtex/Breast.Mammary.Tissue/pearson_RNAseq_samples_Breast.Mammary.Tissue_size_440_rep_5.net -l /work/ccnr/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/TCGA-BRCA_Breast.Mammary.Tissue/diffanalysis_edges_pearson_tcga_TCGA-BRCA_size_440_rep_5.net___pearson_RNAseq_samples_Breast.Mammary.Tissue_size_440_rep_5.net_pval_0.05.txt -v /work/ccnr/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/TCGA-BRCA_Breast.Mammary.Tissue/diffanalysis_nodes_pearson_tcga_TCGA-BRCA_size_440_rep_5.net___pearson_RNAseq_samples_Breast.Mammary.Tissue_size_440_rep_5.net_pval_0.05.txt -t 0.05 -s -f 



#----------------#
# READ ARGUMENTS #
#----------------#

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$coexpression_network_file_D) | is.null(opt$coexpression_network_file_N)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coexpression_network_file_D).n", call.=FALSE)
}

coexpression_network_file_D = opt$coexpression_network_file_D
coexpression_network_file_N = opt$coexpression_network_file_N
output_edges_file = opt$output_edges_file
output_nodes_file = opt$output_nodes_file
threshold = as.double(opt$threshold)
stretch_normalization = opt$stretch_normalization
filter_by_common_nodes = opt$filter_by_common_nodes

#coexpression_network_file_D = "/work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_tcga/tumor/TCGA-BRCA/pearson_tcga_TCGA-BRCA_size_200_rep_1.net"
#coexpression_network_file_N = "/work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_gtex/Breast.Mammary.Tissue/pearson_RNAseq_samples_Breast.Mammary.Tissue_size_200_rep_1.net"
#threshold=0.05
#stretch_normalization=TRUE
#filter_by_common_nodes=TRUE



#-------------------------#
# DEFINITION OF FUNCTIONS #
#-------------------------#

###############################
## read_coexpression_network ##
###############################
#'  @param coexpression_network_file The path to the co-expression network file
#'  @param method The method used to create the co-expression network. 
#'  @param threshold P-value threshold to filter the co-expression network. 

read_coexpression_network <- function(coexpression_network_file, method, threshold, set_unsignificant_scores_to_zero=TRUE){
  
  # Read file
  coexpression_df = as.data.frame(fread(coexpression_network_file))
  
  # If method is aracne or wgcna, transform matrix to dataframe of pairs of genes
  if((method == "aracne") | (method == "wgcna")){
    rownames(coexpression_df) = colnames(coexpression_df)
    coexpression_df = coexpression_df[order(rownames(coexpression_df)), order(colnames(coexpression_df))]
    coexpression_df = coexpression_df %>% wTO.in.line() %>% rename(score=wTO)
    if(method == "aracne"){
      threshold = 0
    }
    # Filter network by score
    if(set_unsignificant_scores_to_zero == TRUE){
      # If set_unsignificant_scores_to_zero is TRUE, set scores below threshold to 0
      coexpression_df$score = ifelse(abs(coexpression_df$score) <= threshold, 0, coexpression_df$score)
    } else {
      # If set_unsignificant_scores_to_zero is FALSE, filter out scores below threshold
      coexpression_df = coexpression_df %>% filter(abs(score) > threshold)
    }
  } else {
    if(method == "wto"){
      coexpression_df = coexpression_df %>% rename("score"= "wTO")
    }
    # Filter network by p-value
    if(set_unsignificant_scores_to_zero == TRUE){
      # If set_unsignificant_scores_to_zero is TRUE, set scores of unsignificant edges to 0
      coexpression_df = coexpression_df %>% dplyr::select(Node.1, Node.2, score, pval.adj)
      coexpression_df$score = ifelse(coexpression_df$pval.adj >= threshold, 0, coexpression_df$score)
    } else {
      # If set_unsignificant_scores_to_zero is FALSE, filter out unsignificant edges
      coexpression_df = coexpression_df %>% dplyr::select(Node.1, Node.2, score, pval.adj) %>% filter(pval.adj < threshold)
    }
  }
  
  return(coexpression_df)
}  

###################################
## filter_network_by_correlation ##
###################################
#'  @param network_df Co-expression network in dataframe format
#'  @param upper_threshold Upper threshold of correlation (absolute). If NA, no upper limit is defined.
#'  @param lower_threshold Lower threshold of correlation (absolute).

filter_network_by_correlation <- function(network_df, upper_threshold=NA, lower_threshold=0){
  # Filter network dataframe  
  if(is.na(upper_threshold)){
    network_filt_df = network_df %>% filter(abs(score) >= abs(lower_threshold))
  }else{
    network_filt_df = network_df %>% filter((abs(score) >= abs(lower_threshold)) & (abs(score) < abs(upper_threshold)))
  }
  return(network_filt_df)
}



#--------------#
# READ NETWORK #
#--------------#

# Find file names and method names
file_name_D = tail(strsplit(coexpression_network_file_D, split="/")[[1]], n=1)
file_name_N = tail(strsplit(coexpression_network_file_N, split="/")[[1]], n=1)
method_D = strsplit(file_name_D, split="_")[[1]][1]
method_N = strsplit(file_name_N, split="_")[[1]][1]

# Read networks and filter them by p-value
start_time <- Sys.time()
print("Reading co-expression networks")
coexpression_df_D = read_coexpression_network(coexpression_network_file=coexpression_network_file_D, method=method_D, threshold=threshold, set_unsignificant_scores_to_zero=T)
coexpression_df_N = read_coexpression_network(coexpression_network_file=coexpression_network_file_N, method=method_N, threshold=threshold, set_unsignificant_scores_to_zero=T)
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

#coexpression_df_D_dum = coexpression_df_D[1:10000,]
#coexpression_df_N_dum = coexpression_df_N[1:10000,]

# Filter genes that are not in both networks (if user wants)
if(filter_by_common_nodes == TRUE){
  nodes_D = intersect(unique(coexpression_df_D$Node.1), unique(coexpression_df_D$Node.2))
  nodes_N = intersect(unique(coexpression_df_N$Node.1), unique(coexpression_df_N$Node.2))
  common_nodes = intersect(nodes_D, nodes_N)
  rm(nodes_D)
  rm(nodes_N)
  coexpression_df_D = coexpression_df_D %>% filter((Node.1 %in% common_nodes) & (Node.2 %in% common_nodes))
  coexpression_df_N = coexpression_df_N %>% filter((Node.1 %in% common_nodes) & (Node.2 %in% common_nodes))
  rm(common_nodes)
}



#-----------------------------------------------#
# CALCULATE DIFFERENTIAL CO-EXPRESSION ANALYSIS #
#-----------------------------------------------#

# Calculate differential co-expression analysis
start_time <- Sys.time()
print("Calculating differentially co-expressed edges")
DiffNet = tryCatch(MakeDiffNet(Data = list(coexpression_df_N, coexpression_df_D), # The first network is the baseline (i.e., normal), the second the one that differs (i.e., disease)
                               Code = c("N", "D"), 
                               stretch = stretch_normalization), # stretch is used when the networks come from different studies
                   error=function(cond){
                     print("Not enough nodes to make differential co-expression analysis")
                     # Empty table
                     cols = c("Node.1", "Node.2", "D", "N", "Phi", "Phi_tilde", "Group", "Score_center", "Score_Phi", "Score_Phi_tilde", "Score_internal", "Score_ratio")
                     DiffNet = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
                     return(DiffNet)
                    })
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# Select only strong and well classified edges
DiffNet_cl = DiffNet %>%
  filter(Score_ratio > 1)

# Categorize the Nodes into different roles
start_time <- Sys.time()
print("Calculating differentially co-expressed nodes")
Nodes = tryCatch(ClusterNodes(DiffNet = DiffNet_cl, cutoff.external = 0, cutoff.internal = 1),
                 error=function(cond){
                   print("Not enough nodes to cluster")
                   # Empty table
                   cols = c("Node", "Phi", "pval_Phi", "Group_pval_Phi", "Phi_tilde", "pval_Phi_Tilde", "Group_pval_Phi_tilde")
                   Nodes = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
                   return(Nodes)
                 })
# Correct p-value for multiple testing
Nodes$pval_Phi_Tilde.adj.bonf = p.adjust(Nodes$pval_Phi_Tilde, method="bonferroni")
Nodes$pval_Phi_Tilde.adj.fdr = p.adjust(Nodes$pval_Phi_Tilde, method="fdr")
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# Write final results
DiffNet_cl %>% as.data.frame() %>% fwrite(output_edges_file)
Nodes %>% as.data.frame() %>% fwrite(output_nodes_file)



#----------------------------------------------#
# CALCULATE RESULTS FOR DIFFERENT CORRELATIONS #
#----------------------------------------------#

# # Calculate results for different levels of correlation
# type_correlation_df = data.frame(row.names=c("all", "very weak", "weak", "moderate", "strong", "very strong", "strong-very strong", "moderate-strong-very strong", "weak-moderate-strong-very strong"), lower_val=c(NA,0,0.2,0.4,0.6,0.8,0.6,0.4,0.2), upper_val=c(NA,0.2,0.4,0.6,0.8,NA,NA,NA,NA))
# cols = c("Node.1", "Node.2", "D", "N", "Phi", "Phi_tilde", "Group", "Score_center", "Score_Phi", "Score_Phi_tilde", "Score_internal", "Score_ratio", "type_correlation")
# DiffNet_cl_all = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
# cols = c("Node", "Phi", "pval_Phi", "Group_pval_Phi", "Phi_tilde", "pval_Phi_Tilde", "Group_pval_Phi_tilde", "pval_Phi_Tilde.adj.bonf", "pval_Phi_Tilde.adj.fdr", "type_correlation")
# Nodes_all = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
# 
# for(type_correlation_selected in row.names(type_correlation_df)){
# 
#   print(paste("Correlation type: ", type_correlation_selected, sep=""))
# 
#   #-------------------------------#
#   # FILTER NETWORK BY CORRELATION #
#   #-------------------------------#
# 
#   # Filter network by correlation
#   if(type_correlation_selected == "all"){
#     # same data frame without filtering
#     coexpression_filt_df_D = data.frame(coexpression_df_D)
#     coexpression_filt_df_N = data.frame(coexpression_df_N)
#   } else {
#     coexpression_filt_df_D = filter_network_by_correlation(network_df=coexpression_df_D, upper_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation_selected,]$upper_val, lower_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation_selected,]$lower_val)
#     coexpression_filt_df_N = filter_network_by_correlation(network_df=coexpression_df_N, upper_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation_selected,]$upper_val, lower_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation_selected,]$lower_val)
#   }
# 
#   #-----------------------------------------------#
#   # CALCULATE DIFFERENTIAL CO-EXPRESSION ANALYSIS #
#   #-----------------------------------------------#
# 
#   if((nrow(coexpression_filt_df_D) > 0) & (nrow(coexpression_filt_df_N) > 0)){
#     # Calculate differential co-expression analysis
#     start_time <- Sys.time()
#     print("Calculating differentially co-expressed edges")
#     DiffNet = tryCatch(MakeDiffNet(Data = list(coexpression_filt_df_N, coexpression_filt_df_D), # The first network is the baseline (i.e., normal), the second the one that differs (i.e., disease)
#                                    Code = c("N", "D"), 
#                                    stretch = stretch_normalization), # stretch is used when the networks come from different studies
#                        error=function(cond){
#                          print("Not enough nodes to make differential co-expression analysis")
#                          # Empty table
#                          cols = c("Node.1", "Node.2", "D", "N", "Phi", "Phi_tilde", "Group", "Score_center", "Score_Phi", "Score_Phi_tilde", "Score_internal", "Score_ratio")
#                          DiffNet = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
#                          return(DiffNet)
#                        })
#     end_time <- Sys.time()
#     time_diff = end_time - start_time
#     print(time_diff)
# 
#     # Select only strong and well classified edges
#     DiffNet_cl = DiffNet %>%
#       filter(Score_ratio > 1)
# 
#     # Categorize the Nodes into different roles
#     start_time <- Sys.time()
#     print("Calculating differentially co-expressed nodes")
#     Nodes = tryCatch(ClusterNodes(DiffNet = DiffNet_cl, cutoff.external = 0, cutoff.internal = 1),
#                      error=function(cond){
#                        print("Not enough nodes to cluster")
#                        # Empty table
#                        cols = c("Node", "Phi", "pval_Phi", "Group_pval_Phi", "Phi_tilde", "pval_Phi_Tilde", "Group_pval_Phi_tilde")
#                        Nodes = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
#                        return(Nodes)
#                       })
#     # Correct p-value for multiple testing
#     Nodes$pval_Phi_Tilde.adj.bonf = p.adjust(Nodes$pval_Phi_Tilde, method="bonferroni")
#     Nodes$pval_Phi_Tilde.adj.fdr = p.adjust(Nodes$pval_Phi_Tilde, method="fdr")
#     end_time <- Sys.time()
#     time_diff = end_time - start_time
#     print(time_diff)
#   } else {
#     # Empty tables
#     cols = c("Node.1", "Node.2", "D", "N", "Phi", "Phi_tilde", "Group", "Score_center", "Score_Phi", "Score_Phi_tilde", "Score_internal", "Score_ratio")
#     DiffNet_cl = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
#     cols = c("Node", "Phi", "pval_Phi", "Group_pval_Phi", "Phi_tilde", "pval_Phi_Tilde", "Group_pval_Phi_tilde", "pval_Phi_Tilde.adj.bonf", "pval_Phi_Tilde.adj.fdr", "type_correlation")
#     Nodes = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
#   }
#   # Save results
#   DiffNet_cl = DiffNet_cl %>% as.data.frame() %>% mutate(type_correlation = type_correlation_selected)
#   DiffNet_cl_all = rbind(DiffNet_cl_all, DiffNet_cl)
#   Nodes = Nodes %>% as.data.frame() %>% mutate(type_correlation = type_correlation_selected)
#   Nodes_all = rbind(Nodes_all, Nodes)
# 
#   rm(coexpression_filt_df_D)
#   rm(coexpression_filt_df_N)
#   rm(DiffNet_cl)
#   rm(Nodes)
# }
# 
# # Write final results
# DiffNet_cl_all %>% fwrite(output_edges_file)
# Nodes_all %>% fwrite(output_nodes_file)
