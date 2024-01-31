#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
require(magrittr)
library(wTO)
set.seed(1510)



option_list = list(
  make_option(c("-l", "--list_coexpression_networks_file"), action="store", type="character", 
              help="File with the paths to the co-expression network files to merge", metavar="character"),
  make_option(c("-n", "--consensus_network_file"), action="store", type="character", 
              help="Output consensus co-expression network file", metavar="character"),
  make_option(c("-t", "--threshold"), action="store", type="double", default = 0.05,
              help="P-value threshold", metavar="double"),
  make_option(c("-m", "--method"), action="store", type="character", default = "pearson",
              help="Co-expression method", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_consensus_gene_coexpression_network.R -l /home/j.aguirreplans/Projects/Scipher/SampleSize/inputs/network_replicates_pearson_gtex_Breast.Mammary.Tissue_size_100.txt -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/reads/Breast.Mammary.Tissue/consensus/pearson_RNAseq_samples_Breast.Mammary.Tissue_size_100_consensus.net -m pearson -t 0.05



#----------------#
# READ ARGUMENTS #
#----------------#

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$list_coexpression_networks_file) | is.null(opt$consensus_network_file)){
  print_help(opt_parser)
  stop("Argument must be supplied (list_coexpression_networks_file, consensus_network_file).n", call.=FALSE)
}

list_coexpression_networks_file = opt$list_coexpression_networks_file
consensus_network_file = opt$consensus_network_file
threshold = as.double(opt$threshold)
method = opt$method

#list_coexpression_networks_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/inputs/network_replicates_pearson_tcga_TCGA-ESCA_size_100.txt"
#consensus_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/reads/tumor/TCGA-ESCA/consensus/pearson_RNAseq_samples_TCGA-ESCA_size_100_consensus.net"
#method = "pearson"
#threshold = 0.05


#-------------------------#
# DEFINITION OF FUNCTIONS #
#-------------------------#

###############################
## read_coexpression_network ##
###############################
# Read co-expression network, keeping filtering by p-value adjusted and keeping p-value
#'  @param coexpression_network_file The path to the co-expression network file
#'  @param method The method used to create the co-expression network. 
#'  @param threshold P-value threshold to filter the co-expression network. 

read_coexpression_network_for_consensus <- function(coexpression_network_file, method, threshold, set_unsignificant_scores_to_zero=TRUE){
  
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
      coexpression_df = coexpression_df %>% dplyr::select(Node.1, Node.2, score, pval, pval.adj)
      coexpression_df$score = ifelse(coexpression_df$pval.adj >= threshold, 0, coexpression_df$score)
      coexpression_df = coexpression_df %>% dplyr::select(-pval.adj)
    } else {
      # If set_unsignificant_scores_to_zero is FALSE, filter out unsignificant edges
      coexpression_df = coexpression_df %>% dplyr::select(Node.1, Node.2, score, pval, pval.adj) %>% filter(pval.adj < threshold) %>% dplyr::select(-pval.adj)
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

###################
## wTO.Consensus ##
###################
#'  @param data 

wTO.Consensus = function(data){
  if (!is.list(data)){
    stop("data must be a list of data.frames.")
  }
  ### Weight
  #weight = pval = nodes = list()
  weight = pval = data.frame()
  
  for ( i in 1:length(data)){
    cols = colnames(data[[i]])
    if(nrow(weight)==0){
      weight=data[[i]][,1:3] %>% rename(!!(paste(cols[3], i, sep="")) := cols[3])
      pval=data[[i]][,c(1:2,4)] %>% rename(!!(paste(cols[4], i, sep="")) := cols[4])
    }else{
      weight=weight %>% inner_join((data[[i]][,1:3] %>% rename(!!(paste(cols[3], i, sep="")) := cols[3])), by=c("Node.1", "Node.2"))
      pval=pval %>% inner_join((data[[i]][,c(1:2,4)] %>% rename(!!(paste(cols[4], i, sep="")) := cols[4])), by=c("Node.1", "Node.2"))
    }
    #weight[[i]] = data[[i]][,1:3]
    #pval[[i]] = data[[i]][,c(1:2,4)]
    #ID = unique(c(as.character(data[[i]]$Node.1), as.character(data[[i]]$Node.2)))
    #nodes[[i]] = data.frame(ID =  ID)
    #names(weight[[i]])[3] = paste0(names(weight[[i]])[3], i)
    #names(pval[[i]])[3] = paste0(names(pval[[i]])[3], i)
  }
  
  #weight = plyr::join_all(weight, type = 'full', by=c("Node.1", "Node.2"))
  #pval = plyr::join_all(pval, type = 'full', by=c("Node.1", "Node.2"))
  #nodes = plyr::join_all(nodes, type = 'inner', by=c("ID"))
  
  #message(paste('Total common nodes:', nrow(nodes)))
  #weight = subset(weight, weight$Node.1 %in% nodes$ID & weight$Node.2 %in% nodes$ID)
  #pval = subset(pval, pval$Node.1 %in% nodes$ID & pval$Node.2 %in% nodes$ID)
  pval[is.na(pval)] <- 1
  weight[is.na(weight)] <- 0.01
  
  wTOCN = CN_aux(weight[, -c(1:2)])
  pvalue_fisher = fishermethod(pval[, -c(1:2)])
  
  Out = data.frame(Node.1 = pval[,1], Node.2 = pval[,2],
                   CN = wTOCN, pval.fisher = pvalue_fisher)
  return(Out)
}

##################
## fishermethod ##
##################
#'  @param data_x 

fishermethod = function(data_x){
  chi = rowSums(log(data_x))*-2
  pval = sapply(chi, function(x) stats::pchisq(x, 2*ncol(data_x), lower.tail = FALSE))
  return(pval)
}

############
## CN_aux ##
############
#'  @param data_x 

CN_aux = function(data_x){
  abs_x = apply(data_x, 2, abs)
  sum_abs_x = apply(abs_x, 1, sum)
  div = (abs_x/sum_abs_x) * data_x
  wTO_cons = apply(div, 1, sum)
  return(wTO_cons)
}


#---------------#
# READ NETWORKS #
#---------------#

# Read list of networks
list_coexpression_networks = unique((fread(list_coexpression_networks_file, header=F))$V1)

# # Read networks
# networks_to_merge = list()
# for (x in 1:length(list_coexpression_networks)){
#   network_file = list_coexpression_networks[x]
#   networks_to_merge[[x]] = read_coexpression_network_for_consensus(coexpression_network_file=list_coexpression_networks[x], method=method, threshold=threshold, set_unsignificant_scores_to_zero=F)
#   #networks_to_merge[[x]] = (fread(list_coexpression_networks[x]))[,1:4]
# }

# Read networks and create consensus network in chunks
chunk_size = 5000000
connections = list()
for(x in 1:length(list_coexpression_networks)){
  network_file = list_coexpression_networks[x]
  connections[[x]] = file(network_file, "r")
}

consensus_network = data.frame()
x=0
chunk_count = 0
repeat {
  subnetworks_to_merge = list()
  x=x+1
  chunk_count=chunk_count+chunk_size
  message(paste('Chunk:', chunk_count))
  skip_param=0
  if(x==1){
    skip_param=1
  }
  for(x in 1:length(connections)){
    subnetworks_to_merge[[x]] <- (read.table(connections[[x]], nrows=chunk_size, skip=skip_param, header=FALSE, fill = TRUE, sep=",", col.names=c("Node.1", "Node.2", "score", "pval", "pval.adj"), colClasses=c("character", "character", "double", "double", "double"))) %>% dplyr::select(Node.1, Node.2, score, pval)
    #subnetworks_to_merge[[x]] = (readr::read_delim(connections[[x]], delim = ",", n_max=chunk_size, skip=skip_param, col_names=c("Node.1", "Node.2", "score", "pval", "pval.adj"), col_types=c("c", "c", "d", "d", "d"))) %>% dplyr::select(Node.1, Node.2, score, pval)
  }
  consensus_subnetwork = wTO.Consensus(data = subnetworks_to_merge)
  if(nrow(consensus_network) == 0){
    consensus_network = data.frame(consensus_subnetwork)
  } else {
    consensus_network = rbind(consensus_network, consensus_subnetwork)
  }
  if(nrow(subnetworks_to_merge[[1]]) < chunk_size){
    for(x in 1:length(connections)){
      close(connections[[x]])
    }
    break
  }
  rm(subnetworks_to_merge)
  rm(consensus_subnetwork)
}


#--------------------------#
# CREATE CONSENSUS NETWORK #
#--------------------------#

#consensus_network = wTO.Consensus(data = networks_to_merge)
consensus_network %>% fwrite(consensus_network_file)

