#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(wTO)
set.seed(1510)
library(NetSci)

option_list = list(
  make_option(c("-c", "--coexpression_network_file"), type="character", 
              help="Co-expression network file", metavar="character"),
  make_option(c("-o", "--output_results_dir"), type="character", default="./results", 
              help="Output results directory [default= %default]", metavar="character"),
  make_option(c("-s", "--output_subgraphs_dir"), type="character", default="./subgraphs", 
              help="Output networks directory [default= %default]", metavar="character"),
  make_option(c("-f", "--file_name"), type="character", default="./subgraphs", 
              help="File name [default= %default]", metavar="character"),
  make_option(c("-t", "--threshold"), type="double", default = 0.05,
              help="P-value threshold", metavar="double"),
  make_option(c("-p", "--ppi_file"), type="character", default = "/home/j.aguirreplans/data/PPI/interactome_2019_merged_symbols.csv",
              help="File containing the protein-protein interaction network", metavar="character"),
  make_option(c("-d", "--disease_genes_file"), type="character", default = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_scipher.csv",
              help="File containing the disease genes in the expression dataset", metavar="character"),
  make_option(c("-e", "--essential_genes_file"), type="character", default = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/essential_genes/OGEE_essential_genes_scipher.csv",
              help="File containing the essential genes in the expression dataset", metavar="character"),
  make_option(c("-g", "--genes_dataset_file"), type="character", default = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/Dec2021/00_data/scipher_rnaseq_gene_info.csv",
              help="File containing the genes in the expression dataset", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_network_by_significant_edges.R -c /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv 


#----------------#
# READ ARGUMENTS #
#----------------#

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$coexpression_network_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coexpression_network_file, ppi_file).n", call.=FALSE)
}

coexpression_network_file = opt$coexpression_network_file
output_results_dir = opt$output_results_dir
output_subgraphs_dir = opt$output_subgraphs_dir
file_name = opt$file_name
threshold = as.double(opt$threshold)
ppi_file = opt$ppi_file
disease_genes_file = opt$disease_genes_file
essential_genes_file = opt$essential_genes_file
genes_dataset_file = opt$genes_dataset_file

#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA/pearson_RNAseq_samples_TCGA_size_100_rep_1.net"
#output_results_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_scipher/all_samples"
#output_subgraphs_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/subgraphs_scipher/all_samples"
#file_name = "spearman_scipher_all_samples_size_100_rep_1_pvalue_0.05"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_scipher/all_samples/analysis_by_top_scoring_edges_wgcna_scipher_all_samples_size_100_rep_1.txt"
#ppi_file = "/home/j.aguirreplans/data/PPI/interactome_2019_merged_symbols.csv"
#disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_scipher.csv"
#disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_tcga.csv"
#disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_gtex.csv"
#essential_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/essential_genes/OGEE_essential_genes_scipher.csv"
#genes_dataset_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/Dec2021/00_data/scipher_rnaseq_gene_info.csv"
#genes_dataset_file = "/home/j.aguirreplans/Databases/TCGA/2022-03-28-Dataset/TCGA/out/tcga_rnaseq_gene_info.csv"


#-------------------------#
# DEFINITION OF FUNCTIONS #
#-------------------------#

###############################
## read_coexpression_network ##
###############################
#'  @param coexpression_network_file The path to the co-expression network file
#'  @param method The method used to create the co-expression network. 
#'  @param threshold P-value threshold to filter the co-expression network. 

read_coexpression_network <- function(coexpression_network_file, method, threshold){
  
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
    coexpression_df = coexpression_df %>% filter(abs(score) > threshold)
  } else {
    if(method == "wto"){
      coexpression_df = coexpression_df %>% rename("score"= "wTO")
    }
    coexpression_df = coexpression_df %>% dplyr::select(Node.1, Node.2, score, pval.adj) %>% filter(pval.adj < threshold)
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


###################################
## calculate_topology_parameters ##
###################################
#'  @param network_graph Co-expression network in igraph format
#'  @param output_main_core_subgraph_file Path to the output file of the main core subgraph. If NA, no file is written.

calculate_topology_parameters <- function(network_graph, output_main_core_subgraph_file=NA){
  
  start_time_topology <- Sys.time()
  print("Calculating topology")
  
  # Calculate number of nodes and edges
  num_nodes = gorder(network_graph)
  num_edges = gsize(network_graph)
  
  # Calculate average degree
  start_time <- Sys.time()
  print("Calculating average degree")
  av_degree = mean(degree(network_graph))
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  # Calculate average path length (it might take long time)    
  start_time <- Sys.time()
  print("Calculating average path length")
  #av_path_length = mean_distance(network_graph)
  av_path_length = NA
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  # Calculate average clustering coefficient
  start_time <- Sys.time()
  print("Calculating average clustering coefficient")
  #av_clustering_coef = transitivity(network_graph)
  av_clustering_coef = NA
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  # Calculate number of components, nodes and edges in the LCC, LCC significance (it might take long time)
  start_time <- Sys.time()
  print("Calculating LCC")
  components_net = components(network_graph)
  lcc = induced.subgraph(network_graph, vids = V(network_graph)[components_net$membership == which.max(components_net$csize)] )
  num_lcc_nodes = gorder(lcc)
  num_lcc_edges = gsize(lcc)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  start_time <- Sys.time()
  print("Calculating LCC significance")
  # LCC significance is the calculation that takes more time!!! Several hours depending on the size of the network
  # With bin=1, the computation time is lower but the degree distribution is not maintained and the hypothesis changes
  #if (length(V(lcc)$name) > 0){
  #  lcc_sig = LCC_Significance(N=1000, Targets=V(lcc)$name, G=network_graph,  bins=1)
  #} else {
  #  lcc_sig = list(Z=NA, emp_p=NA) # If empty network, leave result as NA
  #}
  lcc_sig = list(Z=NA, emp_p=NA)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  # K-core analysis
  start_time <- Sys.time()
  print("Calculating K core")
  k_core = coreness(network_graph)
  max_k = max(k_core)
  main_core = induced.subgraph(network_graph, vids = names(k_core[k_core==max_k]))
  num_main_core_nodes = gorder(main_core)
  num_main_core_edges = gsize(main_core)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  # Output main core subgraph
  if(!(is.na(output_main_core_subgraph_file))){
    if(!(file.exists(output_main_core_subgraph_file))){
      write_graph(main_core, output_main_core_subgraph_file, "graphml")
    }
  }
  
  # Create results table
  topology_results_df = data.frame(num_nodes=num_nodes, num_edges=num_edges, av_degree=av_degree, av_path_length=av_path_length, av_clustering_coef=av_clustering_coef, num_components=components_net$no, num_lcc_nodes=num_lcc_nodes, num_lcc_edges=num_lcc_edges, lcc_z=lcc_sig$Z, lcc_pvalue=lcc_sig$emp_p, max_k=max_k, num_main_core_nodes=num_main_core_nodes, num_main_core_edges=num_main_core_edges)
  
  end_time_topology <- Sys.time()
  time_diff = end_time_topology - start_time_topology
  print(time_diff)
  
  return(topology_results_df)
}

##############################
## calculate_ppi_parameters ##
##############################
#'  @param network_df Co-expression network in dataframe format
#'  @param ppi_df Protein-protein interactions network in dataframe format
#'  @param output_ppi_subgraph_file Path to the output file of the co-expression - PPI subgraph. If NA, no file is written.
#'  @param output_ppi_main_core_subgraph_file Path to the output file of the co-expression - PPI main core subgraph. If NA, no file is written.

calculate_ppi_parameters <- function(network_df, ppi_df, output_ppi_subgraph_file=NA, output_ppi_main_core_subgraph_file=NA){
  
  start_time_ppi <- Sys.time()
  print("Calculating PPI")
  
  # Select PPIs in co-expression network
  coexpression_ppi_df = inner_join(network_df, ppi_df, by=c("Node.1"="proteinA_symbol", "Node.2"="proteinB_symbol"))
  coexpression_ppi_net = graph_from_data_frame(coexpression_ppi_df, directed=F) %>% simplify()
  # Select main core PPIs in co-expression network
  coexpression_ppi_main_core_df = coexpression_ppi_df %>% filter((Node.1 %in% V(ppi_main_core)$name) & (Node.2 %in% V(ppi_main_core)$name))
  coexpression_ppi_main_core_net = graph_from_data_frame(coexpression_ppi_main_core_df, directed=F) %>% simplify()
  
  # Write result networks
  if(!(is.na(output_ppi_subgraph_file))){
    if(!(file.exists(output_ppi_subgraph_file))){
      coexpression_ppi_df %>% fwrite(output_ppi_subgraph_file)
    }
  }
  if(!(is.na(output_ppi_main_core_subgraph_file))){
    if(!(file.exists(output_ppi_main_core_subgraph_file))){
      coexpression_ppi_main_core_df %>% fwrite(output_ppi_main_core_subgraph_file)
    }
  }
  
  # Create results table
  ppi_results_df = data.frame(num_ppi_nodes=gorder(coexpression_ppi_net), num_ppi_edges=gsize(coexpression_ppi_net), fraction_ppi_nodes=gorder(coexpression_ppi_net)/gorder(ppi_net), fraction_ppi_edges=gsize(coexpression_ppi_net)/gsize(ppi_net), num_ppi_main_core_nodes=gorder(coexpression_ppi_main_core_net), num_ppi_main_core_edges=gsize(coexpression_ppi_main_core_net), fraction_ppi_main_core_nodes=gorder(coexpression_ppi_main_core_net)/gorder(ppi_main_core), fraction_ppi_main_core_edges=gsize(coexpression_ppi_main_core_net)/gsize(ppi_main_core))
  
  end_time_ppi <- Sys.time()
  time_diff = end_time_ppi - start_time_ppi
  print(time_diff)
  
  return(ppi_results_df)
}

#######################################
## calculate_disease_gene_parameters ##
#######################################
#'  @param network_df Co-expression network in dataframe format
#'  @param network_graph Co-expression network in igraph format
#'  @param disease_genes_df Disease-gene associations dataframe
#'  @param output_disease_genes_subgraph_dir Path to the output directory of the co-expression - disease genes subgraphs. If NA, no file is written.

calculate_disease_gene_parameters <- function(network_df, network_graph, disease_genes_df, output_disease_genes_subgraph_dir=NA){
  
  start_time_disease <- Sys.time()
  print("Calculating disease genes")
  
  # Create directory to store subgraphs (if necessary)
  if(!(is.na(output_disease_genes_subgraph_dir))){
    dir.create(output_disease_genes_subgraph_dir, showWarnings = FALSE)
  }
  
  # Calculate analysis for a list of selected diseases
  cols = c("disease", "disease_class", "num_disease_genes", "num_disease_edges", "fraction_disease_genes", "num_disease_components", "num_disease_lcc_nodes", "num_disease_lcc_edges", "fraction_disease_lcc_nodes", "disease_lcc_z", "disease_lcc_pvalue")
  disease_gene_results_df <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
  selected_diseases = c("alzheimer disease", "arthritis rheumatoid", "cardiomyopathies", "diabetes mellitus type 2")
  #for (disease in sort(unique(disease_genes_df$DiseaseName))){
  for (disease in selected_diseases){
    # Select genes associated with the disease
    disease_genes_selected_df = disease_genes_df %>% filter(DiseaseName == !!disease)
    disease_genes_selected = unique(disease_genes_selected_df$HGNC_Symbol)
    disease_genes_selected_in_network = disease_genes_selected[disease_genes_selected %in% V(network_graph)$name]
    # Get disease subgraph
    disease_subgraph = induced.subgraph(network_graph, vids=disease_genes_selected_in_network)
    disease_net_df = network_df %>% filter((Node.1 %in% V(disease_subgraph)$name) & (Node.2 %in% V(disease_subgraph)$name))
    # Get calculations of the disease subgraph (number of genes, components, LCC, LCC significance)
    num_disease_genes_subgraph = length(unique(c(disease_net_df$Node.1, disease_net_df$Node.2)))
    disease_components = igraph::components(disease_subgraph)
    disease_lcc = igraph::induced.subgraph(disease_subgraph, vids = V(disease_subgraph)[disease_components$membership == which.max(disease_components$csize)] )
    num_disease_lcc_nodes = gorder(disease_lcc)
    num_disease_lcc_edges = gsize(disease_lcc)
    if (length(V(disease_lcc)$name) > 0){
      disease_lcc_sig = LCC_Significance(N = 1000, Targets = V(disease_lcc)$name, G = network_graph, bins=1) # With bin=1, the degree distribution is not maintained and the hypothesis changes
    } else {
      disease_lcc_sig = list(Z=NA, emp_p=NA) # If empty network, leave result as NA
    }
    #disease_lcc_sig = list(Z=NA, emp_p=NA)
    # Write result network
    if(!(is.na(output_disease_genes_subgraph_dir))){
      disease_name_no_sp_char = unique((disease_genes_df %>% filter(DiseaseName == !!disease))$DiseaseName.no.sp.char)
      output_disease_genes_network_file = paste(output_disease_genes_subgraph_dir, paste(disease_name_no_sp_char, ".txt", sep=""), sep="/")
      if(!(file.exists(output_disease_genes_network_file))){
        write_graph(disease_subgraph, output_disease_genes_network_file, "graphml")
      }
    }
    # Create results table
    disease_selected_classes = strsplit(unique(disease_genes_selected_df$DescriptorName), split="|", fixed=TRUE)[[1]]
    for (disease_selected_class in disease_selected_classes){
      disease_gene_results_df = rbind(disease_gene_results_df, data.frame(disease=disease, disease_class=disease_selected_class, 
                                                                          num_disease_genes=num_disease_genes_subgraph, num_disease_edges=gsize(disease_subgraph), 
                                                                          fraction_disease_genes=num_disease_genes_subgraph/length(disease_genes_selected), 
                                                                          num_disease_components=disease_components$no, num_disease_lcc_nodes=num_disease_lcc_nodes, num_disease_lcc_edges=num_disease_lcc_edges, 
                                                                          fraction_disease_lcc_nodes=num_disease_lcc_nodes/length(disease_genes_selected),
                                                                          disease_lcc_z=disease_lcc_sig$Z, disease_lcc_pvalue=disease_lcc_sig$emp_p))
    }
  }
  
  end_time_disease <- Sys.time()
  time_diff = end_time_disease - start_time_disease
  print(time_diff)
  
  return(disease_gene_results_df)
}

#########################################
## calculate_essential_gene_parameters ##
#########################################
#'  @param network_df Co-expression network in dataframe format
#'  @param network_graph Co-expression network in igraph format
#'  @param essential_genes Dataframe of essential genes
#'  @param output_essential_genes_subgraph Path to the output file of the co-expression - essential genes subgraph. If NA, no file is written.

calculate_essential_gene_parameters <- function(network_df, network_graph, essential_genes, output_essential_genes_subgraph=NA){
  
  start_time_essential <- Sys.time()
  print("Calculating essential genes")
  
  # Get essential genes in network
  essential_genes_in_network = essential_genes[essential_genes %in% V(network_graph)$name]
  # Get essential genes subgraph
  essential_subgraph = induced.subgraph(network_graph, vids=essential_genes_in_network)
  essential_subgraph_df = network_df %>% filter((Node.1 %in% V(essential_subgraph)$name) & (Node.2 %in% V(essential_subgraph)$name))
  num_essential_genes_subgraph = length(unique(c(essential_subgraph_df$Node.1, essential_subgraph_df$Node.2)))
  rm(essential_subgraph_df)
  # Calculate components formed by essential genes, LCC, LCC significance
  essential_components = components(essential_subgraph)
  essential_lcc = induced.subgraph(essential_subgraph, vids = V(essential_subgraph)[essential_components$membership == which.max(essential_components$csize)] )
  num_essential_lcc_nodes = gorder(essential_lcc)
  num_essential_lcc_edges = gsize(essential_lcc)
  #if (length(V(essential_lcc)$name) > 0){
  #  essential_lcc_sig = LCC_Significance(N=1000, Targets=V(essential_lcc)$name, G=network_graph, bins=1)
  #} else {
  #  essential_lcc_sig = list(Z=NA, emp_p=NA) # If empty network, leave result as NA
  #}
  essential_lcc_sig = list(Z=NA, emp_p=NA)
  
  # Create results table
  essential_gene_results_df = data.frame(num_essential_genes=num_essential_genes_subgraph, num_essential_edges=gsize(essential_subgraph), fraction_essential_genes=num_essential_genes_subgraph/length(essential_genes), num_components=essential_components$no, num_lcc_nodes=num_essential_lcc_nodes, num_lcc_edges=num_essential_lcc_edges, fraction_essential_lcc_nodes=num_essential_lcc_nodes/length(essential_genes), lcc_z=essential_lcc_sig$Z, lcc_pvalue=essential_lcc_sig$emp_p)
  
  # Output essential genes subgraph
  if(!(is.na(output_essential_genes_subgraph))){
    if(!(file.exists(output_essential_genes_subgraph))){
      write_graph(essential_subgraph, output_essential_genes_subgraph, "graphml")
    }
  }
  
  end_time_essential <- Sys.time()
  time_diff = end_time_essential - start_time_essential
  print(time_diff)
  
  return(essential_gene_results_df)
}



#--------------#
# READ NETWORK #
#--------------#

# Find method name in file name
file_split = strsplit(file_name, split="_")[[1]]
method = file_split[1]

start_time <- Sys.time()
print("Reading co-expression network")
coexpression_df = read_coexpression_network(coexpression_network_file=coexpression_network_file, method=method, threshold=threshold)
genes_dataset = (fread(genes_dataset_file) %>% filter(enough_counts == TRUE))$HGNC_Symbol
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# Output network filtered by p-value
output_subgraph = paste(output_subgraphs_dir, "/subgraphs/", file_name, "_subgraph.txt", sep="")
if(!(file.exists(output_subgraph))){
  coexpression_df %>% fwrite(output_subgraph)
}


#---------------------#
# READ OTHER DATASETS #
#---------------------#

# Read disease genes file and filter by diseases with at least 20 genes
disease_genes_df = fread(disease_genes_file) %>% filter(HGNC_Symbol %in% genes_dataset) %>% unique() %>% 
  group_by(DiseaseName) %>%
  mutate(TotalDiseaseGenesDataset = n()) %>%
  filter(TotalDiseaseGenesDataset > 19) %>%
  ungroup()

# Read essential genes file
essential_genes = unique((fread(essential_genes_file) %>% filter(HGNC_Symbol %in% genes_dataset))$HGNC_Symbol)

# Read PPI
ppi_df = fread(ppi_file) %>% dplyr::select(proteinA_symbol, proteinB_symbol) %>% filter((proteinA_symbol %in% genes_dataset) & (proteinB_symbol %in% genes_dataset))
ppi_net = graph_from_data_frame(ppi_df, directed=F) %>% simplify()
# Calculate k-core
k_core = coreness(ppi_net)
max_k = max(k_core)
ppi_main_core = induced.subgraph(ppi_net, vids = names(k_core[k_core==max_k]))


#-------------------------------#
# CALCULATE TOPOLOGY PARAMETERS #
#-------------------------------#

output_topology_file = paste(output_results_dir, "/", file_name, "_analysis_topology.txt", sep="")
if(!(file.exists(output_topology_file))){
  
  cols = c("type_correlation", "threshold", "num_nodes", "num_edges", "av_degree", "av_path_length", "av_clustering_coef", "num_components", "num_lcc_nodes", "num_lcc_edges", "lcc_z", "lcc_pvalue", "max_k", "num_main_core_nodes", "num_main_core_edges")
  topology_results_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
  
  # Calculate results for different levels of correlation
  type_correlation_df = data.frame(row.names=c("all", "very weak", "weak", "moderate", "strong", "very strong", "strong-very strong", "moderate-strong-very strong", "weak-moderate-strong-very strong"), lower_val=c(NA,0,0.2,0.4,0.6,0.8,0.6,0.4,0.2), upper_val=c(NA,0.2,0.4,0.6,0.8,NA,NA,NA,NA))
  for(type_correlation in row.names(type_correlation_df)){
    
    print(paste("Correlation type: ", type_correlation, sep=""))
    
    # Filter network by correlation
    if(type_correlation == "all"){
      coexpression_filt_df = data.frame(coexpression_df) # same data frame without filtering
      output_main_core_subgraph_file = paste(output_subgraphs_dir, "/main_core/", file_name, "_main_core_subgraph.txt", sep="")
    } else {
      coexpression_filt_df = filter_network_by_correlation(network_df=coexpression_df, upper_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$upper_val, lower_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$lower_val)
      output_main_core_subgraph_file = NA # we do not plot the subgraph for filtered networks
    }
    
    # Get filtered network in igraph format
    coexpression_filt_net = graph_from_data_frame((coexpression_filt_df %>% rename("weight"="score")), directed=F) %>% simplify()
    
    # Calculate parameters
    topology_correlation_results_df = calculate_topology_parameters(network_graph=coexpression_filt_net, output_main_core_subgraph_file=output_main_core_subgraph_file)
    topology_results_df = rbind(topology_results_df, cbind(data.frame(type_correlation=rep(type_correlation, nrow(topology_correlation_results_df)), threshold=rep(threshold, nrow(topology_correlation_results_df))), topology_correlation_results_df))
    rm(coexpression_filt_df)
    rm(coexpression_filt_net)
    
  }
  
  topology_results_df %>% fwrite(output_topology_file)
  
}


#--------------------------#
# CALCULATE PPI PARAMETERS #
#--------------------------#

output_ppi_file = paste(output_results_dir, "/", file_name, "_analysis_ppi.txt", sep="")
if(!(file.exists(output_ppi_file))){
  
  cols = c("type_correlation", "threshold", "num_ppi_nodes", "num_ppi_edges", "fraction_ppi_nodes", "fraction_ppi_edges", "num_ppi_main_core_nodes", "num_ppi_main_core_edges", "fraction_ppi_main_core_nodes", "fraction_ppi_main_core_edges")
  ppi_results_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))

  # Calculate results for different levels of correlation
  type_correlation_df = data.frame(row.names=c("all", "very weak", "weak", "moderate", "strong", "very strong", "strong-very strong", "moderate-strong-very strong", "weak-moderate-strong-very strong"), lower_val=c(NA,0,0.2,0.4,0.6,0.8,0.6,0.4,0.2), upper_val=c(NA,0.2,0.4,0.6,0.8,NA,NA,NA,NA))
  for(type_correlation in row.names(type_correlation_df)){
    
    print(paste("Correlation type: ", type_correlation, sep=""))
    
    # Filter network by correlation
    if(type_correlation == "all"){
      coexpression_filt_df = data.frame(coexpression_df) # same data frame without filtering
      output_ppi_subgraph_file = paste(output_subgraphs_dir, "/ppi/", file_name, "_ppi_subgraph.txt", sep="")
      output_ppi_main_core_subgraph_file = paste(output_subgraphs_dir, "/ppi_main_core/", file_name, "_ppi_main_core_subgraph.txt", sep="")
    } else {
      coexpression_filt_df = filter_network_by_correlation(network_df=coexpression_df, upper_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$upper_val, lower_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$lower_val)
      output_ppi_subgraph_file = NA # we do not plot the subgraphs for filtered networks
      output_ppi_main_core_subgraph_file = NA
    }
    
    # Calculate parameters
    ppi_correlation_results_df = calculate_ppi_parameters(network_df=coexpression_filt_df, ppi_df=ppi_df, output_ppi_subgraph_file=output_ppi_subgraph_file, output_ppi_main_core_subgraph_file=output_ppi_main_core_subgraph_file)
    ppi_results_df = rbind(ppi_results_df, cbind(data.frame(type_correlation=rep(type_correlation, nrow(ppi_correlation_results_df)), threshold=rep(threshold, nrow(ppi_correlation_results_df))), ppi_correlation_results_df))
    rm(coexpression_filt_df)

  }
  
  ppi_results_df %>% fwrite(output_ppi_file)
  
}


#------------------------------------#
# CALCULATE DISEASE GENES PARAMETERS #
#------------------------------------#

output_disease_genes_file = paste(output_results_dir, "/", file_name, "_analysis_disease_genes.txt", sep="")
if(!(file.exists(output_disease_genes_file))){
  
  cols = c("type_correlation", "threshold", "disease", "disease_class", "num_disease_genes", "num_disease_edges", "fraction_disease_genes", "num_disease_components", "num_disease_lcc_nodes", "num_disease_lcc_edges", "fraction_disease_lcc_nodes", "disease_lcc_z", "disease_lcc_pvalue")
  disease_gene_results_df <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
  
  # Calculate results for different levels of correlation
  type_correlation_df = data.frame(row.names=c("all", "very weak", "weak", "moderate", "strong", "very strong", "strong-very strong", "moderate-strong-very strong", "weak-moderate-strong-very strong"), lower_val=c(NA,0,0.2,0.4,0.6,0.8,0.6,0.4,0.2), upper_val=c(NA,0.2,0.4,0.6,0.8,NA,NA,NA,NA))
  for(type_correlation in row.names(type_correlation_df)){
    
    print(paste("Correlation type: ", type_correlation, sep=""))
    
    # Filter network by correlation
    if(type_correlation == "all"){
      coexpression_filt_df = data.frame(coexpression_df) # same data frame without filtering
      output_disease_genes_subgraph_dir = paste(output_subgraphs_dir, "/disease_genes/", file_name, "_disease_genes_subgraphs", sep="")
    } else {
      coexpression_filt_df = filter_network_by_correlation(network_df=coexpression_df, upper_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$upper_val, lower_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$lower_val)
      output_disease_genes_subgraph_dir = NA # we do not plot subgraphs for filtered networks
    }
    
    # Get filtered network in igraph format
    coexpression_filt_net = graph_from_data_frame((coexpression_filt_df %>% rename("weight"="score")), directed=F) %>% simplify()
    
    # Calculate parameters
    disease_gene_correlation_results_df = calculate_disease_gene_parameters(network_df=coexpression_filt_df, network_graph=coexpression_filt_net, disease_genes_df=disease_genes_df, output_disease_genes_subgraph_dir=output_disease_genes_subgraph_dir)
    disease_gene_results_df = rbind(disease_gene_results_df, cbind(data.frame(type_correlation=rep(type_correlation, nrow(disease_gene_correlation_results_df)), threshold=rep(threshold, nrow(disease_gene_correlation_results_df))), disease_gene_correlation_results_df))
    rm(coexpression_filt_df)
    rm(coexpression_filt_net)
    
  }
  
  disease_gene_results_df %>% fwrite(output_disease_genes_file)
  
}


#--------------------------------------#
# CALCULATE ESSENTIAL GENES PARAMETERS #
#--------------------------------------#

output_essential_genes_file = paste(output_results_dir, "/", file_name, "_analysis_essential_genes.txt", sep="")
if(!(file.exists(output_essential_genes_file))){
  
  cols = c("type_correlation", "threshold", "num_essential_genes", "num_essential_edges", "fraction_essential_genes", "num_components", "num_lcc_nodes", "num_lcc_edges", "fraction_essential_lcc_nodes", "lcc_z", "lcc_pvalue")
  essential_gene_results_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
  
  # Calculate results for different levels of correlation
  type_correlation_df = data.frame(row.names=c("all", "very weak", "weak", "moderate", "strong", "very strong", "strong-very strong", "moderate-strong-very strong", "weak-moderate-strong-very strong"), lower_val=c(NA,0,0.2,0.4,0.6,0.8,0.6,0.4,0.2), upper_val=c(NA,0.2,0.4,0.6,0.8,NA,NA,NA,NA))
  for(type_correlation in row.names(type_correlation_df)){
    
    print(paste("Correlation type: ", type_correlation, sep=""))
    
    # Filter network by correlation
    if(type_correlation == "all"){
      coexpression_filt_df = data.frame(coexpression_df) # same data frame without filtering
      output_essential_genes_subgraph = paste(output_subgraphs_dir, "/essential_genes/", file_name, "_essential_genes_subgraph.txt", sep="")
    } else {
      coexpression_filt_df = filter_network_by_correlation(network_df=coexpression_df, upper_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$upper_val, lower_threshold=type_correlation_df[row.names(type_correlation_df) == type_correlation,]$lower_val)
      output_essential_genes_subgraph = NA # we do not plot subgraphs for filtered networks
    }
    
    # Get filtered network in igraph format
    coexpression_filt_net = graph_from_data_frame((coexpression_filt_df %>% rename("weight"="score")), directed=F) %>% simplify()
    
    # Calculate parameters
    essential_gene_correlation_results_df = calculate_essential_gene_parameters(network_df=coexpression_filt_df, network_graph=coexpression_filt_net, essential_genes=essential_genes, output_essential_genes_subgraph=output_essential_genes_subgraph)
    essential_gene_results_df = rbind(essential_gene_results_df, cbind(data.frame(type_correlation=rep(type_correlation, nrow(essential_gene_correlation_results_df)), threshold=rep(threshold, nrow(essential_gene_correlation_results_df))), essential_gene_correlation_results_df))
    rm(coexpression_filt_df)
    rm(coexpression_filt_net)
    
  }
  
  essential_gene_results_df %>% fwrite(output_essential_genes_file)
  
}


