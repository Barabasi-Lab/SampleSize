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

#### READ ARGUMENTS ####
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
  make_option(c("-i", "--disparity_threshold"), type="character", default=NA, 
              help="Disparity p-value [default= %default]", metavar="character"),
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
disparity_threshold = opt$disparity_threshold
ppi_file = opt$ppi_file
disease_genes_file = opt$disease_genes_file
essential_genes_file = opt$essential_genes_file
genes_dataset_file = opt$genes_dataset_file

#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/complete_dataset/aracne_scipher_complete_dataset_size_100_rep_1.net"
#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/all_samples/wgcna_scipher_all_samples_size_100_rep_1.net"
#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/all_samples/spearman_scipher_all_samples_size_100_rep_1.net"
#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Whole.Blood/aracne_scipher_all_samples_size_100_rep_1.net"
#output_results_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_scipher/all_samples"
#output_subgraphs_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/subgraphs_scipher/all_samples"
#file_name = "spearman_scipher_all_samples_size_100_rep_1_pvalue_0.05_disp_None"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_scipher/all_samples/analysis_by_top_scoring_edges_wgcna_scipher_all_samples_size_100_rep_1.txt"
#ppi_file = "/home/j.aguirreplans/data/PPI/interactome_2019_merged_symbols.csv"
#disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_scipher.csv"
#essential_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/essential_genes/OGEE_essential_genes_scipher.csv"
#genes_dataset_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/Dec2021/00_data/scipher_rnaseq_gene_info.csv"


# Define output result files
output_topology_file = paste(output_results_dir, "/", file_name, "_analysis_topology.txt", sep="")
output_ppi_file = paste(output_results_dir, "/", file_name, "_analysis_ppi.txt", sep="")
output_disease_genes_file = paste(output_results_dir, "/", file_name, "_analysis_disease_genes.txt", sep="")
output_essential_genes_file = paste(output_results_dir, "/", file_name, "_analysis_essential_genes.txt", sep="")

# Define output subgraphs
output_subgraph = paste(output_subgraphs_dir, "/subgraphs/", file_name, "_subgraph.txt", sep="")
output_main_core_subgraph = paste(output_subgraphs_dir, "/main_core/", file_name, "_main_core_subgraph.txt", sep="")
output_ppi_subgraph = paste(output_subgraphs_dir, "/ppi/", file_name, "_ppi_subgraph.txt", sep="")
output_ppi_main_core_subgraph = paste(output_subgraphs_dir, "/ppi_main_core/", file_name, "_ppi_main_core_subgraph.txt", sep="")
output_disease_genes_subgraph_dir = paste(output_subgraphs_dir, "/disease_genes/", file_name, "_disease_genes_subgraphs", sep="")
output_essential_genes_subgraph = paste(output_subgraphs_dir, "/essential_genes/", file_name, "_essential_genes_subgraph.txt", sep="")


# Find method name in file name
file_split = strsplit(coexpression_network_file, split="/")[[1]]
file_name = file_split[length(file_split)]
file_split = strsplit(file_name, split="_")[[1]]
method = file_split[1]


# Read gene co-expression network
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

start_time <- Sys.time()
print("Reading co-expression network")
coexpression_df = read_coexpression_network(coexpression_network_file=coexpression_network_file, method=method, threshold=threshold)
genes_dataset = (fread(genes_dataset_file) %>% filter(enough_counts == TRUE))$symbol
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)


# Read other datasets
# Select diseases with more than 19 genes
disease_genes_df = fread(disease_genes_file) %>% filter(hgnc_symbol %in% genes_dataset) %>% unique() %>% 
  group_by(DiseaseName) %>%
  mutate(TotalDiseaseGenesDataset = n()) %>%
  filter(TotalDiseaseGenesDataset > 19) %>%
  ungroup()
essential_genes = unique((fread(essential_genes_file) %>% filter(gene %in% genes_dataset))$gene)
ppi_df = fread(ppi_file) %>% dplyr::select(proteinA_symbol, proteinB_symbol) %>% filter((proteinA_symbol %in% genes_dataset) & (proteinB_symbol %in% genes_dataset))
ppi_net = graph_from_data_frame(ppi_df, directed=F) %>% simplify()
k_core = coreness(ppi_net)
max_k = max(k_core)
ppi_main_core = induced.subgraph(ppi_net, vids = names(k_core[k_core==max_k]))


# Get network in igraph format
coexpression_net = graph_from_data_frame((coexpression_df %>% rename("weight"="score")), directed=F) %>% simplify()


# Define result dataframes
topology_df = data.frame(matrix(ncol=15,nrow=0, dimnames=list(NULL, c("threshold", "disparity", "num_nodes", "num_edges", "av_degree", "av_path_length", "av_clustering_coef", "num_components", "num_lcc_nodes", "num_lcc_edges", "lcc_z", "lcc_pvalue", "max_k", "num_main_core_nodes", "num_main_core_edges"))))
ppi_results_df = data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("threshold", "disparity", "num_ppi_nodes", "num_ppi_edges", "fraction_ppi_nodes", "fraction_ppi_edges", "num_ppi_main_core_nodes", "num_ppi_main_core_edges", "fraction_ppi_main_core_nodes", "fraction_ppi_main_core_edges"))))
disease_gene_results_df <- setNames(data.frame(matrix(ncol = 13, nrow = 0)), c('threshold', 'disparity', 'disease', 'disease_class', 'num_disease_genes', 'num_disease_edges', 'fraction_disease_genes', 'num_disease_components', 'num_disease_lcc_nodes', 'num_disease_lcc_edges', 'fraction_disease_lcc_nodes', "disease_lcc_z", "disease_lcc_pvalue"))
essential_gene_results_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("threshold", "disparity", "num_essential_genes", "num_essential_edges", "fraction_essential_genes", "num_components", "num_lcc_nodes", "num_lcc_edges", "fraction_essential_lcc_nodes", "lcc_z", "lcc_pvalue"))))


# Calculate topology analysis
start_time_topology <- Sys.time()
print("Calculating topology")
if(!(file.exists(output_topology_file))){
  
  # Output network
  coexpression_df %>% fwrite(output_subgraph)
  
  # Calculate topological parameters
  num_nodes = gorder(coexpression_net)
  num_edges = gsize(coexpression_net)
  
  start_time <- Sys.time()
  print("Calculating average degree")
  av_degree = mean(degree(coexpression_net))
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  start_time <- Sys.time()
  print("Calculating average path length")
  av_path_length = mean_distance(coexpression_net)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  start_time <- Sys.time()
  print("Calculating average clustering coefficient")
  av_clustering_coef = transitivity(coexpression_net)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  # Components analysis
  start_time <- Sys.time()
  print("Calculating LCC")
  components_net = components(coexpression_net)
  lcc = induced.subgraph(coexpression_net, vids = V(coexpression_net)[components_net$membership == which.max(components_net$csize)] )
  num_lcc_nodes = gorder(lcc)
  num_lcc_edges = gsize(lcc)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  start_time <- Sys.time()
  print("Calculating LCC significance")
  # LCC significance is the calculation that takes more time!!! Several hours depending on the size of the network
  # So maybe I will skip it by now
  #lcc_sig = LCC_Significance(N = 1000, 
  #                           Targets = V(lcc)$name, 
  #                           G = coexpression_net)
  lcc_sig = list(Z=NA, emp_p=NA)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  # K-core analysis
  start_time <- Sys.time()
  print("Calculating K core")
  k_core = coreness(coexpression_net)
  max_k = max(k_core)
  main_core = induced.subgraph(coexpression_net, vids = names(k_core[k_core==max_k]))
  num_main_core_nodes = gorder(main_core)
  num_main_core_edges = gsize(main_core)
  end_time <- Sys.time()
  time_diff = end_time - start_time
  print(time_diff)
  
  topology_df = rbind(topology_df, data.frame(threshold=threshold, disparity=disparity_threshold, num_nodes=num_nodes, num_edges=num_edges, av_degree=av_degree, av_path_length=av_path_length, av_clustering_coef=av_clustering_coef, num_components=components_net$no, num_lcc_nodes=num_lcc_nodes, num_lcc_edges=num_lcc_edges, lcc_z=lcc_sig$Z, lcc_pvalue=lcc_sig$emp_p, max_k=max_k, num_main_core_nodes=num_main_core_nodes, num_main_core_edges=num_main_core_edges))
  
  # Output main core subgraph
  write_graph(main_core, output_main_core_subgraph, "graphml")
  
  # Create output topology file
  topology_df %>% fwrite(output_topology_file)
  
}
end_time_topology <- Sys.time()
time_diff = end_time_topology - start_time_topology
print(time_diff)


# Calculate PPI analysis
start_time <- Sys.time()
print("Calculating PPI")
if(!(file.exists(output_ppi_file))){
  
  # Select PPIs in co-expression network
  coexpression_ppi_df = inner_join(coexpression_df, ppi_df, by=c("Node.1"="proteinA_symbol", "Node.2"="proteinB_symbol"))
  coexpression_ppi_main_core_df = coexpression_ppi_df %>% filter((Node.1 %in% V(ppi_main_core)$name) & (Node.2 %in% V(ppi_main_core)$name))
  coexpression_ppi_net = graph_from_data_frame(coexpression_ppi_df, directed=F) %>% simplify()
  coexpression_ppi_main_core_net = graph_from_data_frame(coexpression_ppi_main_core_df, directed=F) %>% simplify()
  
  # Write result networks
  coexpression_ppi_df %>% fwrite(output_ppi_subgraph)
  coexpression_ppi_main_core_df %>% fwrite(output_ppi_main_core_subgraph)

  # Write result dataframes
  ppi_results_df = rbind(ppi_results_df, data.frame(threshold=threshold, disparity=disparity_threshold, num_ppi_nodes=gorder(coexpression_ppi_net), num_ppi_edges=gsize(coexpression_ppi_net), fraction_ppi_nodes=gorder(coexpression_ppi_net)/gorder(ppi_net), fraction_ppi_edges=gsize(coexpression_ppi_net)/gsize(ppi_net), num_ppi_main_core_nodes=gorder(coexpression_ppi_main_core_net), num_ppi_main_core_edges=gsize(coexpression_ppi_main_core_net), fraction_ppi_main_core_nodes=gorder(coexpression_ppi_main_core_net)/gorder(ppi_main_core), fraction_ppi_main_core_edges=gsize(coexpression_ppi_main_core_net)/gsize(ppi_main_core)))
  ppi_results_df %>% fwrite(output_ppi_file)
}
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)


# Calculate disease gene analysis
start_time <- Sys.time()
print("Calculating disease genes")
if(!(file.exists(output_disease_genes_file))){
  
  # Disease genes analysis
  dir.create(output_disease_genes_subgraph_dir, showWarnings = FALSE)
  selected_diseases = c("alzheimer disease", "arthritis rheumatoid", "cardiomyopathies", "diabetes mellitus type 2")
  #for (disease in sort(unique(disease_genes_df$DiseaseName))){
  for (disease in selected_diseases){
    disease_selected_df = disease_genes_df %>% filter(DiseaseName == !!disease)
    disease_selected_genes = unique(disease_selected_df$hgnc_symbol)
    disease_selected_genes_in_network = disease_selected_genes[disease_selected_genes %in% V(coexpression_net)$name]
    disease_subgraph = induced.subgraph(coexpression_net, vids=disease_selected_genes_in_network)
    disease_net_df = coexpression_df %>% filter((Node.1 %in% V(disease_subgraph)$name) & (Node.2 %in% V(disease_subgraph)$name))
    num_disease_genes_subgraph = length(unique(c(disease_net_df$Node.1, disease_net_df$Node.2)))
    disease_components = igraph::components(disease_subgraph)
    disease_lcc = igraph::induced.subgraph(disease_subgraph, vids = V(disease_subgraph)[disease_components$membership == which.max(disease_components$csize)] )
    num_disease_lcc_nodes = gorder(disease_lcc)
    num_disease_lcc_edges = gsize(disease_lcc)
    #disease_lcc_sig = LCC_Significance(N = 1000, Targets = V(disease_lcc)$name, G = coexpression_net)
    disease_lcc_sig = list(Z=NA, emp_p=NA)
    disease_selected_classes = strsplit(unique(disease_selected_df$DescriptorName), split="|", fixed=TRUE)[[1]]
    
    disease_name_no_sp_char = unique((disease_genes_df %>% filter(DiseaseName == !!disease))$DiseaseName.no.sp.char)
    output_disease_genes_network_file = paste(output_disease_genes_subgraph_dir, paste(disease_name_no_sp_char, ".txt", sep=""), sep="/")
    write_graph(disease_subgraph, output_disease_genes_network_file, "graphml")
    
    for (disease_selected_class in disease_selected_classes){
      disease_gene_results_df = rbind(disease_gene_results_df, data.frame(threshold=threshold, disparity=disparity_threshold, disease=disease, disease_class=disease_selected_class, 
                                                                          num_disease_genes=num_disease_genes_subgraph, num_disease_edges=gsize(disease_subgraph), 
                                                                          fraction_disease_genes=num_disease_genes_subgraph/length(disease_selected_genes), 
                                                                          num_disease_components=disease_components$no, num_disease_lcc_nodes=num_disease_lcc_nodes, num_disease_lcc_edges=num_disease_lcc_edges, 
                                                                          fraction_disease_lcc_nodes=num_disease_lcc_nodes/length(disease_selected_genes),
                                                                          disease_lcc_z=disease_lcc_sig$Z, disease_lcc_pvalue=disease_lcc_sig$emp_p))
    }
  }
  
  # Create output disease genes file
  disease_gene_results_df %>% fwrite(output_disease_genes_file)
  
}
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)


start_time <- Sys.time()
print("Calculating essential genes")
if(!(file.exists(output_essential_genes_file))){
  
  # Essential genes analysis
  essential_genes_in_network = essential_genes[essential_genes %in% V(coexpression_net)$name]
  essential_subgraph = induced.subgraph(coexpression_net, vids=essential_genes_in_network)
  essential_subgraph_df = coexpression_df %>% filter((Node.1 %in% V(essential_subgraph)$name) & (Node.2 %in% V(essential_subgraph)$name))
  num_essential_genes_subgraph = length(unique(c(essential_subgraph_df$Node.1, essential_subgraph_df$Node.2)))
  rm(essential_subgraph_df)
  essential_components = components(essential_subgraph)
  essential_lcc = induced.subgraph(essential_subgraph, vids = V(essential_subgraph)[essential_components$membership == which.max(essential_components$csize)] )
  num_essential_lcc_nodes = gorder(essential_lcc)
  num_essential_lcc_edges = gsize(essential_lcc)
  #essential_lcc_sig = LCC_Significance(N = 1000, Targets = V(essential_lcc)$name, G = coexpression_net)
  essential_lcc_sig = list(Z=NA, emp_p=NA)
  essential_gene_results_df = rbind(essential_gene_results_df, data.frame(threshold=threshold, disparity=disparity_threshold, num_essential_genes=num_essential_genes_subgraph, num_essential_edges=gsize(essential_subgraph), fraction_essential_genes=num_essential_genes_subgraph/length(essential_genes), num_components=essential_components$no, num_lcc_nodes=num_essential_lcc_nodes, num_lcc_edges=num_essential_lcc_edges, fraction_essential_lcc_nodes=num_essential_lcc_nodes/length(essential_genes), lcc_z=essential_lcc_sig$Z, lcc_pvalue=essential_lcc_sig$emp_p))

  # Output essential genes subgraph
  write_graph(essential_subgraph, output_essential_genes_subgraph, "graphml")
  
  # Create output essential genes file
  essential_gene_results_df %>% fwrite(output_essential_genes_file)
  
}
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

