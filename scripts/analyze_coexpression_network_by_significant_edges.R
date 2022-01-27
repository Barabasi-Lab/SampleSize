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
  make_option(c("-a", "--coexpression_network_file"), type="character", 
              help="Co-expression network file", metavar="character"),
  make_option(c("-b", "--pvalue_threshold"), type="double", default = 0.05,
              help="P-value threshold", metavar="double"),
  make_option(c("-c", "--disease_genes_file"), type="character", default = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info.csv",
              help="File containing the disease genes in the expression dataset", metavar="character"),
  make_option(c("-d", "--essential_genes_file"), type="character", default = "/home/j.aguirreplans/Databases/OGEE/OGEE_esential_genes.csv",
              help="File containing the essential genes in the expression dataset", metavar="character"),
  make_option(c("-e", "--output_subgraph"), type="character", default="network.txt", 
              help="Output network filtered [default= %default]", metavar="character"),
  make_option(c("-f", "--output_main_core_subgraph"), type="character", default="main_core.txt", 
              help="Output main core subgraph [default= %default]", metavar="character"),
  make_option(c("-g", "--output_disease_genes_subgraph_dir"), type="character", default="disease_subgraphs", 
              help="Output directory to store disease module subraphs [default= %default]", metavar="character"),
  make_option(c("-i", "--output_essential_genes_subgraph"), type="character", default="main_core.txt", 
              help="Output essential genes subgraph [default= %default]", metavar="character"),
  make_option(c("-j", "--output_topology_file"), type="character", default="analysis_topology.txt", 
              help="Output topology file [default= %default]", metavar="character"),
  make_option(c("-k", "--output_disease_genes_file"), type="character", default="analysis_essential_genes.txt", 
              help="Output disease genes file [default= %default]", metavar="character"),
  make_option(c("-l", "--output_essential_genes_file"), type="character", default="analysis_disease_genes.txt", 
              help="Output essential genes file [default= %default]", metavar="character")
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
pvalue_threshold = as.double(opt$pvalue_threshold)
disease_genes_file = opt$disease_genes_file
essential_genes_file = opt$essential_genes_file
output_subgraph = opt$output_subgraph
output_main_core_subgraph = opt$output_main_core_subgraph
output_essential_genes_subgraph = opt$output_essential_genes_subgraph
output_disease_genes_subgraph_dir = opt$output_disease_genes_subgraph_dir
output_topology_file = opt$output_topology_file
output_disease_genes_file = opt$output_disease_genes_file
output_essential_genes_file = opt$output_essential_genes_file

#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/all_samples/wgcna_scipher_all_samples_size_100_rep_1.net"
#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/all_samples/spearman_scipher_all_samples_size_100_rep_1.net"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_scipher/all_samples/analysis_by_top_scoring_edges_wgcna_scipher_all_samples_size_100_rep_1.txt"
#disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info.csv"
#essential_genes_file = "/home/j.aguirreplans/Databases/OGEE/OGEE_esential_genes.csv"

read_coexpression_network <- function(coexpression_network_file, pvalue_threshold=0.05, RsquaredCut=0.80){
  cut = pvalue_threshold
  coexpression_network_df = as.data.frame(fread(coexpression_network_file))
  if(grepl('wto', coexpression_network_file, fixed = TRUE)){
    coexpression_network_df = coexpression_network_df %>% rename("score"="wTO")
    coexpression_network_df = coexpression_network_df %>% filter(pval.adj < cut) %>% dplyr::select(Node.1, Node.2, score)
  } else if((grepl('spearman', coexpression_network_file, fixed = TRUE))){
    coexpression_network_df = coexpression_network_df %>% rename("score"="spearman")
    coexpression_network_df = coexpression_network_df %>% filter(pvalue < cut) %>% dplyr::select(Node.1, Node.2, score)
  } else if((grepl('pearson', coexpression_network_file, fixed = TRUE))){
    coexpression_network_df = coexpression_network_df %>% rename("score"="pearson")
    coexpression_network_df = coexpression_network_df %>% filter(pvalue < cut) %>% dplyr::select(Node.1, Node.2, score)
  } else if((grepl('wgcna', coexpression_network_file, fixed = TRUE))){
    library(WGCNA)
    rownames(coexpression_network_df) = colnames(coexpression_network_df)
    coexpression_network_df <- coexpression_network_df[order(rownames(coexpression_network_df)), order(colnames(coexpression_network_df))]
    hard_threshold = WGCNA::pickHardThreshold(coexpression_network_df, RsquaredCut = RsquaredCut)
    if (is.na(hard_threshold$cutEstimate)){
      cut = hard_threshold$fitIndices[hard_threshold$fitIndices$SFT.R.sq == max(hard_threshold$fitIndices$SFT.R.sq),]$Cut
    } else {
      cut = hard_threshold$cutEstimate
    }
    coexpression_network_df = coexpression_network_df %>% wTO.in.line() %>% rename(score=wTO)
    coexpression_network_df = coexpression_network_df %>% filter(abs(score) > cut)
  } else if((grepl('aracne', coexpression_network_file, fixed = TRUE))){
    cut = 0
    rownames(coexpression_network_df) = colnames(coexpression_network_df)
    coexpression_network_df <- coexpression_network_df[order(rownames(coexpression_network_df)), order(colnames(coexpression_network_df))]
    coexpression_network_df = coexpression_network_df %>% wTO.in.line() %>% rename(score=wTO)
    coexpression_network_df = coexpression_network_df %>% filter(abs(score) > cut)
  } else {
    stop("Unknown co-expression network method!")
  }
  return(list("network" = coexpression_network_df, "cut" = cut))
}


coexpression_df = read_coexpression_network(coexpression_network_file=coexpression_network_file, pvalue_threshold=pvalue_threshold, RsquaredCut=0.80)
cut = coexpression_df$cut
coexpression_df = coexpression_df$network
genes_network = unique(c(as.vector(coexpression_df$Node.1), as.vector(coexpression_df$Node.2)))
disease_genes_df = fread(disease_genes_file)
disease_genes_df = disease_genes_df %>% filter(gene %in% genes_network)
essential_genes = unique(fread(essential_genes_file)$gene)
essential_genes = essential_genes[essential_genes %in% genes_network]


# Define result dataframes
topology_df = data.frame(matrix(ncol=14,nrow=0, dimnames=list(NULL, c("cut", "num_nodes", "num_edges", "av_degree", "av_path_length", "av_clustering_coef", "num_components", "num_lcc_nodes", "num_lcc_edges", "lcc_z", "lcc_pvalue", "max_k", "num_main_core_nodes", "num_main_core_edges"))))
essential_gene_results_df = data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("cut", "num_essential_genes", "fraction_essential_genes", "num_components", "num_lcc_nodes", "fraction_essential_lcc_nodes", "num_lcc_edges", "lcc_z", "lcc_pvalue"))))
disease_gene_results_df <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c('cut', 'disease', 'disease_class', 'num_disease_genes', 'fraction_disease_genes', 'num_disease_components', 'num_disease_lcc_nodes', 'fraction_disease_lcc_nodes', 'num_disease_lcc_edges', "disease_lcc_z", "disease_lcc_pvalue"))

# Get network in igraph format
coexpression_net = graph_from_data_frame((coexpression_df %>% rename("weight"="score")), directed=F) %>% simplify()

# Output network
coexpression_df %>% fwrite(output_subgraph)

# Calculate topological parameters
num_nodes = gorder(coexpression_net)
num_edges = gsize(coexpression_net)
av_degree = mean(degree(coexpression_net))
av_path_length = mean_distance(coexpression_net)
av_clustering_coef = transitivity(coexpression_net)

# Components analysis
components_net = components(coexpression_net)
lcc = induced.subgraph(coexpression_net, vids = V(coexpression_net)[components_net$membership == which.max(components_net$csize)] )
num_lcc_nodes = gorder(lcc)
num_lcc_edges = gsize(lcc)
lcc_sig = LCC_Significance(N = 1000, 
                           Targets = V(lcc)$name, 
                           G = coexpression_net)

# K-core analysis
k_core = coreness(coexpression_net)
max_k = max(k_core)
main_core = induced.subgraph(coexpression_net, vids = names(k_core[k_core==max_k]))
num_main_core_nodes = gorder(main_core)
num_main_core_edges = gsize(main_core)

topology_df = rbind(topology_df, data.frame(cut=cut, num_nodes=num_nodes, num_edges=num_edges, av_degree=av_degree, av_path_length=av_path_length, av_clustering_coef=av_clustering_coef, num_components=components_net$no, num_lcc_nodes=num_lcc_nodes, num_lcc_edges=num_lcc_edges, lcc_z=lcc_sig$Z, lcc_pvalue=lcc_sig$emp_p, max_k=max_k, num_main_core_nodes=num_main_core_nodes, num_main_core_edges=num_main_core_edges))

# Output main core subgraph
write_graph(main_core, output_main_core_subgraph, "graphml")

# Essential genes analysis
essential_genes_in_network = essential_genes[essential_genes %in% V(coexpression_net)$name]
essential_subgraph = induced.subgraph(coexpression_net, vids=essential_genes_in_network)
essential_components = components(essential_subgraph)
essential_lcc = induced.subgraph(essential_subgraph, vids = V(essential_subgraph)[essential_components$membership == which.max(essential_components$csize)] )
num_essential_lcc_nodes = gorder(essential_lcc)
num_essential_lcc_edges = gsize(essential_lcc)
essential_lcc_sig = LCC_Significance(N = 1000, Targets = V(essential_lcc)$name, G = coexpression_net)
essential_gene_results_df = rbind(essential_gene_results_df, data.frame(cut=cut, num_essential_genes=length(essential_genes_in_network), fraction_essential_genes=length(essential_genes_in_network)/length(essential_genes), num_components=essential_components$no, num_essential_lcc_nodes=num_essential_lcc_nodes, fraction_essential_lcc_nodes=num_essential_lcc_nodes/length(essential_genes), num_lcc_edges=num_essential_lcc_edges, lcc_z=essential_lcc_sig$Z, lcc_pvalue=essential_lcc_sig$emp_p))

# Output essential genes subgraph
write_graph(essential_subgraph, output_essential_genes_subgraph, "graphml")

# Disease genes analysis
dir.create(output_disease_genes_subgraph_dir, showWarnings = FALSE)
for (disease in sort(unique(disease_genes_df$disease))){
  disease_selected_df = disease_genes_df %>% filter(disease == !!disease)
  disease_selected_genes = unique(disease_selected_df$gene)
  disease_selected_genes_in_network = disease_selected_genes[disease_selected_genes %in% V(coexpression_net)$name]
  if (length(disease_selected_genes) > 19){
    disease_subgraph = induced.subgraph(coexpression_net, vids=disease_selected_genes_in_network)
    disease_components = components(disease_subgraph)
    disease_lcc = induced.subgraph(disease_subgraph, vids = V(disease_subgraph)[disease_components$membership == which.max(disease_components$csize)] )
    num_disease_lcc_nodes = gorder(disease_lcc)
    num_disease_lcc_edges = gsize(disease_lcc)
    disease_lcc_sig = LCC_Significance(N = 1000, Targets = V(disease_lcc)$name, G = coexpression_net)
    disease_selected_classes = unique(disease_selected_df$root.concept)
    
    disease_net_df = coexpression_df %>% filter((Node.1 %in% V(disease_subgraph)$name) & (Node.2 %in% V(disease_subgraph)$name))
    disease_name_no_sp_char = unique((disease_genes_df %>% filter(disease == !!disease))$disease.no.sp.char)
    output_disease_genes_network_file = paste(output_disease_genes_subgraph_dir, paste(disease_name_no_sp_char, ".txt", sep=""), sep="/")
    write_graph(disease_subgraph, output_disease_genes_network_file, "graphml")
    
    for (disease_selected_class in disease_selected_classes){
      disease_gene_results_df = rbind(disease_gene_results_df, data.frame(cut=cut, disease=disease, disease_class=disease_selected_class, num_disease_genes=length(disease_selected_genes_in_network), fraction_disease_genes=length(disease_selected_genes_in_network)/length(disease_selected_genes), num_disease_components=disease_components$no, num_disease_lcc_nodes=num_disease_lcc_nodes, fraction_disease_lcc_nodes=num_disease_lcc_nodes/length(disease_selected_genes), num_disease_lcc_edges=num_disease_lcc_edges, disease_lcc_z=disease_lcc_sig$Z, disease_lcc_pvalue=disease_lcc_sig$emp_p))
    }
  }
}


# Create output file
topology_df %>% fwrite(output_topology_file)
essential_gene_results_df %>% fwrite(output_essential_genes_file)
disease_gene_results_df %>% fwrite(output_disease_genes_file)


