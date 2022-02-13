#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(wTO)
set.seed(1510)

#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-a", "--coexpression_network_file"), type="character", 
              help="Co-expression network file", metavar="character"),
  make_option(c("-b", "--output_dir"), type="character", default=".", 
              help="Output directory [default= %default]", metavar="character"),
  make_option(c("-c", "--pvalue_threshold"), type="double", default = 0.05,
              help="P-value threshold", metavar="double"),
  make_option(c("-d", "--disparity_pvalue_threshold"), type="character", default=NA, 
              help="Disparity p-value [default= %default]", metavar="character"),
  make_option(c("-e", "--ppi_file"), type="character", default = "/home/j.aguirreplans/data/PPI/interactome_2019_merged_symbols.csv",
              help="File containing the protein-protein interaction network", metavar="character"),
  make_option(c("-f", "--disease_genes_file"), type="character", default = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022.csv",
              help="File containing the disease genes in the expression dataset", metavar="character"),
  make_option(c("-g", "--essential_genes_file"), type="character", default = "/home/j.aguirreplans/Databases/OGEE/OGEE_esential_genes.csv",
              help="File containing the essential genes in the expression dataset", metavar="character")
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
output_dir = opt$output_dir
pvalue_threshold = as.double(opt$pvalue_threshold)
disparity_pvalue_threshold = opt$disparity_pvalue_threshold
ppi_file = opt$ppi_file
disease_genes_file = opt$disease_genes_file
essential_genes_file = opt$essential_genes_file

#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/all_samples/wgcna_scipher_all_samples_size_100_rep_1.net"
#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/all_samples/spearman_scipher_all_samples_size_100_rep_1.net"
#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Whole.Blood/aracne_RNAseq_samples_Whole.Blood_size_90_rep_9.net"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_scipher/all_samples/analysis_by_top_scoring_edges_wgcna_scipher_all_samples_size_100_rep_1.txt"
#ppi_file = "/home/j.aguirreplans/data/PPI/interactome_2019_merged_symbols.csv"
#disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022.csv"
#essential_genes_file = "/home/j.aguirreplans/Databases/OGEE/OGEE_esential_genes.csv"
#pvalue_threshold = 0.05
#disparity_pvalue_threshold = NA

# Find method name in file name
file_split = strsplit(coexpression_network_file, split="/")[[1]]
file_name = file_split[length(file_split)]
file_split = strsplit(file_name, split="_")[[1]]
method = file_split[1]

read_coexpression_network <- function(coexpression_network_file, pvalue_threshold=0.05, disparity_pvalue_threshold=NA, RsquaredCut=0.80){
  cut = pvalue_threshold
  coexpression_network_df = as.data.frame(fread(coexpression_network_file))
  if(grepl('wto', coexpression_network_file, fixed = TRUE)){
    genes_network = unique(c(as.vector(coexpression_network_df$Node.1), as.vector(coexpression_network_df$Node.2)))
    coexpression_network_df = coexpression_network_df %>% rename("score"="wTO")
    coexpression_network_df = coexpression_network_df %>% filter(pval.adj < cut) %>% dplyr::select(Node.1, Node.2, score)
  } else if((grepl('spearman', coexpression_network_file, fixed = TRUE))){
    genes_network = unique(c(as.vector(coexpression_network_df$Node.1), as.vector(coexpression_network_df$Node.2)))
    coexpression_network_df = coexpression_network_df %>% rename("score"="spearman")
    if(!(is.na(disparity_pvalue_threshold))){
      if(as.numeric(pvalue_threshold) == 0.05){
        coexpression_network_df = coexpression_network_df %>% filter(disp.pvalue.after.filt.0.05 < as.numeric(disparity_pvalue_threshold)) %>% dplyr::select(Node.1, Node.2, score)
      } else if(as.numeric(pvalue_threshold) == 0.01){
        coexpression_network_df = coexpression_network_df %>% filter(disp.pvalue.after.filt.0.01 < as.numeric(disparity_pvalue_threshold)) %>% dplyr::select(Node.1, Node.2, score)
      } else {
        coexpression_network_df = coexpression_network_df %>% filter(disp.pvalue < cut) %>% dplyr::select(Node.1, Node.2, score)
      }
    } else {
      coexpression_network_df = coexpression_network_df %>% filter(p.adj < cut) %>% dplyr::select(Node.1, Node.2, score)
    }
  } else if((grepl('pearson', coexpression_network_file, fixed = TRUE))){
    genes_network = unique(c(as.vector(coexpression_network_df$Node.1), as.vector(coexpression_network_df$Node.2)))
    coexpression_network_df = coexpression_network_df %>% rename("score"="pearson")
    if(!(is.na(disparity_pvalue_threshold))){
      if(as.numeric(pvalue_threshold) == 0.05){
        coexpression_network_df = coexpression_network_df %>% filter(disp.pvalue.after.filt.0.05 < as.numeric(disparity_pvalue_threshold)) %>% dplyr::select(Node.1, Node.2, score)
      } else if(as.numeric(pvalue_threshold) == 0.01){
        coexpression_network_df = coexpression_network_df %>% filter(disp.pvalue.after.filt.0.01 < as.numeric(disparity_pvalue_threshold)) %>% dplyr::select(Node.1, Node.2, score)
      } else {
        coexpression_network_df = coexpression_network_df %>% filter(disp.pvalue < cut) %>% dplyr::select(Node.1, Node.2, score)
      }
    } else {
      coexpression_network_df = coexpression_network_df %>% filter(p.adj < cut) %>% dplyr::select(Node.1, Node.2, score)
    }
  } else if((grepl('wgcna', coexpression_network_file, fixed = TRUE))){
    library(WGCNA)
    genes_network = unique(names(coexpression_network_df))
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
    genes_network = unique(names(coexpression_network_df))
    rownames(coexpression_network_df) = colnames(coexpression_network_df)
    coexpression_network_df <- coexpression_network_df[order(rownames(coexpression_network_df)), order(colnames(coexpression_network_df))]
    coexpression_network_df = coexpression_network_df %>% wTO.in.line() %>% rename(score=wTO)
    coexpression_network_df = coexpression_network_df %>% filter(abs(score) > cut)
  } else {
    stop("Unknown co-expression network method!")
  }
  return(list("network" = coexpression_network_df, "genes_network" = genes_network, "cut" = cut))
}


coexpression_df = read_coexpression_network(coexpression_network_file=coexpression_network_file, pvalue_threshold=pvalue_threshold, disparity_pvalue_threshold=disparity_pvalue_threshold, RsquaredCut=0.80)
cut = coexpression_df$cut
genes_network = coexpression_df$genes_network
coexpression_df = coexpression_df$network
ppi_df = fread(ppi_file) %>% dplyr::select(proteinA_symbol, proteinB_symbol)
ppi_net = graph_from_data_frame(ppi_df, directed=F) %>% simplify()
disease_genes_df = fread(disease_genes_file)
disease_genes_df = disease_genes_df %>% filter(gene %in% genes_network)
essential_genes = unique(fread(essential_genes_file)$gene)
essential_genes = essential_genes[essential_genes %in% genes_network]

# Get PPI main core
k_core = coreness(ppi_net)
max_k = max(k_core)
ppi_main_core = induced.subgraph(ppi_net, vids = names(k_core[k_core==max_k]))

# Select PPIs in co-expression network
coexpression_ppi_df = inner_join(coexpression_df, ppi_df, by=c("Node.1"="proteinA_symbol", "Node.2"="proteinB_symbol"))
rm(coexpression_df)
coexpression_ppi_main_core_df = coexpression_ppi_df %>% filter((Node.1 %in% V(ppi_main_core)$name) & (Node.2 %in% V(ppi_main_core)$name))
coexpression_ppi_net = graph_from_data_frame(coexpression_ppi_df, directed=F) %>% simplify()
coexpression_ppi_main_core_net = graph_from_data_frame(coexpression_ppi_main_core_df, directed=F) %>% simplify()

#output_subgraph = "/home/j.aguirreplans/coexpression_ppi_test.csv"
coexpression_ppi_df %>% fwrite(output_subgraph)

# Write result dataframes
ppi_results_df = data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("cut", "disparity", "num_ppi_nodes", "num_ppi_edges", "faction_ppi_nodes", "fraction_ppi_edges", "num_ppi_main_core_nodes", "num_ppi_main_core_edges", "fraction_ppi_main_core_nodes", "fraction_ppi_main_core_edges"))))
ppi_results_df = rbind(ppi_results_df, data.frame(cut=cut, disparity=disparity_pvalue_threshold, num_ppi_nodes=gorder(coexpression_ppi_net), num_ppi_edges=gsize(coexpression_ppi_net), faction_ppi_nodes=gorder(coexpression_ppi_net)/gorder(ppi_net), fraction_ppi_edges=gsize(coexpression_ppi_net)/gsize(ppi_net), num_ppi_main_core_nodes=gorder(coexpression_ppi_main_core_net), num_ppi_main_core_edges=gsize(coexpression_ppi_main_core_net), fraction_ppi_main_core_nodes=gorder(coexpression_ppi_main_core_net)/gorder(ppi_main_core), fraction_ppi_main_core_edges=gsize(coexpression_ppi_main_core_net)/gsize(ppi_main_core)))

ppi_disease_gene_results_df <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c('cut', 'disparity', 'disease', 'disease_class', 'num_disease_genes', 'fraction_disease_genes', 'num_disease_gene_edges', 'fraction_disease_gene_edges'))




