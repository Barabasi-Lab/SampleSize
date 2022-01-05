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
  make_option(c("-c", "--coexpression_network_file"), type="character", 
              help="Co-expression network file", metavar="character"),
  make_option(c("-p", "--ppi_file"), type="character", 
              help="Protein-protein interaction network file.", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="comparison_coexpression_ppi.txt", 
              help="Output file [default= %default]", metavar="character")
  #make_option(c("-o", "--output_dir"), type="character", default="comparison_coexpression_to_ppi", 
  #            help="Output dir [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_topology_coexpression_network.R -c /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -p /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_all_samples.net/comparison_coexpression_ppi.txt

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if ((is.null(opt$coexpression_network_file)) | (is.null(opt$ppi_file))){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coexpression_network_file, ppi_file).n", call.=FALSE)
}
tissue_sex_network_file = opt$coexpression_network_file
ppi_tissue_sex_filt_file = opt$ppi_file
#output_dir = opt$output_dir
output_file = opt$output_file

#tissue_sex_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/aracne_RNAseq_samples_Spleen_female_all_samples.net"
#ppi_tissue_sex_filt_file = "/home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv"
#output_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/aracne_RNAseq_samples_Spleen_female_all_samples.net"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/comparison_coexpression_ppi_aracne_RNAseq_samples_Spleen_female_all_samples.net"

# Read PPI network
ppi_tissue_sex_filt_df = fread(ppi_tissue_sex_filt_file) %>% select(proteinA_symbol, proteinB_symbol)
ppi_tissue_sex_filt_net = graph_from_data_frame(ppi_tissue_sex_filt_df, directed=F) %>% simplify()

# Read tissue-sex network for all samples
tissue_sex_network_df = as.data.frame(fread(tissue_sex_network_file))
if (!(grepl('wto', tissue_sex_network_file, fixed = TRUE))){
  # Pass the matrix to the wTO in-line format
  rownames(tissue_sex_network_df) = colnames(tissue_sex_network_df)
  tissue_sex_network_df <- tissue_sex_network_df[order(rownames(tissue_sex_network_df)), order(colnames(tissue_sex_network_df))]
  tissue_sex_network_df = tissue_sex_network_df %>% wTO.in.line() %>% rename(score=wTO)
} else {
  tissue_sex_network_df = tissue_sex_network_df %>% rename(score=wTO)
}

# Calculate topological metrics of co-expression network
topology_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("nodes", "edges", "av_degree", "av_path_length", "av_clustering_coef", "num_components", "size_lcc"))))



topology_df %>% fwrite(output_file)


