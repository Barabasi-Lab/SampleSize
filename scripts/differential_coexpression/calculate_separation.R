#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
library(igraph)
require(magrittr)
library(NetSci)

option_list = list(
  make_option(c("-a", "--input_file"), action="store", type="character", 
              help="Input dataframe (compete path) containing Category and Node", metavar="character"),
  make_option(c("-b", "--ppi_file"), action="store", type="character", 
              help="Protein-protein interactions network file", metavar="character"),
  make_option(c("-c", "--output_file"), action="store", type="character", 
              help="Output proximity results file (complete path)", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_separation.R 


# Define inputs
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_file = opt$input_file
ppi_file = opt$ppi_file
output_file = opt$output_file

# Read PPI
ppi_df = fread(ppi_file) %>% dplyr::select(HGNC_Symbol.1, HGNC_Symbol.2) %>% dplyr::rename("Node.1"="HGNC_Symbol.1", "Node.2"="HGNC_Symbol.2")
ppi_net = graph_from_data_frame(ppi_df, directed=F) %>% simplify()

# Read input data and filter by genes in the PPI network
All = fread(input_file) %>% dplyr::rename(Category=1, Node=2) %>% dplyr::filter(Node %in% V(ppi_net)$name)

# Calculate Jaccard
Jac = Jaccard(All)

# Calculate separation
S = separation_Significance(G = ppi_net,
                            ST = All,
                            correct_by_target = FALSE,
                            Threads = 1)

# Merge results
S_result = Jac %>% full_join(S, by=c("Node.1"="x", "Node.2"="y")) %>% dplyr::rename("Category.1"="Node.1", "Category.2"="Node.2")

# Output file
S_result %>% fwrite(output_file)
