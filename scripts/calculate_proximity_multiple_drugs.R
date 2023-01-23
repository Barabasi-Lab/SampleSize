#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
source("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/Proximity_optimized.R")
library(optparse)
library(data.table)
library(dplyr)
library(igraph)
require(magrittr)

option_list = list(
  make_option(c("-a", "--case_genes_file"), action="store", type="character", 
              help="Input file (compete path) containing case genes", metavar="character"),
  make_option(c("-b", "--drug_targets_file"), action="store", type="character", 
              help="Input file (compete path) containing a dataframe with drugs (ID) and targets (Target)", metavar="character"),
  make_option(c("-c", "--ppi_file"), action="store", type="character", 
              help="Protein-protein interactions network file", metavar="character"),
  make_option(c("-d", "--ppi_distances_file"), action="store", type="character", default = NA,
              help="Protein-protein interactions distances file", metavar="character"),
  make_option(c("-e", "--output_file"), action="store", type="character", 
              help="Output proximity results file (complete path)", metavar="character"),
  make_option(c("-f", "--n_iterations"), action="store", type="integer", default=1000,
              help="Number of iterations", metavar="integer")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_proximity_multiple_drugs.R


# Define inputs
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

case_genes_file = opt$case_genes_file
drug_targets_file = opt$drug_targets_file
ppi_file = opt$ppi_file
ppi_distances_file = opt$ppi_distances_file
output_file = opt$output_file
n_iterations = as.numeric(opt$n_iterations)

# Example inputs
# case_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/nodes_by_type_and_category_and_size/codina_ppi_lcc_nodes_type_all_cat_disease-specific_size_20.txt"
# drug_targets_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/drug_targets/drugbank_targets_breast.neoplasm.txt" 
# ppi_file = "/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv"
# ppi_distances_file = "/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022_distances_matrix.txt"
# output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/proximity_by_type_and_category_and_size/codina_ppi_lcc_proximity_type_all_cat_disease-specific_size_20.txt"
# n_iterations = 1000

# Read PPI
ppi_df = fread(ppi_file) %>% dplyr::select(HGNC_Symbol.1, HGNC_Symbol.2) %>% dplyr::rename("Node.1"="HGNC_Symbol.1", "Node.2"="HGNC_Symbol.2")
ppi_net = graph_from_data_frame(ppi_df, directed=F) %>% simplify()
if ((is.na(ppi_distances_file)) | (is.null(ppi_distances_file))){
  ppi_distances = distances(ppi_net) %>% as.data.frame()
} else {
  ppi_distances = fread(ppi_distances_file) %>% as.data.frame()
}

# Read drug targets and filter by targets in the PPI network
drug_targets_df = fread(drug_targets_file) %>% 
  filter(Target %in% V(ppi_net)$name) %>% 
  dplyr::select(ID, Target) %>%
  unique()

# Read case genes and filter by genes in the PPI network
case_genes = unique(fread(case_genes_file)[,1][[1]])
case_genes = case_genes[case_genes %in% V(ppi_net)$name]

# Calculate proximity
codina_proximity_df = proximity_sign(N = n_iterations, 
                                     disease_genes = case_genes, # Disease genes
                                     Drugs = drug_targets_df, # Dataframe of drugs (ID) and targets (Target)
                                     distance = ppi_distances, # PPI network distances matrix
                                     disease_name = "disease",
                                     possible_targets = unique(drug_targets_df$Target)) # Targets
codina_proximity_df = codina_proximity_df %>% dplyr::select(-Disease)

codina_proximity_df %>% fwrite(output_file)
