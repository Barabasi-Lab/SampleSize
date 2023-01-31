#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
source("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/GO.R")
library(optparse)
library(data.table)
library(dplyr)
require(magrittr)
set.seed(1510)

option_list = list(
  make_option(c("-a", "--case_genes_file"), action="store", type="character", 
              help="Input file (compete path) containing case genes", metavar="character"),
  make_option(c("-b", "--background_genes_file"), action="store", type="character", 
              help="Input file (compete path) containing background genes", metavar="character"),
  make_option(c("-c", "--output_file"), action="store", type="character", 
              help="Output enrichment results file (complete path)", metavar="character"),
  make_option(c("-d", "--ontology"), action="store", type="character", default="BP",
              help="Ontology used for the GO enrichment (BP, MF, CC)", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_go_enrichment.R 


# Define inputs
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

case_genes_file = opt$case_genes_file
background_genes_file = opt$background_genes_file
output_file = opt$output_file
ontology_selected = opt$ontology

#case_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/nodes_by_type_and_category_and_size/codina_ppi_lcc_nodes_type_all_cat_common_size_20.txt"
#background_genes_file = "/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022_nodes.csv"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/functions_by_type_and_category_and_size/codina_ppi_lcc_go_enrichment_type_all_cat_common_size_20.txt"
#ontology_selected = "BP"

# Read input files and get unique values of first column
case_genes = unique(fread(case_genes_file)[,1]) %>% pull()
background_genes = unique(fread(background_genes_file)[,1]) %>% pull()

# Calculate GO enrichment
GO_enrichment = GO(ID_type = "symbol", 
                   g = case_genes,
                   ONTO = ontology_selected, 
                   bg = background_genes)

# Write output
GO_enrichment$Res %>% fwrite(output_file)
