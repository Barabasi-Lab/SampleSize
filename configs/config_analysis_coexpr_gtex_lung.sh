#!/bin/bash

input_folder="/scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung"
output_results_dir="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_gtex/reads/Lung"
output_subgraphs_dir="/scratch/j.aguirreplans/Scipher/SampleSizeProject/subgraphs_gtex/reads/Lung"
genes_dataset_file="/work/ccnr/j.aguirreplans/Databases/GTEx/v8/reads/genes_from_rnaseq_filtered_files_by_tissue/gtex_rnaseq_gene_info_Lung.csv"
ppi_file="/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv"
disease_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_gtex.csv"
essential_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/essential_genes/OGEE_essential_genes_gtex.csv"
threshold=0.05
method="spearman"
