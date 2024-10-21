#!/bin/bash

input_folder="/scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tissue/Lung"
output_results_dir="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_tcga/reads/tissue/Lung"
output_subgraphs_dir="/scratch/j.aguirreplans/Scipher/SampleSizeProject/subgraphs_tcga/reads/tissue/Lung"
genes_dataset_file="/work/ccnr/j.aguirreplans/Databases/TCGA/2022-11-18-Dataset/TCGA/out/reads/tissue/filter_genes_low_counts/genes_filtered_files_by_sample_group/rnaseq_gene_info_Lung.csv"
ppi_file="/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv"
disease_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_tcga.csv"
essential_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/essential_genes/OGEE_essential_genes_tcga.csv"
threshold=0.05
method="spearman"
