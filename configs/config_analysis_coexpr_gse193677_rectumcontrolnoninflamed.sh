#!/bin/bash

input_folder="/scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_Control_noninflamed"
output_results_dir="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_Control_noninflamed"
output_subgraphs_dir="/scratch/j.aguirreplans/Scipher/SampleSizeProject/subgraphs_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_Control_noninflamed"
genes_dataset_file="/work/ccnr/j.aguirreplans/Databases/GEO/GSE193677/out/reads/filter_genes_low_counts/rnaseq_gene_info_GSE193677.csv"
ppi_file="/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv"
disease_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_GSE193677.csv"
essential_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/essential_genes/OGEE_essential_genes_GSE193677.csv"
threshold=0.05
method="spearman"
