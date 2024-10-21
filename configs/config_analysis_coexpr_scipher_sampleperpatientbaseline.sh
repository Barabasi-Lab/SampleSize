#!/bin/bash

input_folder="/scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_scipher/reads/scipher.sample.per.patient.baseline"
output_results_dir="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_scipher/reads/scipher.sample.per.patient.baseline"
output_subgraphs_dir="/scratch/j.aguirreplans/Scipher/SampleSizeProject/subgraphs_scipher/reads/scipher.sample.per.patient.baseline"
genes_dataset_file="/work/ccnr/j.aguirreplans/data/Scipher/Dec2021/00_data/genes_from_rnaseq_filtered_files_by_group/scipher_rnaseq_gene_info_scipher.sample.per.patient.baseline.csv"
ppi_file="/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv"
disease_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022_scipher.csv"
essential_genes_file="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/essential_genes/OGEE_essential_genes_scipher.csv"
threshold=0.05
method="spearman"
