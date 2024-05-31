#!/bin/bash

samples_folder="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/TCGA/2022-11-18-Dataset/sampling_with_repetition/tissue/Breast_female"
rnaseq_file="/work/ccnr/j.aguirreplans/Databases/TCGA/2022-11-18-Dataset/TCGA/out/reads/tissue/filter_genes_low_counts/rnaseq_filtered_files_by_sample_group/tcga_rnaseq_Breast_female.csv"
output_folder="/scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tissue/Breast_female"
dataset="tcga_Breast_female"
metric="aracne"
max_size=1200
