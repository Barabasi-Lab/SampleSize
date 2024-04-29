#!/bin/bash

samples_folder="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GEO/GSE193677/sampling_with_repetition/disease_inflamed_vs_control/Rectum_CD_inflamed"
rnaseq_file="/work/ccnr/j.aguirreplans/Databases/GEO/GSE193677/out/reads/filter_genes_low_counts/rnaseq_GSE193677.txt"
output_folder="/scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_CD_inflamed"
dataset="GSE193677_Rectum_CD_inflamed"
metric="aracne"
max_size=1000
