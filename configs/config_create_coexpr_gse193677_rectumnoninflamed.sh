#!/bin/bash

samples_folder="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GEO/GSE193677/sampling_with_repetition/inflamed_vs_noninflamed/Rectum_noninflamed"
rnaseq_file="/work/ccnr/j.aguirreplans/Databases/GEO/GSE193677/out/reads/filter_genes_low_counts/rnaseq_GSE193677.txt"
output_folder="/scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_noninflamed"
dataset="GSE193677_Rectum_noninflamed"
metric="aracne"
max_size=1000
