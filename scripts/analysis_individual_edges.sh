#!/bin/bash
#SBATCH --mem=200000
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH --constraint="[cascadelake|zen2]"
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analysis_individual_edges.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analysis_individual_edges.sh.err
module load R/4.0.3
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_individual_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/gtex_Whole.Blood_pearson_combined.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots -f gtex_Whole.Blood_pearson
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_individual_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/tcga_TCGA_pearson_combined.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots -f tcga_TCGA_pearson
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_individual_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/scipher_complete.dataset_pearson_combined.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots -f scipher_complete.dataset_pearson
