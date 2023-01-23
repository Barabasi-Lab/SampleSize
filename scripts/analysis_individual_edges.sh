#!/bin/bash
#SBATCH --mem=200000
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analysis_individual_edges.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analysis_individual_edges.sh.err
module load R/4.0.3
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_individual_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/reads/gtex_Whole.Blood_pearson_combined.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots -f gtex_Whole.Blood_pearson
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_individual_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/reads/tumor/tcga_TCGA-BRCA_pearson_combined.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots -f tcga_TCGA-BRCA_pearson
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_individual_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/reads/scipher_scipher.complete.dataset_pearson_combined.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots -f scipher_scipher.complete.dataset_pearson
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_sd_coexpression_weight.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/reads/gtex_Whole.Blood_pearson_combined_with_5_reps.txt -t /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/analysis_sd_coexpression_weight.txt -p /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/analysis_sd_coexpression_weight.png
