#!/bin/bash
#SBATCH --mem=200000
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH --constraint="[cascadelake|zen2]"
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_merge_networks.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_merge_networks.sh.err
module load R/4.0.3
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/tcga_TCGA_pearson_combined.txt -s 200 -r 1
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/complete.dataset -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/scipher_complete.dataset_pearson_combined.txt -s 20 -r 1
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/Whole.Blood -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/gtex_Whole.Blood_pearson_combined.txt -s 20 -r 1
