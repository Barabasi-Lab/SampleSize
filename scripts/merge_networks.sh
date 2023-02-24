#!/bin/bash
#SBATCH --mem=200000
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH --constraint="[cascadelake|zen2]"
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_merge_networks.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_merge_networks.sh.err
module load R/4.0.3
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_scipher/reads/scipher.complete.dataset -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_scipher/reads/scipher_scipher.complete.dataset_pearson_combined.txt -s 20 -r 1
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Whole.Blood -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Whole.Blood_pearson_combined.txt -s 20 -r 1
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Whole.Blood -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Whole.Blood_pearson_combined_with_5_reps.txt -s 100 -r 5
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_with_5_reps.txt -s 100 -r 5
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-LUAD -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-LUAD_pearson_combined_with_5_reps.txt -s 100 -r 5
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_with_small_sizes.txt -s 20 -r 5
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-LUAD -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-LUAD_pearson_combined_with_small_sizes.txt -s 20 -r 5
