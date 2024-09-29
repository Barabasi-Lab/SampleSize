#!/bin/bash
#SBATCH --nodes=2
#SBATCH --mem=200000
#SBATCH --job-name=analysis_sd_coexpression_weight
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_analysis_sd_coexpression_weight.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_analysis_sd_coexpression_weight.sh.err
#SBATCH -p short

source /etc/profile.d/modules.sh
module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"

singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_sd_coexpression_weight.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_spearman_combined_filtered_by_strong_correlations.txt -t /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/analysis_sd_coexpression_weight_filtered_by_strong_correlations_tcga_TCGA-BRCA_female_spearman.txt  -p /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/analysis_sd_coexpression_weight_filtered_by_strong_correlations_tcga_TCGA-BRCA_female.png
