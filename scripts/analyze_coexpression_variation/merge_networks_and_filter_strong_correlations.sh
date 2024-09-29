#!/bin/bash
#SBATCH --nodes=2
#SBATCH --mem=200000
#SBATCH --job-name=merge_networks_strong_corr
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_merge_networks_strong_corr.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_merge_networks_strong_corr.sh.err
#SBATCH -p short

source /etc/profile.d/modules.sh
module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"

singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_filtered_by_strong_correlations.txt -m pearson -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_UC_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_pearson_combined_filtered_by_strong_correlations.txt -m pearson -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_pearson_combined_filtered_by_strong_correlations.txt -m pearson -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_pearson_combined_filtered_by_strong_correlations.txt -m pearson -s 20 -r 5 -t 0.6
