#!/bin/bash
#SBATCH --nodes=2
#SBATCH --mem=200000
#SBATCH --job-name=merge_networks
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_merge_networks.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_merge_networks.sh.err
#SBATCH -p netsi_standard
# Possible queues ==> short, netsi_standard

source /etc/profile.d/modules.sh
module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"

# Pearson
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_scipher/reads/scipher.complete.dataset -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_scipher/reads/scipher_scipher.complete.dataset_pearson_combined.txt -s 20 -r 1

#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Whole.Blood -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Whole.Blood_pearson_combined_20-400.txt -s 100 -r 5

#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_with_5_reps.txt -s 100 -r 5
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-LUAD -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-LUAD_pearson_combined_with_5_reps.txt -s 100 -r 5
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_with_small_sizes.txt -s 20 -r 5
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-LUAD -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-LUAD_pearson_combined_with_small_sizes.txt -s 20 -r 5

# Aracne
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_UC_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_aracne_combined.txt -m aracne -s 20 -r 5

# Genie3
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_UC_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_genie3_combined.txt -m genie3 -s 20 -r 5
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_genie3_combined.txt -m genie3 -s 20 -r 5
singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_genie3_combined.txt -m genie3 -s 20 -r 5
