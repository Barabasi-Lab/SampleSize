#!/bin/bash
#SBATCH --nodes=2
#SBATCH --mem=200000
#SBATCH --job-name=merge_networks_strong_corr
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_merge_networks_strong_corr.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_merge_networks_strong_corr.sh.err
#SBATCH -p netsi_standard
# Possible queues ==> short, netsi_standard

source /etc/profile.d/modules.sh
module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"

# Pearson
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_filtered_by_threshold_0.6.txt -m pearson -a 20 -b 120 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_UC_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_pearson_combined_filtered_by_threshold_0.6.txt -m pearson -a 20 -b 120 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_pearson_combined_filtered_by_threshold_0.6.txt -m pearson -a 20 -b 120 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_pearson_combined_filtered_by_threshold_0.6.txt -m pearson -a 20 -b 120 -s 20 -r 5 -t 0.6
singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_pearson_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_pearson_combined_filtered_by_threshold_0.6.txt -m pearson -a 20 -b 120 -s 20 -r 5 -t 0.6

#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_pearson_combined_filtered_by_threshold_0.6_140-240.txt -m pearson -a 140 -b 240 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_pearson_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_pearson_combined_filtered_by_threshold_0.6_140-240.txt -m pearson -a 140 -b 240 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_pearson_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_pearson_combined_filtered_by_threshold_0.6_140-240.txt -m pearson -a 140 -b 240 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_pearson_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_pearson_combined_filtered_by_threshold_0.6_140-240.txt -m pearson -a 140 -b 240 -s 20 -r 5 -t 0.6

# Spearman
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_spearman_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_spearman_combined_filtered_by_threshold_0.6.txt -m spearman -a 20 -b 120 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_spearman_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_spearman_combined_filtered_by_threshold_0.6.txt -m spearman -a 20 -b 120 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_UC_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_spearman_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_spearman_combined_filtered_by_threshold_0.6.txt -m spearman -a 20 -b 120 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_spearman_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_spearman_combined_filtered_by_threshold_0.6.txt -m spearman -a 20 -b 120 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_spearman_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_spearman_combined_filtered_by_threshold_0.6.txt -m spearman -a 20 -b 120 -s 20 -r 5 -t 0.6

#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_spearman_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_spearman_combined_filtered_by_threshold_0.6_140-240.txt -m spearman -a 140 -b 240 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_spearman_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_spearman_combined_filtered_by_threshold_0.6_140-240.txt -m spearman -a 140 -b 240 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_spearman_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_spearman_combined_filtered_by_threshold_0.6_140-240.txt -m spearman -a 140 -b 240 -s 20 -r 5 -t 0.6
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_spearman_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_spearman_combined_filtered_by_threshold_0.6_140-240.txt -m spearman -a 140 -b 240 -s 20 -r 5 -t 0.6

# Aracne
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_UC_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_aracne_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_aracne_combined_filtered_by_threshold_1.txt -m aracne -a 20 -b 120 -s 20 -r 5 -t 1
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_aracne_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_aracne_combined_filtered_by_threshold_1.txt -m aracne -a 20 -b 120 -s 20 -r 5 -t 1
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_aracne_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_aracne_combined_filtered_by_threshold_1.txt -m aracne -a 20 -b 120 -s 20 -r 5 -t 1
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_aracne_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_aracne_combined_filtered_by_threshold_1.txt -m aracne -a 20 -b 120 -s 20 -r 5 -t 1
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_aracne_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_aracne_combined_filtered_by_threshold_1.txt -m aracne -a 20 -b 120 -s 20 -r 5 -t 1

#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_aracne_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_aracne_combined_filtered_by_threshold_1_140-240.txt -m aracne -a 140 -b 240 -s 20 -r 5 -t 1
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_aracne_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_aracne_combined_filtered_by_threshold_1_140-240.txt -m aracne -a 140 -b 240 -s 20 -r 5 -t 1
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_aracne_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_aracne_combined_filtered_by_threshold_1_140-240.txt -m aracne -a 140 -b 240 -s 20 -r 5 -t 1
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_aracne_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_aracne_combined_filtered_by_threshold_1_140-240.txt -m aracne -a 140 -b 240 -s 20 -r 5 -t 1

# Genie3
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/Rectum_UC_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_genie3_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/disease_inflamed_vs_control/gse193677_Rectum_UC_inflamed_genie3_combined_filtered_by_threshold_0.01.txt -m genie3 -a 20 -b 120 -s 20 -r 5 -t 0.01
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_genie3_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_genie3_combined_filtered_by_threshold_0.01.txt -m genie3 -a 20 -b 120 -s 20 -r 5 -t 0.01
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_genie3_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_genie3_combined_filtered_by_threshold_0.01.txt -m genie3 -a 20 -b 120 -s 20 -r 5 -t 0.01
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_genie3_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_genie3_combined_filtered_by_threshold_0.01.txt -m genie3 -a 20 -b 120 -s 20 -r 5 -t 0.01

#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/TCGA-BRCA_female -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_genie3_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_tcga/reads/tumor/tcga_TCGA-BRCA_female_genie3_combined_filtered_by_threshold_0.01_140-240.txt -m genie3 -a 140 -b 240 -s 20 -r 5 -t 0.01
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Breast.Mammary.Tissue -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_genie3_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Breast.Mammary.Tissue_genie3_combined_filtered_by_threshold_0.01_140-240.txt -m genie3 -a 140 -b 240 -s 20 -r 5 -t 0.01
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/Lung -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_genie3_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_gtex/reads/gtex_Lung_genie3_combined_filtered_by_threshold_0.01_140-240.txt -m genie3 -a 140 -b 240 -s 20 -r 5 -t 0.01
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_variation/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/Rectum_inflamed -o /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_genie3_combined_140-240.txt -f /scratch/j.aguirreplans/Scipher/SampleSizeProject/networks_geo/GSE193677/reads/inflamed_vs_noninflamed/gse193677_Rectum_inflamed_genie3_combined_filtered_by_threshold_1_140-240.txt -m genie3 -a 140 -b 240 -s 20 -r 5 -t 0.01
