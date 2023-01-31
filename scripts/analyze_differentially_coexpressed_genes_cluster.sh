#!/bin/bash
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH --mem=60000
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_analyze_differentially_coexpressed_genes.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%j_analyze_differentially_coexpressed_genes.sh.err
module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"

singularity run -B "/scratch:/scratch,/work:/work" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_differentially_coexpressed_genes.R -a /scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-BRCA_female___TCGA-Breast_female/consensus -b tcga.brca.female -c tcga.breast.female -d /home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022.csv -e breast.neoplasms -f /work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv -g /work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022_distances_matrix.txt -i /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female -j /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female -k 0.00001 -l pval_Phi_Tilde.adj.bonf -m /home/j.aguirreplans/Projects/Scipher/SampleSize/data/nodes_to_follow/breast_cancer.txt -n /home/j.aguirreplans/Projects/Scipher/SampleSize/data/drug_targets/drugbank_targets_breast.neoplasm.txt
