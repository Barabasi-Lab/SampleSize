#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --job-name=go_enrichment_coexpr
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.err
#SBATCH -p short

module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"
#singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_go_enrichment.R -a /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/nodes_by_type_and_category_and_size/codina_coexpr_nodes_type_"$type"_cat_"$cat"_size_"$size".txt -b /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/nodes_by_type_and_category_and_size/codina_coexpr_background_nodes.txt -c /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/functions_by_type_and_category_and_size/codina_coexpr_go_enrichment_type_"$type"_cat_"$cat"_size_"$size".txt -d BP
singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_go_enrichment.R -a /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/GSE193677-inflamed___GSE193677-noninflamed/nodes_by_type_and_category_and_size/codina_coexpr_nodes_type_"$type"_cat_"$cat"_size_"$size".txt -b /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/GSE193677-inflamed___GSE193677-noninflamed/nodes_by_type_and_category_and_size/codina_coexpr_background_nodes.txt -c /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/GSE193677-inflamed___GSE193677-noninflamed/functions_by_type_and_category_and_size/codina_coexpr_go_enrichment_type_"$type"_cat_"$cat"_size_"$size".txt -d BP
