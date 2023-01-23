#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --job-name=separation
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.err
#SBATCH -p short

module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"
singularity run -B "/scratch:/scratch,/work:/work" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_separation.R -a /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/nodes_by_type_and_category_and_size/codina_ppi_lcc_nodes_separation_type_"$type"_size_"$size".txt -b /work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv -c /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female/separation_by_type_and_size/codina_ppi_lcc_separation_type_"$type"_size_"$size".txt

