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

# Create dirs if necessary
mkdir -p /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/"$folder"/nodes_by_type_and_category_and_size/
mkdir -p /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/"$folder"/functions_by_type_and_category_and_size/

# Run script
singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_go_enrichment.R -a /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/"$folder"/nodes_by_type_and_category_and_size/codina_coexpr_nodes_type_"$type"_cat_"$cat"_size_"$size".txt -b /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/"$folder"/nodes_by_type_and_category_and_size/codina_coexpr_background_nodes.txt -c /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/"$folder"/functions_by_type_and_category_and_size/codina_coexpr_go_enrichment_type_"$type"_cat_"$cat"_size_"$size".txt -d BP
