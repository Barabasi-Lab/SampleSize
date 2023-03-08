#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --job-name=proximity
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.err
#SBATCH -p short

module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"
singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_proximity_multiple_drugs.R -a /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/"$folder"/nodes_by_type_and_category_and_size/codina_ppi_lcc_nodes_type_"$type"_cat_"$cat"_size_"$size".txt -b /home/j.aguirreplans/Projects/Scipher/SampleSize/data/drug_targets/drugbank_targets_inflammatory.bowel.disease.txt -c /work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv -d /work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022_distances_matrix.txt -e /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/"$folder"/proximity_by_type_and_category_and_size/codina_ppi_lcc_proximity_type_"$type"_cat_"$cat"_size_"$size".txt -f 1000
