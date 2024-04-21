#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --job-name=gene_coexpression
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.err
#SBATCH -p short

source /etc/profile.d/modules.sh
module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"

# Create dirs if necessary
mkdir -p "$output_folder"

# Run script
singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/gene_coexpression_networks/create_gene_coexpression_network_from_samples_list.R -s "$samples_folder"/"$dataset"_size_"$size"_rep_"$rep".txt -f "$rnaseq_file" -o "$output_folder"/"$metric"_"$dataset"_size_"$size"_rep_"$rep".net -m "$metric" -n "$wto_n" -d "$wto_delta" -p "$wgcna_power" -t "$wgcna_type" -e "$mi_estimator" -a "$aracne_eps" -c "$correction_method"
