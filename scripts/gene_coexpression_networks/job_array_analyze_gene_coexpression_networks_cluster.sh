#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=60000
#SBATCH --job-name=analyze_gene_coexpression
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/../logs/%x.%j.sh.err
#SBATCH -p short

source /etc/profile.d/modules.sh
module load singularity/3.5.3
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"

# Create dirs if necessary
mkdir -p "$output_results_dir"
mkdir -p "$output_subgraphs_dir"

# Run script
singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/gene_coexpression_networks/analyze_coexpression_network_by_significant_edges.R -c "$coexpression_network_file" -o "$output_results_dir" -s "$output_subgraphs_dir" -f "$dataset_name" -t "$threshold" -p "$ppi_file" -d "$disease_genes_file" -e "$essential_genes_file" -g "$genes_dataset_file"
