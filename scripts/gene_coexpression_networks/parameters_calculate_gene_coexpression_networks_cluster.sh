#!/bin/bash

config_file="$1"

# Check if the config file exists
if [ ! -f "$config_file" ]; then
    echo "Config file not found: $config_file"
    exit 1
fi

# Read variables from the config file
source "$config_file"

# Define fixed variables
min_size=20
step_size=20
repetitions=(1 2 3 4 5)
wto_n=100
wto_delta=0.2
wgcna_power=6
wgcna_type="signed"
mi_estimator="pearson"
aracne_eps=0
correction_method="bonferroni"

sizes=()
for (( i=min_size; i<=max_size; i+=step_size )); do
    sizes+=("$i")
done

for size in "${sizes[@]}"; do
    for rep in "${repetitions[@]}"; do
        samples_file="$samples_folder"/"$dataset"_size_"$size"_rep_"$rep".txt
        output_file="$output_folder"/"$metric"_"$dataset"_size_"$size"_rep_"$rep".net
        if [ -f "$samples_file" ]; then
            echo "Samples file exists: $samples_file"
            if [ -f "$output_file" ]; then
                echo "Output file exists: $output_file"
            else
                sbatch --export=size=$size,rep=$rep,samples_folder=$samples_folder,rnaseq_file=$rnaseq_file,output_folder=$output_folder,dataset=$dataset,metric=$metric,wto_n=$wto_n,wto_delta=$wto_delta,wgcna_power=$wgcna_power,wgcna_type=$wgcna_type,mi_estimator=$mi_estimator,aracne_eps=$aracne_eps,correction_method=$correction_method /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/gene_coexpression_networks/job_array_calculate_gene_coexpression_networks_cluster.sh
                #echo "Output file does not exist: $output_file"
            fi
        else
            echo "Samples file does not exist: $samples_file"
            break  # Exit the loop if RNAseq file doesn't exist
        fi
    done
done
