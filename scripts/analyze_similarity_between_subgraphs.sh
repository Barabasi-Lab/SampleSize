#!/bin/bash
#SBATCH --mem=50000
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH --constraint="[cascadelake|zen2]"
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analyze_similarity_between_subgraphs.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analyze_similarity_between_subgraphs.sh.err
module load R/4.0.3
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_similarity_between_subgraphs.R

