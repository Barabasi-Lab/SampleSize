#!/bin/bash
#SBATCH --mem=100000
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analysis_specific_edges.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_analysis_specific_edges.sh.err
module load R/4.0.3
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_specific_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/reads/Whole.Blood -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/analysis_specific_edges_pearson_gtex_Whole.Blood.txt -l 100 -m pearson
