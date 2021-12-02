#!/bin/bash
#SBATCH --mem=50000
#SBATCH -p debug
#SBATCH --time=0:20:00
#SBATCH --constraint="[cascadelake|zen2]"
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_generate_cn_test.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_generate_cn_test.sh.err
module load R/4.0.3
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/01_Generate_CN_with_wTO.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/Meta_RNA_NonResponders.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/wTO/NonResponders_Corr_wTO.txt
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_pearson/example.net -m pearson
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_wgcna/example.net -m wgcna
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_wto/example.net -m wto
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network_from_samples_list.R -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Liver/RNAseq_samples_Liver_female_size_10_rep_1.txt -f /home/j.aguirreplans/Databases/GTEx/v8/GTEx_RNASeq_gene_tpm_filtered_t.gct -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/RNAseq_Liver_female_size_10_rep_1_genes_100.net -m wto -n 10
