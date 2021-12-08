#!/bin/bash
#SBATCH --mem=50000
#SBATCH -p short
#SBATCH --time=24:00:00
#SBATCH --constraint="[cascadelake|zen2]"
#SBATCH -o /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_generate_cn_test.sh.out
#SBATCH -e /home/j.aguirreplans/Projects/Scipher/SampleSize/logs/%j_generate_cn_test.sh.err
module load R/4.0.3
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/01_Generate_CN_with_wTO.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/Meta_RNA_NonResponders.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/wTO/NonResponders_Corr_wTO.txt
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_pearson/example.net -m pearson
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_wgcna/example.net -m wgcna
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network.R -f /home/j.aguirreplans/Projects/Scipher/SampleSize/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_wto/example.net -m wto
#Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network_from_samples_list.R -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Liver/RNAseq_samples_Liver_female_size_10_rep_1.txt -f /home/j.aguirreplans/Databases/GTEx/v8/GTEx_RNASeq_gene_tpm_filtered_t.gct -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/RNAseq_Liver_female_size_10_rep_1_genes_100.net -m wto -n 10
Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network_with_wto_finding_pvalue_convergence.R -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Spleen_female/RNAseq_samples_Spleen_female_size_50_rep_5.txt -f /home/j.aguirreplans/Databases/GTEx/v8/tpm_filtered_files_by_tissue/GTEx_RNASeq_Spleen_female.gct -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/wto_RNAseq_samples_Spleen_female_size_50_rep_5_2000_genes.net -b /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/wto_bootstrap_RNAseq_samples_Spleen_female_size_50_rep_5_2000_genes.csv -u /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/pearson_RNAseq_samples_Spleen_female_size_50_rep_5_2000_genes.csv -n 300 -e 50 -d 0.2 -p 10 -c 0.01

