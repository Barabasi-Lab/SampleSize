#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
require(dplyr)
require(magrittr)
require(data.table)
require(wTO)
set.seed(1510)

#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-s", "--samples_file"), type="character", 
              help="Samples list file. A file with the samples that will be used from the RNAseq file, each one in a different line.", metavar="character"),
  make_option(c("-f", "--rnaseq_file"), type="character", 
              help="RNAseq file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="gene_coexpression_network_wto.net", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-m", "--metric"), type="character", default="wto", 
              help="metric (e.g., pearson, spearman, mutual_information, wto, wgcna, aracne) [default= %default]", metavar="character"),
  make_option(c("-n", "--wto_n"), type="integer", default=100, 
              help="Number of wTO bootstrap repetitions to calculate the p-value [default= %default]", metavar="integer"),
  make_option(c("-d", "--wto_delta"), type="double", default=0.05,
              help="Value that defines the interval of confidence from which the p-values of the bootstrap repetitions are calculated [default= %default]", metavar="double"),
  make_option(c("-p", "--wgcna_power"), type="integer", default=6, 
              help="Number of wTO bootstrap repetitions to calculate the p-value [default= %default]", metavar="integer"),
  make_option(c("-t", "--wgcna_type"), type="character", default="signed", 
              help="Type of network to create with WGCNA [default= %default]", metavar="character"),
  make_option(c("-e", "--mi_estimator"), type="character", default="pearson", 
              help="Estimator used to calculate the mutual information (if metric is mutual_information or aracne) [default= %default]", metavar="character"),
  make_option(c("-a", "--aracne_eps"), type="double", default=0, 
              help="Eps is a numeric value indicating the threshold used when removing an edge in ARACNE [default= %default]", metavar="double"),
  make_option(c("-c", "--correction_method"), type="character", default="bonferroni", 
              help="Method to correct for multiple testing if we use pearson or spearman as metric [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network_from_samples_list.R -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Liver/RNAseq_samples_Liver_female_size_10_rep_1.txt -f /home/j.aguirreplans/Databases/GTEx/v8/GTEx_RNASeq_gene_tpm_filtered_t.gct -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Liver/RNAseq_Liver_female_size_10_rep_1.net -m wgcna

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$rnaseq_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (rnaseq file).n", call.=FALSE)
}
samples_file = opt$samples_file
rnaseq_file = opt$rnaseq_file
output_file = opt$output_file
metric = opt$metric
wto_n = strtoi(opt$wto_n)
wto_delta = as.double(opt$wto_delta)
wgcna_power = strtoi(opt$wgcna_power)
wgcna_type = opt$wgcna_type
mi_estimator = opt$mi_estimator
aracne_eps = as.double(opt$aracne_eps)
correction_method = opt$correction_method

#samples_file = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Whole.Blood_female/RNAseq_samples_Whole.Blood_female_size_10_rep_1.txt'
#rnaseq_file = '/home/j.aguirreplans/Databases/GTEx/v8/tpm_filtered_files_by_tissue/GTEx_RNASeq_Whole.Blood_female.gct'
#output_file = '/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/wto_RNAseq_Whole.Blood_female_size_10_rep_1_100_genes.net'
#scripts.dir = '/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts'
#samples_file = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/Scipher/Dec2021/sampling_with_repetition/complete_dataset/scipher_complete_dataset_size_20_rep_1.txt'
#rnaseq_file = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/Dec2021/00_data/scipher_rnaseq_counts_processed.csv'
#metric = 'wgcna'
#wto_n = 100
#wto_delta = 0.2
#wgcna_power = 6
#wgcna_type = "signed"
#mi_estimator = "pearson"
#aracne_eps = 0.1

# Get scripts dir
initial.options <- commandArgs(trailingOnly = FALSE)
scripts.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

# Source needed files
coexpression.functions.file <- paste(scripts.dir, "coexpression_functions.R", sep="/")
source(coexpression.functions.file)


#### READ DATASET ####

# Read samples
subsample = fread(samples_file)[, 1][[1]]

# Read RNAseq dataset (rows = genes, columns = samples)
rnaseq = fread(rnaseq_file)

# Subset gene expression datast by samples in the samples file
rnaseq = rnaseq %>% dplyr::select(c(colnames(rnaseq)[1], all_of(subsample)))

# Remove the samples column and include it as "rownames"
gene.ids <- rnaseq[, 1][[1]]
rnaseq <- as.matrix(rnaseq[, -c(1)]) %>% as.data.frame()
rownames(rnaseq) <- gene.ids

#rnaseq = rnaseq[row.names(rnaseq) %in% sample(row.names(rnaseq), size=200, replace=FALSE),] # Check example with less genes

# Get transposed dataframe
rnaseq.t <- t(as.matrix(rnaseq)) %>% as.data.frame()
colnames(rnaseq.t) <- rownames(rnaseq)
rownames(rnaseq.t) <- colnames(rnaseq)


#### RUN CO-EXPRESSION ####

if((metric == 'spearman') | (metric == 'pearson')){
  
  #calculate_correlation(rnaseq, output_file, cor_method=metric, disparity_filter=TRUE, corr_threshold=NA, pval_threshold=NA)
  calculate_correlation_and_pvalue(rnaseq.t, output_file, cor_method=metric, correction_method=correction_method, disparity_filter=FALSE)

} else if((metric == 'mi') | (metric == 'mutual_information')){
  
  calculate_mutual_information(rnaseq.t, output_file, estimator=mi_estimator)
  
} else if(metric == 'wgcna'){
  
  calculate_network_WGCNA(rnaseq.t, output_file, type=wgcna_type, power=wgcna_power)
  
} else if(metric == 'aracne'){
  
  calculate_network_ARACNE(rnaseq.t, output_file, estimator=mi_estimator, eps=aracne_eps)
  
} else if(metric == 'wto'){
  
  # Function to calculate network using wTO faster
  prepare_wTO <- function(input, save, scripts.dir, N=100){
    newwTO.file <- paste(scripts.dir, "99_newwTO.R", sep="/")
    source(newwTO.file)
    wto = wTO.faster(Data = input, n = N)
    wto %>% fwrite(save)
    return(wto)
  }

  # Run wTO fast  
  #wto = wTO.fast(Data = rnaseq, Overlap = row.names(rnaseq), method = "p",
  #               sign = "sign", delta = wto_delta, n = wto_n,
  #               method_resampling = "Bootstrap")
  #wto %>% fwrite(output_file)
  
  # Run wTO faster
  wto = prepare_wTO(input=rnaseq, save=output_file, scripts.dir=scripts.dir, N=wto_n)
  
}




