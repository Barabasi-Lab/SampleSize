#!/usr/bin/env Rscript
renv::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
require(dplyr)
require(magrittr)
require(data.table)
set.seed(1510)

#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-s", "--samples_file"), type="character", 
              help="Samples list file. A file with the samples that will be used from the RNAseq file, each one in a different line.", metavar="character"),
  make_option(c("-f", "--rnaseq_file"), type="character", 
              help="RNAseq file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="gene_coexpression_network_WGCNA.net", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-m", "--metric"), type="character", default="wto", 
              help="metric (e.g., pearson, spearman, wto, wgcna) [default= %default]", metavar="character"),
  make_option(c("-n", "--wto_n"), type="integer", default=100, 
              help="Number of wTO bootstrap repetitions to calculate the p-value [default= %default]", metavar="integer")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network_from_samples_list.R -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Liver/RNAseq_samples_Liver_female_size_10_rep_1.txt -f /home/j.aguirreplans/Databases/GTEx/v8/GTEx_RNASeq_gene_tpm_filtered_t.gct -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Liver/RNAseq_Liver_female_size_10_rep_1.net -m wto

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

samples_file = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Liver/RNAseq_samples_Liver_female_size_10_rep_1.txt'
rnaseq_file = '/home/j.aguirreplans/Databases/GTEx/v8/GTEx_RNASeq_gene_tpm_filtered_t.gct'
output_file = '/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Liver/RNAseq_Liver_female_size_10_rep_1.net'
metric = 'wto'
wto_n = 100

# Get scripts dir
initial.options <- commandArgs(trailingOnly = FALSE)
scripts.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

# Source needed files
newwTO.file <- paste(scripts.dir, "99_newwTO.R", sep="/")
coexpression.functions.file <- paste(scripts.dir, "coexpression_functions.R", sep="/")
source(newwTO.file)
source(coexpression.functions.file)


#### READ DATASET ####

# Read samples
subsample = fread(samples_file)[, 1][[1]]

# Read RNAseq dataset (rows = samples, columns = genes)
rnaseq = fread(rnaseq_file)
#rnaseq <- rnaseq[, 1:100] # Check example with less genes

# Subset gene expression datast by samples in the samples file
rnaseq = rnaseq %>% filter(get(colnames(rnaseq)[1]) %in% subsample)

# Remove the samples column and include it as "rownames"
sample.ids <- rnaseq[, 1][[1]]
rnaseq <- as.matrix(rnaseq[, -c(1)])
rownames(rnaseq) <- sample.ids

# Get transposed dataframe
rnaseq.t <- t(as.matrix(rnaseq)) %>% as.data.frame()
colnames(rnaseq.t) <- rownames(rnaseq)
rownames(rnaseq.t) <- colnames(rnaseq)


#### RUN CO-EXPRESSION ####

if((metric == 'spearman') | (metric == 'pearson')){
  calculate_correlation(rnaseq.t, output_file, cor_method=metric)
} else if(metric == 'wgcna'){
  calculate_correlation_WGCNA(rnaseq, output_file, type="signed", power=12)
} else if(metric == 'wto'){
  # Function to calculate network using wTO faster
  prepare_wTO <- function(input, save, N=1000){
    wto = wTO.faster(Data = input, n = N)
    wto %>% fwrite(save)
    return(wto)
  }
  # Run wTO
  nR_wto <- prepare_wTO(input=rnaseq.t, save=output_file, N=wto_n)
} else {
  stop('The metric introduced is not correct. Please introduce: pearson, spearman, wgcna or wto.')
}

