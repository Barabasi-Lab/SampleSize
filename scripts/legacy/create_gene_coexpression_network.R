#!/usr/bin/env Rscript
renv::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
require(dplyr)
require(magrittr)
require(data.table)
set.seed(1510)

#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-f", "--rnaseq_file"), type="character", 
              help="RNAseq file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="gene_coexpression_network_WGCNA.net", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-m", "--metric"), type="character", default="wto", 
              help="metric (e.g., pearson, spearman, wto, wgcna) [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /Users/quim/Dropbox/Postdoc/Projects/Scipher/SampleSizeEffect/scripts/create_gene_coexpression_network.R -f /Users/quim/Documents/DATA/BioGroup/Scipher/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /Users/quim/Documents/DATA/BioGroup/Scipher/data/bootstrap_spearman/RNAseq_NonResponders_sample_10_rep_1.net -m spearman

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$rnaseq_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (rnaseq file).n", call.=FALSE)
}
rnaseq_file = opt$rnaseq_file
output_file = opt$output_file
metric = opt$metric
#rnaseq_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt'
#output_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/bootstrap_WGCNA/RNAseq_NonResponders_sample_10_rep_1.net'
#metric = 'spearman'

# Get scripts dir
initial.options <- commandArgs(trailingOnly = FALSE)
scripts.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

# Source needed files
newwTO.file <- paste(scripts.dir, "99_newwTO.R", sep="/")
coexpression.functions.file <- paste(scripts.dir, "coexpression_functions.R", sep="/")
source(newwTO.file)
source(coexpression.functions.file)


#### READ DATASET ####

# Read RNAseq dataset (rows = samples, columns = genes)
rnaseq = fread(rnaseq_file)
#rnaseq <- rnaseq[, 1:100] # Check example with less genes
# Remove the subjects column and include it as "rownames"
subject.ids <- rnaseq[, 1][[1]]
rnaseq <- as.matrix(rnaseq[, -c(1)])
rownames(rnaseq) <- subject.ids
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
  nR_wto <- prepare_wTO(input=rnaseq.t, save=output_file, N=100)
} else {
  stop('The metric introduced is not correct. Please introduce: pearson, spearman, wgcna or wto.')
}

