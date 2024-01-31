#!/usr/bin/env Rscript
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
              help="output file [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /Users/quim/Dropbox/Postdoc/Projects/Scipher/SampleSizeEffect/scripts/01_Generate_CN_with_WGCNA.R -f /Users/quim/Documents/DATA/BioGroup/Scipher/data/bootstrap/nonresponders/RNAseq_NonResponders_sample_10_rep_1.txt -o /Users/quim/Documents/DATA/BioGroup/Scipher/data/bootstrap_WGCNA/RNAseq_NonResponders_sample_10_rep_1.net

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
#rnaseq_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/bootstrap/nonresponders/RNAseq_NonResponders_all.txt'
#output_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/out/NonResponders_Corr_Spearman.txt'

# Get scripts dir
initial.options <- commandArgs(trailingOnly = FALSE)
scripts.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

# Source needed files
newwTO.file <- paste(scripts.dir, "99_newwTO.R", sep="/")
source(newwTO.file)


#### READ DATASET ####

# Read RNAseq dataset (rows = genes, columns = samples)
rnaseq = fread(rnaseq_file)
rnaseq <- rnaseq[, 1:100] # Check example with less genes
# Remove the subjects column and include it as "rownames"
subject.ids <- rnaseq[, 1][[1]]
rnaseq <- as.matrix(rnaseq[, -c(1)])
rownames(rnaseq) <- subject.ids
# Get transposed dataframe
rnaseq.t <- t(as.matrix(rnaseq)) %>% as.data.frame() 
colnames(rnaseq.t) <- rownames(rnaseq)
rownames(rnaseq.t) <- colnames(rnaseq)


#### RUN WGCNA ####

calculate_correlation_WGCNA(rnaseq, output_file, type="signed", power=12)

