#!/usr/bin/env Rscript
library(optparse)
require(dplyr)
require(magrittr)
require(data.table)

# Read arguments
option_list = list(
  make_option(c("-f", "--rnaseq_file"), type="character", 
              help="RNAseq file name", metavar="character"),
  make_option(c("-s", "--samples_file"), type="character", 
              help="Name of the file containing the samples", metavar="character"),
  make_option(c("-z", "--sample_size"), type="integer", default=10, 
              help="Size of the sample", metavar="integer"),
  make_option(c("-r", "--repetition_number"), type="integer", default=1,
              help="Repetition number of the subsample", metavar="integer"),
  make_option(c("-o", "--output_dir"), type="character", default=".", 
              help="output directory [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /Users/quim/Dropbox/Postdoc/Projects/Scipher/SampleSizeEffect/scripts/01_Generate_CN_with_wTO.R -f /Users/quim/Documents/DATA/BioGroup/Scipher/data/out/Meta_RNA_Cleaned.csv -s /Users/quim/Documents/DATA/BioGroup/Scipher/data/subsampling/subsampling_random_labels_original_proportion.txt -z 10 -r 1 -o /Users/quim/Documents/DATA/BioGroup/Scipher/data/subsampling/wTO

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rnaseq_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (rnaseq file).n", call.=FALSE)
}
if (is.null(opt$samples_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (samples file).n", call.=FALSE)
}

# Get scripts dir
initial.options <- commandArgs(trailingOnly = FALSE)
scripts.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

# Source needed files
newwTO.file <- paste(scripts.dir, "99_newwTO.R", sep="/")
source(newwTO.file)

# Read RNAseq file
#rnaseq_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/out/Meta_RNA_Cleaned.csv'
rnaseq_file = opt$rnaseq_file
metadata_RNA <- fread(rnaseq_file)
labels_R = metadata_RNA$Subject[(metadata_RNA$acr50_3m == "responder")&(metadata_RNA$acr50_6m == "responder")]
labels_nR = metadata_RNA$Subject[(metadata_RNA$acr50_3m == "nonresponder")&(metadata_RNA$acr50_6m == "nonresponder")]

# Read samples file
#samples_file <- '/Users/quim/Documents/DATA/BioGroup/Scipher/data/subsampling/subsampling_random_labels_original_proportion.txt'
samples_file = opt$samples_file

# Read output directory
#output_dir <- '.'
output_dir = opt$output_dir

# Read parameters
sample_size = opt$sample_size
repetition_number = opt$repetition_number

# Get information from sample
samples_df <- fread(samples_file)
subsample_df <- samples_df %>% dplyr::filter(size==sample_size, rep==repetition_number) %>% dplyr::select(-size, -rep)
subsample <- colnames(subsample_df[, which(subsample_df==1), with=FALSE])
subsample_info_df <- data.frame(patient=subsample[subsample %in% labels_R], response='R')
subsample_info_df <- rbind(subsample_info_df, data.frame(patient=subsample[subsample %in% labels_nR], response='NR'))

# Divide gene expression according to labels (R3M, R6M, NR3M, NR6M)
Responder = metadata_RNA %>% filter(Subject %in% subsample_info_df[subsample_info_df$response=='R',]$patient)
NResponder = metadata_RNA %>% filter(Subject %in% subsample_info_df[subsample_info_df$response=='NR',]$patient)

R = Responder[, -c(1:52)] %>% t %>% as.data.frame() 
NR = NResponder[, -c(1:52)] %>% t %>% as.data.frame() 

wto_N = 1000

prepare_wTO <- function(input, save, N=1000){
  wto = wTO.faster(Data = input, n = N)
  wto %>% fwrite(save)
  return(wto)
}

R_wto <- prepare_wTO(input=R, save=paste(output_dir, paste(paste('Responder_CN', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'), N=wto_N)
nR_wto <- prepare_wTO(input=NR, save=paste(output_dir, paste(paste('NonResponder_CN', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'), N=wto_N)

