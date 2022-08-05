#!/usr/bin/env Rscript
library("optparse")
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


# Read RNAseq file
rnaseq_file = '/Users/quim/Documents/DATA/BioGroup/ScipherMedicine/data/out/Meta_RNA_Cleaned.csv'
#rnaseq_file = opt$rnaseq_file
metadata_RNA <- fread(rnaseq_file)

# Read samples file
samples_file <- '/Users/quim/Documents/DATA/BioGroup/ScipherMedicine/data/subsampling/subsampling_random_labels_equal_proportion.txt'
#samples_file = opt$samples_file

# Read output directory
output_dir <- '.'
#output_dir = opt$output_dir

# Read parameters
sample_size = opt$sample_size
repetition_number = opt$repetition_number

# Get information from sample
samples_df <- fread(samples_file)
subsample_df <- samples_df %>% dplyr::filter(size==sample_size, rep==repetition_number) %>% dplyr::select(-size, -rep)
subsample <- colnames(subsample_df[, which(subsample_df==1), with=FALSE])
subsample_info_df <- data.frame(matrix(ncol=2, nrow=0))
colnames(subsample_info_df) <- c('id', 'label')
for(id in subsample){
  id_info <- unlist(strsplit(id, '_'))
  subsample_info_df <- rbind(subsample_info_df, data.frame(id=id_info[1], label=id_info[2]))
}

# Divide gene expression according to labels (R3M, R6M, NR3M, NR6M)
Responder_3m = metadata_RNA %>% filter(Subject %in% subsample_info_df[subsample_info_df$label=='R3m',]$id)
Responder_6m = metadata_RNA %>% filter(Subject %in% subsample_info_df[subsample_info_df$label=='R6m',]$id)
NResponder_3m = metadata_RNA %>% filter(Subject %in% subsample_info_df[subsample_info_df$label=='NR3m',]$id)
NResponder_6m = metadata_RNA %>% filter(Subject %in% subsample_info_df[subsample_info_df$label=='NR6m',]$id)

R3 = Responder_3m[, -c(1:52)] %>% t %>% as.data.frame() 
R6 = Responder_6m[, -c(1:52)] %>% t %>% as.data.frame() 
nR3 = NResponder_3m[, -c(1:52)] %>% t %>% as.data.frame() 
nR6 = NResponder_6m[, -c(1:52)] %>% t%>% as.data.frame() 

source("./99_newwTO.R")
wto_N = 1000

prepare_wTO <- function(input, save, N=1000){
  wto = wTO.faster(Data = input, n = N)
  wto %>% fwrite(save)
  return(wto)
}

R3_wto <- prepare_wTO(input=R3, save=paste(output_dir, paste(paste('Responder_3m', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'), N=wto_N)
R6_wto <- prepare_wTO(input=R6, save=paste(output_dir, paste(paste('Responder_6m', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'), N=wto_N)
nR3_wto <- prepare_wTO(input=nR3, save=paste(output_dir, paste(paste('NonResponder_3m', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'), N=wto_N)
nR6_wto <- prepare_wTO(input=nR6, save=paste(output_dir, paste(paste('NonResponder_6m', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'), N=wto_N)


# Calculate consensus networks for responders and non-responders
source("./99_Consensus.R")
R = wTO.Consensus(data = list(R3_wto[R3_wto$pval.adj < 0.001,1:4], 
                              R6_wto[R6_wto$pval.adj < 0.001, 1:4]))
fwrite(R, paste(output_dir, paste(paste('Responder_CN', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'))

nR = wTO.Consensus(data = list(nR3_wto[nR3_wto$pval.adj < 0.001,1:4], 
                               nR6_wto[nR6_wto$pval.adj < 0.001, 1:4]))
fwrite(nR, paste(output_dir, paste(paste('NonResponder_CN', sample_size, repetition_number, sep='_'), '.csv', sep=''), sep='/'))
