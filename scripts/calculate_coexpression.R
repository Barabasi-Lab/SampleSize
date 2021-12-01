library(data.table)
library(dplyr)
library(optparse)
source("/Users/quim/Dropbox/Postdoc/Projects/Scipher/SampleSizeEffect/scripts/coexpression_functions.R")

#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-f", "--rnaseq_file"), type="character", 
              help="RNAseq file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="gene_correlation_pearson.txt", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-c", "--correlation_metric"), type="character", default="pearson", 
              help="correlation metric (i.e., pearson, spearman) [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /Users/quim/Dropbox/Postdoc/Projects/Scipher/SampleSizeEffect/scripts/calculate_coexpression.R -f /Users/quim/Documents/DATA/BioGroup/Scipher/data/out/Meta_RNA_NonResponders.txt -o /Users/quim/Documents/DATA/BioGroup/Scipher/data/out/NonResponders_Corr_Pearson.txt -c pearson

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$rnaseq_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (rnaseq file).n", call.=FALSE)
}
#rnaseq_file = opt$rnaseq_file
#output_file = opt$output_file
#correlation_metric = opt$correlation_metric
rnaseq_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/out/Meta_RNA_NonResponders.txt'
output_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/out/NonResponders_Corr_Spearman.txt'
correlation_metric = 'spearman'


#### READ DATASET ####

# Read RNAseq dataset (rows = genes, columns = samples)
rnaseq = fread(rnaseq_file)
rnaseq <- rnaseq[, 1:100] # Check example with less genes
# Remove the subjects column and include it as "rownames"
subject.ids <- rnaseq[, 1][[1]]
rnaseq <- as.matrix(rnaseq[, -c(1)])
rownames(rnaseq) <- subject.ids
# Get transposed dataframe
rnaseq.t <- t(as.matrix(rnaseq))
colnames(rnaseq.t) <- rownames(rnaseq)
rownames(rnaseq.t) <- colnames(rnaseq)


#### CALCULATE CORRELATION ####

# Calculate all possible combinations of pairs: n*(n-1)/2
# Calculate correlation and output result
output_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/out/NonResponders_Corr_Spearman.txt'
calculate_correlation(rnaseq, output_file, cor_method=correlation_metric)


#### CALCULATION USING WGCNA ####

# Load the WGCNA package
output_file = '/Users/quim/Documents/DATA/BioGroup/Scipher/data/out/NonResponders_Corr_WGCNA.txt'
calculate_correlation_WGCNA(rnaseq.t, out_name=output_file, type="signed", power=12)

