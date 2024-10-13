library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
require(magrittr)
library(optparse)
library(tidyr)
set.seed(1510)
`%ni%` <- Negate(`%in%`)


#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-n", "--networks_dir"), type="character", 
              help="Co-expression networks directory", metavar="character"),
  make_option(c("-m", "--method"), type="character", default="pearson", 
              help="Co-expression network method [default= %default]", metavar="character"),
  make_option(c("-o", "--output_merged_network"), type="character",
              help="Output file to save the merged network without filtering. If it exists, it will be used for the filtering.", metavar="character"),
  make_option(c("-f", "--output_filtered_network"), type="character",
              help="Output file to save the filtered merged network", metavar="character"),
  make_option(c("-a", "--min_size"), type="integer", default=20,
              help="Minimum sample size", metavar="integer"),
  make_option(c("-b", "--max_size"), type="integer", default=120,
              help="Maximum sample size", metavar="integer"),
  make_option(c("-s", "--step"), type="integer", default=20,
              help="Step between sample sizes", metavar="integer"),
  make_option(c("-r", "--max_rep"), type="integer", default=1,
              help="Maximum number of repetitions", metavar="integer"),
  make_option(c("-t", "--score_threshold"), type="double", default=0.6,
              help="Score threshold to select (strong) correlations", metavar="double")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks_and_filter_strong_correlations.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA_combined.txt -f /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA_combined_filtered.txt
# Memory needed for 50 files: 200000

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$networks_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (networks_dir).n", call.=FALSE)
}

networks_dir = opt$networks_dir
method = opt$method
output_merged_network = opt$output_merged_network
output_filtered_network = opt$output_filtered_network
min_size = as.integer(opt$min_size)
max_size = as.integer(opt$max_size)
step = as.integer(opt$step)
max_rep = as.integer(opt$max_rep)
score_threshold = as.numeric(opt$score_threshold)
selected_sizes <- seq(min_size, max_size, step)


#### DEFINE FUNCTIONS ####

###############################
## read_coexpression_network ##
###############################
#'  @param coexpression_network_file The path to the co-expression network file
#'  @param method The method used to create the co-expression network. 

read_coexpression_network <- function(coexpression_network_file, method){
  
  # Read file
  coexpression_df = as.data.frame(data.table::fread(coexpression_network_file))
  
  # If method is aracne, genie3 or wgcna, transform matrix to dataframe of pairs of genes
  if (method %in% c("aracne", "genie3", "wgcna")){
    rownames(coexpression_df) = colnames(coexpression_df)
    coexpression_df = coexpression_df[order(rownames(coexpression_df)), order(colnames(coexpression_df))]
    coexpression_df = coexpression_df %>%
      wTO::wTO.in.line() %>%
      dplyr::rename(score=wTO)
  } else {
    if (method == "wto"){
      coexpression_df = coexpression_df %>%
        dplyr::rename("score"= "wTO")
    }
    coexpression_df = coexpression_df %>%
      dplyr::select(Node.1, Node.2, score)
  }
  
  return(coexpression_df)
}


#### PARSE NETWORKS AND CREATE UNIQUE DATAFRAME ####

if (file.exists(output_merged_network) && !dir.exists(output_merged_network)) {
  # Read merged network
  nets_df <- data.table::fread(output_merged_network, header = TRUE)

} else {

  result_files <- list.files(networks_dir)
  result_files_df <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("file_name", "size", "rep"))))
  
  for (result_file in result_files) {
    file_split <- strsplit(gsub(".net", "", result_file), split="_")[[1]]
    if (file_split[1] == method) {
      if (length(file_split) %in% c(7, 8, 9)) {
        size <- as.integer(file_split[length(file_split)-2])
        rep <- as.integer(file_split[length(file_split)])
        result_files_df <- rbind(
          result_files_df,
          data.frame(
            file_name = result_file,
            name = paste(size,rep,sep="."),
            size = size,
            rep = rep
          )
        )
      }
    }
  }
  
  selected_files_df <- result_files_df %>%
    dplyr::filter((size %in% selected_sizes) & (rep %in% seq(from=1, to=max_rep, by=1))) %>%
    dplyr::arrange(size, rep)

  for (i in 1:nrow(selected_files_df)){
    result_file = selected_files_df[i,]$file_name
    result_name = selected_files_df[i,]$name
    print(result_file)
    
    if (i == 1) {
      nets_df <- read_coexpression_network(paste(networks_dir, result_file, sep="/"), method) %>%
        dplyr::rename(!!result_name := "score")
    } else {
      nets_df <- nets_df %>%
        dplyr::full_join(
          (read_coexpression_network(paste(networks_dir, result_file, sep="/"), method) %>%
             dplyr::rename(!!result_name := "score")),
          by = c("Node.1", "Node.2")
        )
    }
  }
  
  # Write output merged network
  nets_df %>% data.table::fwrite(output_merged_network)
}

# Filter by score
cols <- colnames(nets_df)
score_cols <- cols[!(cols %in% c("Node.1", "Node.2"))]
if (score_threshold == 0) {
  nets_df <- nets_df %>%
    dplyr::filter(dplyr::if_any(score_cols, ~ abs(.) > abs(score_threshold)))
} else {
  nets_df <- nets_df %>%
    dplyr::filter(dplyr::if_any(score_cols, ~ abs(.) >= abs(score_threshold)))
}

# Write output table
nets_df %>% data.table::fwrite(output_filtered_network)
