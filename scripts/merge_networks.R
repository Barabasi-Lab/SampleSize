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
  make_option(c("-o", "--output_dataframe"), type="character",
              help="Output dataframe containing all correlations", metavar="character"),
  make_option(c("-s", "--step"), type="character", default=20,
              help="Step between sample sizes", metavar="integer"),
  make_option(c("-r", "--max_rep"), type="character", default=1,
              help="Maximum number of repetitions", metavar="integer")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/merge_networks.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA -m pearson -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA_combined.txt
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
output_dataframe = opt$output_dataframe
step = as.integer(opt$step)
max_rep = as.integer(opt$max_rep)

# Examples
#networks_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/TCGA"
#networks_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/complete.dataset"
#networks_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/reads/Whole.Blood"
#method = "pearson"
#step = 100
#step = 20
#step = 40
#max_rep = 1
#max_rep = 2
#max_rep = 5
selected_sizes = c(20, 100, 200, 300, 400, 500)

#### PARSE NETWORKS AND CREATE UNIQUE DATAFRAME ####

result_files = list.files(networks_dir)
result_files_df = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("file_name", "size", "rep"))))

for (result_file in result_files){
  file_split = strsplit(gsub(".net", "", result_file), split="_")[[1]]
  if(file_split[1] == method){
    if(length(file_split) == 7){
      size = as.integer(file_split[length(file_split)-2])
      rep = as.integer(file_split[length(file_split)])
      result_files_df = rbind(result_files_df, data.frame(file_name=result_file, name=paste(size,rep,sep="."), size=size, rep=rep))
    }
  }
}

#selected_files_df = result_files_df %>% filter((size %% step == 0) & (rep %in% seq(from=1, to=max_rep, by=1))) %>% arrange(size, rep)
selected_files_df = result_files_df %>% filter((size %in% selected_sizes) & (rep %in% seq(from=1, to=max_rep, by=1))) %>% arrange(size, rep)

for (i in 1:nrow(selected_files_df)){
    result_file = selected_files_df[i,]$file_name
    result_name = selected_files_df[i,]$name
    print(result_file)
    if(i == 1){
      nets_df = fread(paste(networks_dir, result_file, sep="/")) %>% select(Node.1, Node.2, score) %>% rename(!!result_name := "score")
    }else{
      nets_df = cbind(nets_df, (fread(paste(networks_dir, result_file, sep="/")) %>% select(score) %>% rename(!!result_name := "score")))
    }
}

nets_df %>% fwrite(output_dataframe)
