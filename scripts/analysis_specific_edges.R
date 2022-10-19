#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
require(magrittr)
library(optparse)
library(tidyr)
set.seed(1510)
`%ni%` <- Negate(`%in%`)
options(bitmapType='cairo')


#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-n", "--networks_dir"), type="character", 
              help="Co-expression networks directory", metavar="character"),
  make_option(c("-o", "--output_file"), type="character",
              help="Output file with results of selected links", metavar="character"),
  make_option(c("-m", "--method"), type="character", default="pearson",
              help="Gene co-expression network calculation method [default= %default]", metavar="character"),
  make_option(c("-l", "--num_links_selected"), type="integer", default=100,
              help="Number of links selected [default= %default]", metavar="integer")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_specific_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/reads/Whole.Blood

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$networks_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (networks_dir).n", call.=FALSE)
}

networks_dir = opt$networks_dir
output_file = opt$output_file
method = opt$method
num_links_selected = as.numeric(opt$num_links_selected)

#networks_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/reads/Whole.Blood"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/analysis_specific_edges_pearson_gtex_Whole.Blood.txt"
#num_links_selected = 100
#method = "pearson"

# Read network names
cols = c("name", "method", "size", "rep")
network_info_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
network_files = list.files(networks_dir, include.dirs = FALSE)
for (network_file in network_files){
  network_path = paste(networks_dir, network_file, sep="/")
  if(file.exists(network_path) && !dir.exists(network_path)){
    file_split = strsplit(gsub(".net", "", network_file), split="_")[[1]]
    if (("threshold" %in% file_split) | ("edges" %in% file_split)){
      method = file_split[1]
      threshold = as.numeric(file_split[(length(file_split))])
      size = as.numeric(file_split[(length(file_split)-4)])
      r = as.numeric(file_split[(length(file_split))-2])
    } else {
      method = file_split[1]
      size = as.numeric(file_split[(length(file_split)-2)])
      r = as.numeric(file_split[(length(file_split))])
    }
    network_info_df = rbind(network_info_df, data.frame(name=network_file, method=method, size=size, rep=r))
  }
}

# Read largest size network
max_size_network = (network_info_df %>% filter(size == max(size) & rep == "1" & method == method))$name
network_df = fread(paste(networks_dir, max_size_network, sep="/"), header=T)

# Select randomly links from different correlation levels
network_data_df = (network_df %>% filter(abs(score) < 0.2)) %>% sample_n(size=num_links_selected, replace = FALSE)
network_data_df = rbind(network_data_df, (network_df %>% filter((abs(score) >= 0.2) & (abs(score) < 0.4))) %>% sample_n(size=num_links_selected, replace = FALSE))
network_data_df = rbind(network_data_df, (network_df %>% filter((abs(score) >= 0.4) & (abs(score) < 0.6))) %>% sample_n(size=num_links_selected, replace = FALSE))
network_data_df = rbind(network_data_df, (network_df %>% filter((abs(score) >= 0.6) & (abs(score) < 0.8))) %>% sample_n(size=num_links_selected, replace = FALSE))
network_data_df = rbind(network_data_df, (network_df %>% filter(abs(score) >= 0.8)) %>% sample_n(size=num_links_selected, replace = FALSE))
network_data_df$method = method
network_data_df$size = max(network_info_df$size)
network_data_df$rep = 1
rm(network_df)
selected_links = network_data_df %>% select(Node.1, Node.2)
selected_links$link = paste(selected_links$Node.1, selected_links$Node.2)

# Read all networks and get correlation values and p-values from selected links
for(selected_size in sort(unique((network_info_df %>% filter(method == method))$size))){
  for (selected_rep in sort(unique((network_info_df %>% filter(method == method))$rep))){
    if (!((selected_size == max(network_info_df$size)) & (selected_rep == 1))){
      start_time <- Sys.time()
      print(paste("Reading network of size ", selected_size, " and rep ", selected_rep, sep=""))
      network_name = (network_info_df %>% filter(size == selected_size & rep == selected_rep & method == method))$name
      network_filt_df = fread(paste(networks_dir, network_name, sep="/"), header=T) %>% 
        filter((Node.1 %in% selected_links$Node.1) & (Node.2 %in% selected_links$Node.2)) %>% 
        filter(paste(Node.1, Node.2) %in% selected_links$link)
      network_data_df = rbind(network_data_df, cbind(network_filt_df, data.frame(method=method, size=selected_size, rep=selected_rep)))
      rm(network_filt_df)
      end_time <- Sys.time()
      time_diff = end_time - start_time
      print(time_diff)
    }
  }
}

network_data_df %>% fwrite(output_file)


