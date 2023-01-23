library(data.table)
library(dplyr)
library(gghalves)
library(ggplot2)
library(igraph)
require(magrittr)
library(optparse)
library(stringr)
library(tidyr)
set.seed(1510)
`%ni%` <- Negate(`%in%`)
options(bitmapType='cairo')


#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-n", "--networks_file"), type="character", 
              help="Merged co-expression network file", metavar="character"),
  make_option(c("-t", "--output_sd_table_file"), type="character",
              help="Output table file", metavar="character"),
  make_option(c("-p", "--output_sd_plot_file"), type="character", 
              help="Output plot file", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analysis_individual_edges.R -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/gtex_Whole.Blood_pearson_combined.txt -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots -f gtex_Whole.Blood_pearson
# Memory needed: 200000

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$networks_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (networks_file).n", call.=FALSE)
}

networks_file = opt$networks_file
output_sd_table_file = opt$output_sd_table_file
output_sd_plot_file = opt$output_sd_plot_file


#### SD calculation (at once) ####

networks_df = fread(networks_file, header=T)
networks_df = networks_df[row.names(networks_df) %in% sample(row.names(networks_df), size=1000000, replace=FALSE),] # Check example with less edges
networks_df = networks_df %>%
  mutate(gp = paste(Node.1, Node.2, sep = " -- ")) %>%
  select(!(c(Node.1, Node.2))) %>%
  pivot_longer(!gp, names_to = "ss.rep", values_to = "correlation") %>%
  separate("ss.rep", into=c("ss", "rep"), sep="[.]") %>%
  group_by(gp, ss) %>%
  summarise(sd_correlation=sd(correlation)) %>%
  ungroup()
networks_df %>% fwrite(output_sd_table_file)

#### SD calculation (by chunks) ####

# # Connection to read columns
# con <- file(networks_file,"r")
# first_line <- readLines(con,n=1)
# close(con)
# cols = strsplit(first_line, split=",")[[1]]
# 
# # Connection to read chunks
# chunk_size = 5000000
# x=0
# chunk_count = 0
# connection = file(networks_file, "r")
# networks_df = data.frame()
# repeat {
#   x=x+1
#   chunk_count=chunk_count+chunk_size
#   message(paste('Chunk:', chunk_count))
#   skip_param=0
#   if(x==1){
#     # Skip first line (columns)
#     skip_param=1
#   }
#   # Read chunk
#   network_chunk = read.table(connection, nrows=chunk_size, skip=skip_param, header=FALSE, fill = TRUE, sep=",", col.names=cols)
#   # Pivot longer and calculate SD for each edge and size
#   network_chunk = network_chunk %>% 
#     mutate(gp = paste(Node.1, Node.2, sep = " -- ")) %>%
#     select(!(c(Node.1, Node.2))) %>%
#     pivot_longer(!gp, names_to = "ss.rep", values_to = "correlation") %>%
#     separate("ss.rep", into=c("ss", "rep"), sep="[.]") %>%
#     mutate(across('ss', str_replace, 'X', '')) %>%
#     group_by(gp, ss) %>%
#     summarise(sd_correlation=sd(correlation)) %>%
#     ungroup()
#   # Append calculations
#   if(nrow(networks_df) == 0){
#     networks_df = data.frame(network_chunk)
#   } else {
#     networks_df = rbind(networks_df, network_chunk)
#   }
#   # Close connection when necessary
#   if(nrow(network_chunk) < chunk_size){
#     close(connection)
#     break
#   }
#   rm(network_chunk)
# }

# Plot size vs SD distribution
networks_df %>%
  mutate(ss = factor(ss, levels = as.character(sort(as.numeric(unique(networks_df$ss)))))) %>%
  ggplot(aes(x = ss, y = sd_correlation)) +
  geom_half_violin(side = "r", color="orangered4", fill="orangered1") +
  theme_linedraw() +
  theme(aspect.ratio=1, 
        plot.title =  element_text(size = 15, face="bold"), 
        axis.title = element_text(size = 14, face="bold"), 
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13), 
        legend.title=element_text(size=14, face="bold")) +
  labs (x = "Number of samples",
        y = "Co-expression weight SD",
        color = "",
        fill = "")
ggsave(
  output_sd_plot_file,
  dpi = 1200,
  #width = 9300,
  #height = 6000,
  units = c("px")
)

