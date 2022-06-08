library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(igraph)
require(magrittr)
library(optparse)
library(plot3D)
library(som)
library(tidyr)
set.seed(1510)
`%ni%` <- Negate(`%in%`)
options(bitmapType='cairo')


#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-n", "--networks_file"), type="character", 
              help="Co-expression networks directory", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character",
              help="Output directory", metavar="character"),
  make_option(c("-f", "--output_file_name"), type="character", default="plot", 
              help="Start of the name for the output files [default= %default]", metavar="character")
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
output_dir = opt$output_dir
output_file_name = opt$output_file_name

#networks_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_tcga/tcga_TCGA_pearson_combined.txt"
#networks_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/gtex_Whole.Blood_pearson_combined.txt"
#networks_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_gtex/gtex_Whole.Blood_pearson_combined_dummy.txt"
#networks_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/scipher_complete.dataset_pearson_combined.txt"
#output_dir = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots'
#output_file_name = "gtex_dummy_pearson"

som_plot_file = paste(output_dir, paste(output_file_name, "som.png", sep="_"), sep="/")
hist_plot_file = paste(output_dir, paste(output_file_name, "histogram_sd.png", sep="_"), sep="/")
hist_gp_plot_file = paste(output_dir, paste(output_file_name, "histogram_sd_by_size.png", sep="_"), sep="/")
hist_correlation_distribution_file = paste(output_dir, paste(output_file_name, "distribution_score.png", sep="_"), sep="/")
hist_3d_plot_file = paste(output_dir, paste(output_file_name, "histogram_3d.png", sep="_"), sep="/")
heat_2d_plot_file = paste(output_dir, paste(output_file_name, "heatmap_2d.png", sep="_"), sep="/")

# Read networks
start_time <- Sys.time()
print("Reading network")
networks_df = fread(networks_file, header=T)
# Select randomly 1 million edges
print("Selecting random interactions from network")
networks_df = networks_df[row.names(networks_df) %in% sample(row.names(networks_df), size=1000000, replace=FALSE),] # Check example with less edges
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)


#--------------#
# SOM analysis #
#--------------#

# Normalize values
start_time <- Sys.time()
print("Normalize values for SOM algorithm")
norm = som::normalize(networks_df[,!(c("Node.1","Node.2"))], byrow=TRUE)
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# Calculate som
start_time <- Sys.time()
print("Apply SOM algorithm")
som_vars <- som::som(norm, xdim=1, ydim=5)
#max_val = ceiling(max(norm))
#min_val = floor(min(norm))
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# Plot SOM
start_time <- Sys.time()
print("Plot SOM")
som_df = som_vars$code
colnames(som_df) = colnames(networks_df[,!(c("Node.1","Node.2"))])
som_df = as.data.frame(t(som_df))
som_df$ss.rep = rownames(som_df) 
som_df %<>% separate("ss.rep", into=c("ss", "rep"), sep="[.]") %>% select(!(rep))
som_df %<>%
  group_by(ss) %>% 
  pivot_longer(cols=-ss, names_to="group", values_to="score")
som_df$ss=as.numeric(som_df$ss)
som_df$group = sub('V', '', som_df$group)
#som_sum_df = Rmisc::summarySE(som_df, measurevar="score", groupvars=c("ss")) #http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(som_df, aes(x=ss, y=score, col=group)) + 
  geom_line() +
  theme_linedraw() +
  xlab("Number of samples") +
  ylab("Normalized correlation") +
  theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom") +
  guides(col=guide_legend(title="SOM group"))
ggsave(
  som_plot_file,
  dpi = 1200,
  width = 10000,
  height = 6000,
  units = c("px")
)
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)


#-------------#
# SD analysis #
#-------------#

# Transform the networks matrix in long format, where gp is the edge, ss is sample size, score is correlation
start_time <- Sys.time()
print("Transform dataframe to long format")
networks_df = networks_df %>% 
  mutate(gp = paste(Node.1, Node.2, sep = " -- ")) %>%
  select(!(c(Node.1, Node.2))) %>%
  group_by(gp) %>%
  pivot_longer(!gp, names_to = "ss.rep", values_to = "score") %>%
  mutate(sd=sd(score)) %>%
  ungroup() %>% 
  separate("ss.rep", into=c("ss", "rep"), sep="[.]")

networks_df$ss = as.integer(networks_df$ss)
networks_df$rep = as.integer(networks_df$rep)
max_sample_cut = plyr::round_any(max(networks_df$ss), 100, f=ceiling)

end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# ssgr are groups of sample sizes by intervals of 100
# sdgr is the standard deviation of each group
start_time <- Sys.time()
print("Group by intervals")
networks_df %<>% mutate(ssgr=cut(ss, breaks= seq(0, max_sample_cut, by = 100))) %>%
  group_by(gp, ssgr) %>%
  mutate(sdgr=sd(score)) %>%
  ungroup()
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# Plot histogram of standard deviations of each edge
start_time <- Sys.time()
print("Plot histogram of standard deviation 1")
ggplot((networks_df %>% select(gp, sd) %>% unique()), aes(x=sd)) + 
  geom_histogram(binwidth = 0.01, alpha=0.8, color=2, fill=2) +
  theme_linedraw() +
  xlab("Standard deviation") +
  ylab("Number of edges") +
  theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom")
ggsave(
  hist_plot_file,
  dpi = 1200,
  width = 10000,
  height = 6000,
  units = c("px")
)
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

# Plot histogram of standard deviation of each edge by groups of sizes
start_time <- Sys.time()
print("Plot histogram of standard deviation 2")
ggplot((networks_df %>% select(gp, sdgr, ssgr) %>% unique()), aes(x=sdgr, fill=ssgr)) + 
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.01) +
  theme_linedraw() +
  xlab("Standard deviation") +
  ylab("Number of edges") +
  guides(col=guide_legend(title="Group of sizes")) +
  theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom")
ggsave(
  hist_gp_plot_file,
  dpi = 1200,
  width = 10000,
  height = 6000,
  units = c("px")
)
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)


#-----------------------#
# Distribution analysis #
#-----------------------#

# Plot histogram of distribution
start_time <- Sys.time()
print("Plot histogram of distribution")
hist_correlation_distribution = networks_df %>%
  filter(ss %% 100 == 0) %>%
  ggplot(aes(score, color=ss, group=ss)) +
  geom_freqpoly(binwidth=0.1) +
  #geom_histogram(binwidth=0.1) +
  #scale_y_log10() + 
  ggtitle("Distribution of Pearson correlation") +
  labs(x="Pearson correlation", y = "Number of edges", color = "Number of samples") +
  theme_linedraw() + 
  theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.position="bottom")
ggsave(
  filename=hist_correlation_distribution_file,
  plot=hist_correlation_distribution,
  dpi = 1200
)
end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)

start_time <- Sys.time()
print("Create 3D plots")

#  Create cuts
x_c <- cut(networks_df$score, breaks= seq(-1, 1, by = 0.1))
y_c <- cut(networks_df$ss, breaks= seq(0, max_sample_cut, by = 100))

#  Calculate joint counts at cut levels
z <- table(x_c, y_c)

#  Plot as a 3D histogram
png(file=hist_3d_plot_file,
    width=1000, 
    height=1000,
    units = c("px"),
    pointsize=20
)
plot3D::hist3D(z=z, border="black", xlab = "Correlation", ylab = "Number of samples", zlab="Number of edges")
dev.off()

#  Plot as a 2D heatmap
png(file=heat_2d_plot_file,
    width=1000, 
    height=1000,
    units = c("px"),
    pointsize=20
)
plot3D::image2D(z=z, x=seq(-1, 1, 0.1), y=seq(0, max_sample_cut, 100), border="black", xlab = "Correlation", ylab = "Num. samples", clab="Num. edges")
dev.off()

end_time <- Sys.time()
time_diff = end_time - start_time
print(time_diff)
