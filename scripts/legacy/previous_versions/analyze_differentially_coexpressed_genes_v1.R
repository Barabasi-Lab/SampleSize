#!/usr/bin/env Rscript
#packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
library(igraph)
library(ggalluvial)
library(ggplot2)
require(magrittr)
library(tidyr)
set.seed(1510)
options(bitmapType='cairo')
`%ni%` <- Negate(`%in%`)


#---- Define inputs ----#

plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots"
tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables"
disease_genes_info_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022.csv"

# # For Breast Cancer (Female)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/reads/TCGA-BRCA_female___Breast.Mammary.Tissue_female"
# name_D = "tcga.brca.female"
# name_N = "gtex.breast.female"
# disease_name = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.fdr"
# pval_threshold = 0.05
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")

# # For Breast Cancer (both sex)
input_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/reads/TCGA-BRCA___Breast.Mammary.Tissue"
name_D = "tcga.brca"
name_N = "gtex.breast"
disease_name = "breast.neoplasms"
pval_correction_field = "pval_Phi_Tilde.adj.fdr"
pval_threshold = 0.05
nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")

# # For Breast Cancer - NORMAL TISSUE vs. GTEx Breast (both sex)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/reads/TCGA-BRCA-normal___Breast.Mammary.Tissue"
# name_D = "tcga.brca.normal"
# name_N = "gtex.breast"
# disease_name = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.fdr"
# pval_threshold = 0.05
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")

# # For Breast Cancer vs. Breast Cancer - NORMAL TISSUE (both sex)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/reads/TCGA-BRCA___TCGA-BRCA-normal"
# name_D = "tcga.brca"
# name_N = "tcga.brca.normal"
# disease_name = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.fdr"
# pval_threshold = 0.05
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")


# # For RA (complete.dataset)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/reads/scipher.complete.dataset___Whole.Blood"
# name_D = "scipher.complete.dataset"
# name_N = "gtex.whole.blood"
# disease_name = "arthritis.rheumatoid"
# pval_correction_field = "pval_Phi_Tilde.adj.fdr"
# pval_threshold = 0.05
# nodes_to_follow = c("TNF", "IL6ST", "IL6R", "PTPN22", "HOTAIR", "MALAT1")

# # For RA (complete.nonresponder)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/reads/scipher.complete.nonresponder___Whole.Blood"
# name_D = "scipher.complete.nonresponder"
# name_N = "gtex.whole.blood"
# disease_name = "arthritis.rheumatoid"
# pval_correction_field = "pval_Phi_Tilde.adj.fdr"
# pval_threshold = 0.05
# nodes_to_follow = c("TNF", "IL6ST", "IL6R", "PTPN22", "HOTAIR", "MALAT1")


#---- Read disease genes data ----#

GDA = fread(disease_genes_info_file)
GDA_disease = GDA %>% select("DiseaseName.no.sp.char", "HGNC_Symbol") %>% filter(DiseaseName.no.sp.char == disease_name)


#---- Compile files ----#

result_files = list.files(input_dir)
cols = c("file_name", "type_analysis", "size", "rep_D", "rep_N")
file_info_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))

for(file_name in result_files){
  info_file = gsub(".txt","",file_name)
  file_split1 = strsplit(info_file, split="_")[[1]]
  if(file_split1[1] == "diffanalysis"){
    type_analysis = file_split1[2]
    pval = tail(file_split1, n=1)
    info_file = gsub(paste("diffanalysis_", type_analysis, "_", sep=""),"",info_file)
    info_file = gsub(paste("_pval_", pval, sep=""),"",info_file)
    file_split2 = strsplit(info_file, split="___")[[1]]
    file_D_split = strsplit(gsub(".net","",file_split2[1]), split="_")[[1]]
    file_N_split = strsplit(gsub(".net","",file_split2[2]), split="_")[[1]]
    method = file_D_split[1]
    size = as.numeric(file_D_split[(length(file_D_split)-2)])
    r_D = as.numeric(file_D_split[(length(file_D_split))])
    r_N = as.numeric(file_N_split[(length(file_N_split))])
    dataset_D = file_D_split[(length(file_D_split)-4)]
    dataset_N = file_N_split[(length(file_N_split)-4)]
    file_info_df <- rbind(file_info_df, data.frame(file_name=file_name, type_analysis=type_analysis, size=size, rep_D=r_D, rep_N=r_N))
  }
}


#---- Read node files ----#

cols = c("Node", "DiseaseName.no.sp.char", "Phi_tilde", "pval", "size", "rep_D", "rep_N")
nodes_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
missing_files = c()
sizes_list = sort(unique(file_info_df$size))

for(x in seq(sizes_list)){
  size = sizes_list[x]
  #print(size)
  same_size_files_df = file_info_df %>% filter((size==!!size) & (type_analysis=="nodes")) %>% arrange(rep_D, rep_N)
  input_files = same_size_files_df$file_name
  # same_size_info_df = data.frame()
  for(input_file in input_files){
    r_D = (same_size_files_df %>% filter(file_name == input_file))$rep_D
    r_N = (same_size_files_df %>% filter(file_name == input_file))$rep_N
    individual_df = fread(paste(input_dir, input_file, sep="/"))
    if (nrow(individual_df) > 0){
      node_results_filt_df = individual_df %>%
        select(Node, Phi_tilde, !!as.symbol(pval_correction_field)) %>% 
        unique() %>% 
        rename("pval" = !!as.symbol(pval_correction_field)) %>% 
        left_join(GDA_disease, by=c("Node"="HGNC_Symbol"))
      if(nrow(nodes_df) == 0){
        nodes_df = cbind((node_results_filt_df %>% select("Node", "DiseaseName.no.sp.char", "Phi_tilde", "pval")), data.frame(size=size, rep_D=r_D, rep_N=r_N))
      } else {
        nodes_df = nodes_df %>% full_join(cbind((node_results_filt_df %>% select("Node", "DiseaseName.no.sp.char", "Phi_tilde", "pval")), data.frame(size=size, rep_D=r_D, rep_N=r_N)), by=c("Node", "DiseaseName.no.sp.char", "Phi_tilde", "pval", "size", "rep_D", "rep_N"))
      }
    } else {
      missing_files = c(missing_files, input_file)
    }
  }
}

# Define Phi_name
nodes_df$Phi_name = ifelse(nodes_df$Phi_tilde == "a", "common", ifelse(nodes_df$Phi_tilde == "b.D", "different", ifelse(nodes_df$Phi_tilde == "g.D", "disease-specific", ifelse(nodes_df$Phi_tilde == "g.N", "normal-specific", "undefined"))))


#---- Get undefined genes ----#

# We consider as "undefined" both the genes that are classified as U, but also the ones that are not significantly classified in any of the categories. 
# To do so, we need create an empty table with all considered genes for all sizes and repetitions, and fill it with the data that we already have
reps = sort(unique(paste(nodes_df$size, nodes_df$rep_D, nodes_df$rep_N, sep="-")))
nodes_empty_df = expand.grid(sort(unique(nodes_df$Node)), reps)
nodes_empty_df = nodes_empty_df %>% separate("Var2", into=c("size", "rep_D", "rep_N"), sep="-", convert=TRUE)
colnames(nodes_empty_df) = c("Node", "size", "rep_D", "rep_N")
nodes_empty_df = nodes_empty_df %>% left_join((nodes_df %>% select(Node, DiseaseName.no.sp.char) %>% unique()), by="Node") %>% arrange(Node, size, rep_D, rep_N)
nodes_expanded_df = nodes_df %>% full_join(nodes_empty_df, by = c("Node", "DiseaseName.no.sp.char", "size", "rep_D", "rep_N")) %>% arrange(Node, size, rep_D)
nodes_expanded_df = nodes_expanded_df %>%
  mutate(Phi_tilde = replace(Phi_tilde, is.na(Phi_tilde), "U")) %>%
  mutate(Phi_name = replace(Phi_name, is.na(Phi_name), "undefined"))

# Write output file
nodes_expanded_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_nodes", sep=""), sep="/")
nodes_expanded_df %>% fwrite(nodes_expanded_file)


#---- Count categories for each gene and size ----#

# Count for each gene and size the number of times that each category is predicted
nodes_counted_df = nodes_expanded_df %>% 
  group_by(Node, DiseaseName.no.sp.char, Phi_tilde, Phi_name, size) %>% 
  summarize(n_phi = n()) %>% 
  ungroup()  %>% 
  arrange(Node, size)

# Write output file
nodes_counted_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_nodes_counted.txt", sep=""), sep="/")
nodes_counted_df %>% fwrite(nodes_counted_file)


#---- Select predominant category for each gene and size ----#

# Select the for each gene and size the node type that occurs more times in the different repetitions
nodes_filtered_df = nodes_counted_df %>%
  group_by(Node, DiseaseName.no.sp.char, size) %>%
  filter(n_phi == max(n_phi)) %>% 
  arrange(Node, size) %>% 
  sample_n(1) %>% # Choose randomly cases in which there is the same number of node types 
  ungroup()

# Write output file
nodes_filtered_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_nodes_filtered.txt", sep=""), sep="/")
nodes_filtered_df %>% fwrite(nodes_filtered_file)
head(nodes_filtered_df)


#---- Count number of disease genes for each category and size ----#

# We count the number of disease genes in each category for each size
counts_categories_empty_df = expand.grid(sort(unique(nodes_expanded_df$Phi_name)), reps)
counts_categories_empty_df = counts_categories_empty_df %>% separate("Var2", into=c("size", "rep_D", "rep_N"), sep="-", convert=TRUE)
colnames(counts_categories_empty_df) = c("Phi_name", "size", "rep_D", "rep_N")
counts_categories_df = nodes_expanded_df %>% 
  filter(!(is.na(DiseaseName.no.sp.char))) %>% 
  group_by(Phi_name, size, rep_D, rep_N) %>% count() %>% rename("n_disease"="n") %>% arrange(size)
counts_categories_df = counts_categories_df %>% full_join(counts_categories_empty_df, by = c("Phi_name", "size", "rep_D", "rep_N")) %>% arrange(Phi_name, size, rep_D, rep_N)
counts_categories_df$n_disease = ifelse(is.na(counts_categories_df$n_disease), 0, counts_categories_df$n_disease)

# Write output file
counts_categories_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_category_counts.txt", sep=""), sep="/")
counts_categories_df %>% fwrite(counts_categories_file)
head(counts_categories_df)


#---- Plot number of disease genes for category and size ----#

# We plot the number of disease genes in each category

custom_palette = c(
  "#619CFF", #common
  "#F8766D", #disease-specific 
  "#00BA38", #normal-specific
  "#808080", #undefined
  "#C77CFF" #different
)

counts_categories_df %>% 
  group_by(Phi_name, size) %>%
  summarize(mean_disease = mean(n_disease), sd_disease = sd(n_disease), max_disease=max(n_disease), min_disease=min(n_disease)) %>%
  ggplot(aes(x = size, y = mean_disease, fill = Phi_name)) + 
  geom_line(aes(col = Phi_name)) +
  geom_ribbon(aes(ymin=min_disease, ymax=max_disease), alpha=0.2) +
  xlab("Number of samples") +
  ylab("Number of disease genes") +
  scale_fill_manual(values = custom_palette) +
  scale_color_manual(values = custom_palette) +
  guides(col=guide_legend(title="Gene category"), fill=guide_legend(title="Gene category")) +
  theme_linedraw() +
  theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 14, face="bold"), axis.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title=element_text(size=14, face="bold"))

plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_counts_disease_genes.png", sep=""), sep="/")
ggsave(
  plot_file,
  dpi = 1200,
  width = 9000,
  height = 6000,
  units = c("px")
)


#---- Plot flow between gene categories across sample size ----#

# We create the alluvium plot showing the flow between gene categories across sample size

custom_palette = c(
  "#619CFF", #common
  "#F8766D", #disease-specific 
  "#00BA38", #normal-specific
  "#808080", #undefined
  "#C77CFF" #different
)

if (max(sizes_list) <= 160) {
  sizes_to_plot = c(20, 40, 60, 80, 100, 120, 140, 160)
} else if ((max(sizes_list) > 160) & (max(sizes_list) <= 300)) {
  sizes_to_plot = c(20, 60, 100, 140, 180, 220, 260, 300)
} else if ((max(sizes_list) > 300) & (max(sizes_list) <= 440)) {
  sizes_to_plot = c(20, 80, 140, 200, 260, 320, 380, 440)
} else {
  sizes_to_plot = seq(20, max(sizes_list), 100)
}

nodes_filtered_df$Phi_name <- as.factor(nodes_filtered_df$Phi_name)
nodes_filtered_df$size = as.character(nodes_filtered_df$size)
levels(nodes_filtered_df$Phi_name) = sort(unique(nodes_filtered_df$Phi_name))
nodes_filtered_df$size <- factor(nodes_filtered_df$size, levels=as.character(sort(unique(as.integer(nodes_filtered_df$size)))))

nodes_filtered_df %>% 
  filter(size %in% sizes_to_plot) %>%
  filter(!(is.na(DiseaseName.no.sp.char))) %>%
  ggplot(aes(x = size, stratum = Phi_name, alluvium = Node, 
             fill = Phi_name, label = Phi_name)) +
  scale_fill_manual(values = custom_palette) +
  scale_color_manual(values = custom_palette) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum() +
  xlab("Number of samples") +
  ylab("Number of disease genes") +
  guides(fill=guide_legend(title="Gene category")) +
  theme_linedraw() +
  theme(plot.title =  element_text(size = 17, face="bold"), 
        axis.title = element_text(size = 14, face="bold"), 
        axis.text = element_text(size = 12), 
        #legend.position = "bottom",
        legend.text = element_text(size = 12), 
        legend.title=element_text(size=14, face="bold"))

plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_change_disease_genes.png", sep=""), sep="/")
ggsave(
  plot_file,
  dpi = 1200,
  width = 9000,
  height = 6000,
  units = c("px")
)


#---- Plot change of category of relevant genes ----#

for (important_node in nodes_to_follow){
  
  # Skip if node not in analysis
  if (!(important_node %in% unique(nodes_expanded_df$Node))){
    next
  }
  
  # Filter results by important node
  important_df = nodes_expanded_df %>% 
    filter(size %in% sizes_to_plot) %>%
    filter(!(is.na(DiseaseName.no.sp.char))) %>%
    filter(Node == important_node) %>%
    select(Phi_name, size, rep_D, rep_N) %>%
    group_by(Phi_name, size) %>% 
    summarize(n_evidences = n()) %>% 
    ungroup()
  
  # Define levels
  important_df$Phi_name <- as.factor(important_df$Phi_name)
  important_df$size = as.character(important_df$size)
  levels(important_df$Phi_name) = sort(unique(important_df$Phi_name))
  important_df$size <- factor(important_df$size, levels=as.character(sort(unique(as.integer(important_df$size)))))
  
  # Define palette
  if (levels(important_df$Phi_name) == c("common", "disease-specific", "normal-specific", "undefined")){
    personalized_palette = c("#619CFF", "#F8766D", "#00BA38", "#808080")
  } else if (levels(important_df$Phi_name) == c("disease-specific", "normal-specific", "undefined")){
    personalized_palette = c("#F8766D", "#00BA38", "#808080")
  } else if (levels(important_df$Phi_name) == c("common", "normal-specific", "undefined")){
    personalized_palette = c("#619CFF", "#00BA38", "#808080")
  } else if (levels(important_df$Phi_name) == c("common", "disease-specific", "undefined")){
    personalized_palette = c("#619CFF", "#F8766D", "#808080")
  } else if (levels(important_df$Phi_name) == c("normal-specific", "undefined")){
    personalized_palette = c("#00BA38", "#808080")
  } else if (levels(important_df$Phi_name) == c("disease-specific", "undefined")){
    personalized_palette = c("#F8766D", "#808080")
  } else if (levels(important_df$Phi_name) == c("normal-specific")){
    personalized_palette = c("#00BA38")
  } else if (levels(important_df$Phi_name) == c("disease-specific")){
    personalized_palette = c("#F8766D")
  } else {
    personalized_palette = c("#619CFF", "#F8766D", "#00BA38", "#808080", "#C77CFF")
  }
  
  # Plot
  important_node_plot = important_df %>%
    ggplot(aes(fill=Phi_name, y=n_evidences, x=size)) + 
    geom_bar(position="fill", stat="identity", width = 0.7, col="black") +
    #geom_bar(position = position_dodge(), stat="identity", col="black", width = 0.6) +
    xlab("Number of samples") +
    ylab("Fraction of evidences") +
    scale_fill_manual(values = personalized_palette) +
    guides(fill=guide_legend(title="Gene category")) +
    theme_linedraw() +
    theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 14, face="bold"), axis.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title=element_text(size=14, face="bold"))
  
  # Save plot
  plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_important_gene_", important_node, ".png", sep=""), sep="/")
  ggsave(
    plot_file,
    plot=important_node_plot,
    dpi = 1200,
    width = 9000,
    height = 6000,
    units = c("px")
  )
}


#---- Calculate fraction of genes that change or not change category for each sample size ----#

# Define output table
cols = c("Node", "Phi_name.prev", "size.prev", "Phi_name.curr", "size.curr")
change_of_category_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))

# Remove levels
nodes_filtered_df$size = as.numeric(as.character(nodes_filtered_df$size))

# For each size, compare nodes with the ones of previous size
for (x in seq(sizes_list)){
  size = sizes_list[x]
  nodes_filtered_by_selected_size = nodes_filtered_df %>%
    filter((!(is.na(DiseaseName.no.sp.char))) & (size == !!size)) %>%
    select(Node, Phi_name, size)
  if (x == 1){
    previous_nodes_filtered = data.frame(nodes_filtered_by_selected_size)
  } else {
    nodes_filtered_compared_by_sizes = previous_nodes_filtered %>%
      inner_join(nodes_filtered_by_selected_size, by=("Node")) %>%
      rename(c("Phi_name.prev"="Phi_name.x", "Phi_name.curr"="Phi_name.y", "size.prev"="size.x", "size.curr"="size.y"))
    change_of_category_df = rbind(change_of_category_df, nodes_filtered_compared_by_sizes)
    previous_nodes_filtered = data.frame(nodes_filtered_by_selected_size)
  }
}

# Detect changes
change_of_category_df$transition_type = ifelse(change_of_category_df$Phi_name.prev == change_of_category_df$Phi_name.curr, "stable", "change")
change_of_category_df$transition_name = paste(change_of_category_df$Phi_name.prev, change_of_category_df$Phi_name.curr, sep=" / ")

# Print table
change_of_category_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_changes.txt", sep=""), sep="/")
change_of_category_df %>% fwrite(change_of_category_file)


#---- Plot stability in gene categories ----#

# Calculate total of stable genes
total_stable_change_genes_df = change_of_category_df %>%
  group_by(size.prev) %>% 
  mutate(n_total = n()) %>%
  ungroup() %>%
  select(Node, Phi_name.prev, size.prev, transition_type, n_total) %>% 
  group_by(size.prev, transition_type, n_total) %>% 
  summarize(n_genes = n()) %>% 
  ungroup() %>%
  mutate(frac_genes = n_genes / n_total)
total_stable_change_genes_df$Phi_name.prev = "all"

# Calculate stable genes by category
categories_stable_change_genes_df = change_of_category_df %>%
  group_by(size.prev) %>% 
  mutate(n_total = n()) %>%
  ungroup() %>%
  filter(transition_type == "stable") %>%
  select(Node, Phi_name.prev, size.prev, transition_type, n_total) %>% 
  group_by(Phi_name.prev, size.prev, transition_type, n_total) %>% 
  summarize(n_genes = n()) %>% 
  ungroup() %>%
  mutate(frac_genes = n_genes / n_total)

# Plot
rbind(total_stable_change_genes_df, categories_stable_change_genes_df) %>%
  filter(transition_type == "stable") %>%
  ggplot(aes(x = size.prev, y = frac_genes, col = Phi_name.prev)) + 
  geom_line(aes(size = Phi_name.prev)) +
  xlab("Number of samples") +
  ylab("Fraction of stable genes") +
  scale_size_manual(values = c("all" = 1, "common" = 0.5, "disease-specific" = 0.5, "normal-specific" = 0.5, "undefined" = 0.5)) +
  scale_color_manual(values = c("#CD9600", "#619CFF", "#F8766D", "#00BA38", "#808080")) +
  guides(col=guide_legend(title="Gene category"), size=FALSE) +
  theme_linedraw() +
  theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 14, face="bold"), axis.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title=element_text(size=14, face="bold"))

# Save plot
plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_stable_genes.png", sep=""), sep="/")
ggsave(
  plot_file,
  dpi = 1200,
  width = 9000,
  height = 6000,
  units = c("px")
)


#---- Read edge files ----#

cols = c("Node.1", "Node.2", "DiseaseName.no.sp.char.1", "DiseaseName.no.sp.char.2", "Phi_tilde", "Score_Phi_tilde", "Score_ratio", "size", "rep_D", "rep_N")
edges_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
missing_edge_files = c()
sizes_list = sort(unique(file_info_df$size))

for(x in seq(sizes_list)){
  size = sizes_list[x]
  print(size)
  same_size_files_df = file_info_df %>% filter((size==!!size) & (type_analysis=="edges")) %>% arrange(rep_D, rep_N)
  input_files = same_size_files_df$file_name
  # same_size_info_df = data.frame()
  for(input_file in input_files){
    r_D = (same_size_files_df %>% filter(file_name == input_file))$rep_D
    r_N = (same_size_files_df %>% filter(file_name == input_file))$rep_N
    individual_df = fread(paste(input_dir, input_file, sep="/"))
    if (nrow(individual_df) > 0){
      edge_results_filt_df = individual_df %>%
        select(Node.1, Node.2, Phi_tilde, Score_Phi_tilde, Score_ratio) %>% 
        unique() %>% 
        left_join(GDA_disease, by=c("Node.1"="HGNC_Symbol")) %>% rename("DiseaseName.no.sp.char.1" = "DiseaseName.no.sp.char") %>%
        left_join(GDA_disease, by=c("Node.2"="HGNC_Symbol")) %>% rename("DiseaseName.no.sp.char.2" = "DiseaseName.no.sp.char")
      if(nrow(edges_df) == 0){
        edges_df = cbind((edge_results_filt_df %>% select("Node.1", "Node.2", "DiseaseName.no.sp.char.1", "DiseaseName.no.sp.char.2", "Phi_tilde", "Score_Phi_tilde", "Score_ratio")), data.frame(size=size, rep_D=r_D, rep_N=r_N))
      } else {
        edges_df = edges_df %>% full_join(cbind((edge_results_filt_df %>% select("Node.1", "Node.2", "DiseaseName.no.sp.char.1", "DiseaseName.no.sp.char.2", "Phi_tilde", "Score_Phi_tilde", "Score_ratio")), data.frame(size=size, rep_D=r_D, rep_N=r_N)), by=c("Node.1", "Node.2", "DiseaseName.no.sp.char.1", "DiseaseName.no.sp.char.2", "Phi_tilde", "Score_Phi_tilde", "Score_ratio", "size", "rep_D", "rep_N"))
      }
    } else {
      missing_edge_files = c(missing_edge_files, input_file)
    }
  }
}

