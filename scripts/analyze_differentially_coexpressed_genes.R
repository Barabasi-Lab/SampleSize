#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
source("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/GO.R")
source("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/Proximity_optimized.R")
library(optparse)
library(data.table)
library(dplyr)
library(igraph)
library(ggalluvial)
library(ggplot2)
require(magrittr)
library(patchwork)
library(tidyr)
set.seed(1510)
library(NetSci)
options(bitmapType='cairo')
`%ni%` <- Negate(`%in%`)


option_list = list(
  make_option(c("-a", "--input_dir"), action="store", type="character", 
              help="Directory with input files", metavar="character"),
  make_option(c("-b", "--name_disease"), action="store", type="character", 
              help="Disease name", metavar="character"),
  make_option(c("-c", "--name_normal"), action="store", type="character", 
              help="Normal name", metavar="character"),
  make_option(c("-d", "--disease_gene_associations_file"), action="store", type="character", 
              help="File containing disease-gene associations", metavar="character"),
  make_option(c("-e", "--disease_name_in_associations_file"), action="store", type="character", 
              help="Name of the disease in the disease-gene associations file", metavar="character"),
  make_option(c("-f", "--ppi_file"), action="store", type="character", 
              help="Protein-protein interactions network file", metavar="character"),
  make_option(c("-g", "--ppi_distances_file"), action="store", type="character", default = NA,
              help="Protein-protein interactions distances file", metavar="character"),
  make_option(c("-i", "--plots_dir"), action="store", type="character", 
              help="Directory to store the plots", metavar="character"),
  make_option(c("-j", "--tables_dir"), action="store", type="character",
              help="Directory to store the tables", metavar="character"),
  make_option(c("-k", "--pval_threshold"), action="store", type="double", default = 0.00001,
              help="P-value threshold", metavar="double"),
  make_option(c("-l", "--pval_correction"), action="store", type="character", default = "pval_Phi_Tilde.adj.bonf",
              help="P-value threshold", metavar="character"),
  make_option(c("-m", "--nodes_to_follow_file"), action="store", type="character",
              help="File with the nodes to follow", metavar="character"),
  make_option(c("-n", "--drug_targets_file"), action="store", type="character", 
              help="Drug-targets file used for the disease", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_differentially_coexpressed_network.R 


#---- Define inputs ----#

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$input_dir)){
  print_help(opt_parser)
  stop("Argument must be supplied (input_dir, nodes_to_follow_file).n", call.=FALSE)
}

# Read arguments
input_dir = opt$input_dir
name_D = opt$name_disease
name_N = opt$name_normal
disease_gene_associations_file = opt$disease_gene_associations_file
disease_name_in_associations_file = opt$disease_name_in_associations_file
ppi_file = opt$ppi_file
ppi_distances_file = opt$ppi_distances_file
plots_dir = opt$plots_dir
tables_dir = opt$tables_dir
pval_threshold = as.double(opt$pval_threshold)
pval_correction_field = opt$pval_correction
nodes_to_follow_file = opt$nodes_to_follow_file
drug_targets_file = opt$drug_targets_file

# Define other inputs
name2color = data.frame(name=c("common", "disease-specific", "normal-specific", "undefined", "different", "all"), color.codes=c("#619CFF", "#F8766D", "#00BA38", "#808080", "#C77CFF", "#CD9600"))


# PREDEFINED INPUTS
# disease_gene_associations_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info_2022.csv"
# plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression"
# tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression"
# ppi_file = "/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022.csv"
# ppi_distances_file = "/work/ccnr/j.aguirreplans/data/PPI/PPI_2022_04042022_distances_matrix.txt"

# # For Breast Cancer (female) vs. Breast (female)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-BRCA_female___TCGA-Breast_female/consensus"
# plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female"
# tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA_female___TCGA-Breast_female"
# name_D = "tcga.brca.female"
# name_N = "tcga.breast.female"
# disease_name_in_associations_file = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/nodes_to_follow/breast_cancer.txt"
# drug_targets_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/drug_targets/drugbank_targets_breast.neoplasm.txt"

# # For Breast Cancer (Female)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-BRCA_female___Breast.Mammary.Tissue_female"
# name_D = "tcga.brca.female"
# name_N = "gtex.breast.female"
# disease_name_in_associations_file = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")

# # For Breast Cancer (both sex)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-BRCA___Breast.Mammary.Tissue/consensus"
# name_D = "tcga.brca"
# name_N = "gtex.breast"
# disease_name_in_associations_file = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")

# # For Breast Cancer - NORMAL TISSUE vs. GTEx Breast (both sex)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-BRCA-normal___Breast.Mammary.Tissue/consensus"
# name_D = "tcga.brca.normal"
# name_N = "gtex.breast"
# disease_name_in_associations_file = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")

# # For Breast Cancer vs. Breast Cancer - NORMAL TISSUE (both sex)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-BRCA___TCGA-BRCA-normal/consensus"
# name_D = "tcga.brca"
# name_N = "tcga.brca.normal"
# disease_name_in_associations_file = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53")


# # For RA (complete.dataset)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/scipher.complete.dataset___Whole.Blood/consensus"
# name_D = "scipher.complete.dataset"
# name_N = "gtex.whole.blood"
# disease_name_in_associations_file = "arthritis.rheumatoid"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow = c("TNF", "IL6ST", "IL6R", "PTPN22", "HOTAIR", "MALAT1")

# # For RA (complete.nonresponder)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/scipher.complete.nonresponder___Whole.Blood/consensus"
# name_D = "scipher.complete.nonresponder"
# name_N = "gtex.whole.blood"
# disease_name_in_associations_file = "arthritis.rheumatoid"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow = c("TNF", "IL6ST", "IL6R", "PTPN22", "HOTAIR", "MALAT1")

# For KIRC vs. Kidney
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-KIRC___TCGA-Kidney/consensus"
# plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression/TCGA-KIRC___TCGA-Kidney"
# tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-KIRC___TCGA-Kidney"
# name_D = "tcga.kirc"
# name_N = "tcga.kidney"
# disease_name_in_associations_file = "kidney.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001

# For LUAD (Lung Adenocarcinoma) vs. Lung
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-LUAD___TCGA-Lung/consensus"
# plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression/TCGA-LUAD___TCGA-Lung"
# tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-LUAD___TCGA-Lung"
# name_D = "tcga.luad"
# name_N = "tcga.lung"
# disease_name_in_associations_file = "carcinoma.non.small.cell.lung"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001

# For LUSC (Lung Squamous Cell Carcinoma) vs. Lung
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-LUSC___TCGA-Lung/consensus"
# plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression/TCGA-LUSC___TCGA-Lung"
# tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-LUSC___TCGA-Lung"
# name_D = "tcga.lusc"
# name_N = "tcga.lung"
# disease_name_in_associations_file = "carcinoma.non.small.cell.lung"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001

# For LUAD (Lung Adenocarcinoma) vs. LUSC (Lung Squamous Cell Carcinoma)
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-LUAD___TCGA-LUSC/consensus"
# plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression/TCGA-LUAD___TCGA-LUSC"
# tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-LUAD___TCGA-LUSC"
# name_D = "tcga.luad"
# name_N = "tcga.lusc"
# disease_name_in_associations_file = "carcinoma.non.small.cell.lung"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001

# For BRCA.LumA vs. BRCA.LumB
# input_dir = "/scratch/j.aguirreplans/Scipher/SampleSizeProject/differential_coexpression_analysis/reads/TCGA-BRCA.LumA___TCGA-BRCA.LumB/consensus"
# plots_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots/plots_differential_coexpression/TCGA-BRCA.LumA___TCGA-BRCA.LumB"
# tables_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/tables/tables_differential_coexpression/TCGA-BRCA.LumA___TCGA-BRCA.LumB"
# name_D = "tcga.brca.luma"
# name_N = "tcga.brca.lumb"
# disease_name_in_associations_file = "breast.neoplasms"
# pval_correction_field = "pval_Phi_Tilde.adj.bonf"
# pval_threshold = 0.00001
# nodes_to_follow = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11", "TP53", "HER2", "ER", "ER1", "LIV1", "FOXA1", "KI67")


#---- Create directories if necessary ----#
if(!(is.na(plots_dir))){
  dir.create(plots_dir, showWarnings = FALSE)
}
if(!(is.na(tables_dir))){
  dir.create(tables_dir, showWarnings = FALSE)
}
tables_nodes_dir = paste(tables_dir, "nodes_by_type_and_category_and_size", sep="/") 
if(!(is.na(tables_nodes_dir))){
  dir.create(tables_nodes_dir, showWarnings = FALSE)
}


#---- Read disease genes data ----#

GDA = fread(disease_gene_associations_file)
GDA_disease = GDA %>% dplyr::select("DiseaseName.no.sp.char", "HGNC_Symbol") %>% dplyr::filter(DiseaseName.no.sp.char == disease_name_in_associations_file)
#unique(GDA$DiseaseName.no.sp.char[grepl("lung", GDA$DiseaseName.no.sp.char)])


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
    if(file_D_split[(length(file_D_split))] == "consensus"){
      size = as.numeric(file_D_split[(length(file_D_split)-1)])
      r_D = file_D_split[(length(file_D_split))]
      r_N = file_N_split[(length(file_N_split))]
    } else {
      size = as.numeric(file_D_split[(length(file_D_split)-2)])
      r_D = as.numeric(file_D_split[(length(file_D_split))])
      r_N = as.numeric(file_N_split[(length(file_N_split))])
    }
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
sizes_reps = sort(unique(paste(nodes_df$size, nodes_df$rep_D, nodes_df$rep_N, sep="-")))
nodes_empty_df = expand.grid(sort(unique(nodes_df$Node)), sizes_reps)
nodes_empty_df = nodes_empty_df %>% separate("Var2", into=c("size", "rep_D", "rep_N"), sep="-", convert=TRUE)
colnames(nodes_empty_df) = c("Node", "size", "rep_D", "rep_N")
nodes_empty_df = nodes_empty_df %>% left_join((nodes_df %>% select(Node, DiseaseName.no.sp.char) %>% unique()), by="Node") %>% arrange(Node, size, rep_D, rep_N)
# We join the info from the analysis with the empty table, and define as undefined the missing nodes in the analysis
nodes_expanded_df = nodes_df %>% full_join(nodes_empty_df, by = c("Node", "DiseaseName.no.sp.char", "size", "rep_D", "rep_N")) %>% arrange(Node, size, rep_D)
nodes_expanded_df = nodes_expanded_df %>%
  mutate(Phi_tilde = replace(Phi_tilde, is.na(Phi_tilde), "U")) %>%
  mutate(Phi_name = replace(Phi_name, is.na(Phi_name), "undefined"))
# We also consider as undefined nodes with a pvalue above the threshold
nodes_expanded_df = nodes_expanded_df %>%
  mutate(Phi_tilde = replace(Phi_tilde, pval >= pval_threshold, "U")) %>%
  mutate(Phi_name = replace(Phi_name, pval >= pval_threshold, "undefined"))

# Write output file
nodes_expanded_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_nodes.txt", sep=""), sep="/")
nodes_expanded_df %>% fwrite(nodes_expanded_file)

# Write background genes of the co-expression network
codina_coexpr_background_nodes_file = paste(tables_nodes_dir, "/codina_coexpr_background_nodes.txt", sep="")
codina_coexpr_background_nodes = unique(nodes_expanded_df$Node)
data.frame(Node=codina_coexpr_background_nodes) %>% fwrite(codina_coexpr_background_nodes_file)


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

reps = sort(unique(paste(nodes_df$rep_D, nodes_df$rep_N, sep="-")))
if(!(length(reps) == 1 && reps == "consensus-consensus")){
  # Select the for each gene and size the node type that occurs more times in the different repetitions
  nodes_filtered_df = nodes_counted_df %>%
    group_by(Node, DiseaseName.no.sp.char, size) %>%
    filter(n_phi == max(n_phi)) %>% 
    arrange(Node, size) %>% 
    sample_n(1) %>% # Choose randomly cases in which there is the same number of node types 
    ungroup()
} else {
  nodes_filtered_df = data.frame(nodes_counted_df)
}

# Write output file
nodes_filtered_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_nodes_filtered.txt", sep=""), sep="/")
nodes_filtered_df %>% fwrite(nodes_filtered_file)
#head(nodes_filtered_df)


#---- Count number of disease genes for each category and size ----#

# We count the number of disease genes in each category for each size
counts_categories_empty_df = expand.grid(sort(unique(nodes_expanded_df$Phi_name)), sizes_reps)
counts_categories_empty_df = counts_categories_empty_df %>% separate("Var2", into=c("size", "rep_D", "rep_N"), sep="-", convert=TRUE)
colnames(counts_categories_empty_df) = c("Phi_name", "size", "rep_D", "rep_N")
counts_categories_total_df = nodes_expanded_df %>% 
  group_by(Phi_name, size, rep_D, rep_N) %>% count() %>% rename("n_total"="n") %>% ungroup() %>% arrange(size)
counts_categories_disease_df = nodes_expanded_df %>% 
  filter(!(is.na(DiseaseName.no.sp.char))) %>% 
  group_by(Phi_name, size, rep_D, rep_N) %>% count() %>% rename("n_disease"="n") %>% ungroup() %>% arrange(size)
counts_categories_df = (counts_categories_total_df %>% inner_join(counts_categories_disease_df, by=c("Phi_name", "size", "rep_D", "rep_N"))) %>% full_join(counts_categories_empty_df, by = c("Phi_name", "size", "rep_D", "rep_N")) %>% arrange(Phi_name, size, rep_D, rep_N)
counts_categories_df$n_total = ifelse(is.na(counts_categories_df$n_total), 0, counts_categories_df$n_total)
counts_categories_df$n_disease = ifelse(is.na(counts_categories_df$n_disease), 0, counts_categories_df$n_disease)
#counts_categories_df$frac_disease = counts_categories_df$n_disease / counts_categories_df$n_total

# Write output file
counts_categories_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_category_counts.txt", sep=""), sep="/")
counts_categories_df %>% fwrite(counts_categories_file)
#head(counts_categories_df)


#---- Plot number of disease genes for category and size ----#

# We plot the number of disease genes in each category

# Get colors
counts_categories_col_df = counts_categories_df %>%
  left_join(name2color, by=c("Phi_name"="name"))

if(!(length(reps) == 1 && reps == "consensus-consensus")){
  plot_counts_categories = counts_categories_col_df %>% 
    group_by(Phi_name, size) %>%
    summarize(mean_disease = mean(n_disease), sd_disease = sd(n_disease), max_disease=max(n_disease), min_disease=min(n_disease)) %>%
    ggplot(aes(x = size, y = mean_disease, fill = Phi_name)) + 
    geom_line(aes(col = Phi_name)) +
    geom_ribbon(aes(ymin=min_disease, ymax=max_disease), alpha=0.2) +
    xlab("Number of samples") +
    ylab("Number of disease genes") +
    scale_fill_manual(values=setNames(counts_categories_col_df$color.codes, counts_categories_col_df$Phi_name)) +
    scale_color_manual(values=setNames(counts_categories_col_df$color.codes, counts_categories_col_df$Phi_name)) +
    guides(col=guide_legend(title="Gene category"), fill=guide_legend(title="Gene category")) +
    theme_linedraw() +
    theme(aspect.ratio=1,
          plot.title =  element_text(size = 15, face="bold"), 
          axis.title = element_text(size = 14, face="bold"), 
          axis.text = element_text(size = 13),
          panel.grid.major=element_blank(), 
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 13), 
          legend.title=element_text(size=14, face="bold"),
          text = element_text(family = "Helvetica"))
} else {
  plot_counts_categories = counts_categories_col_df %>% 
    group_by(Phi_name, size) %>%
    summarize(mean_disease = mean(n_disease)) %>%
    ggplot(aes(x = size, y = mean_disease)) + 
    geom_line(aes(col = Phi_name)) +
    xlab("Number of samples") +
    ylab("Number of disease genes") +
    scale_color_manual(values=setNames(counts_categories_col_df$color.codes, counts_categories_col_df$Phi_name)) +
    guides(col=guide_legend(title="Gene category"), fill=guide_legend(title="Gene category")) +
    theme_linedraw() +
    theme(aspect.ratio=1,
          plot.title =  element_text(size = 15, face="bold"), 
          axis.title = element_text(size = 14, face="bold"), 
          axis.text = element_text(size = 13),
          panel.grid.major=element_blank(), 
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 13), 
          legend.title=element_text(size=14, face="bold"),
          text = element_text(family = "Helvetica"))
}
plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_counts_disease_genes.png", sep=""), sep="/")
ggsave(
  plot_file,
  plot = plot_counts_categories,
  dpi = 1200,
  width = 9000,
  height = 6000,
  units = c("px")
)


#---- Plot flow between gene categories across sample size ----#

# We create the alluvium plot showing the flow between gene categories across sample size

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

# Get colors
nodes_filtered_col_df = nodes_filtered_df %>%
  left_join(name2color, by=c("Phi_name"="name"))

# Plot
plot_flow = nodes_filtered_col_df %>% 
  filter(size %in% sizes_to_plot) %>%
  filter(!(is.na(DiseaseName.no.sp.char))) %>%
  ggplot(aes(x = size, stratum = Phi_name, alluvium = Node, 
             fill = Phi_name, label = Phi_name)) +
  scale_fill_manual(values=setNames(nodes_filtered_col_df$color.codes, nodes_filtered_col_df$Phi_name)) +
  scale_color_manual(values=setNames(nodes_filtered_col_df$color.codes, nodes_filtered_col_df$Phi_name)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum() +
  xlab("Number of samples") +
  ylab("Number of disease genes") +
  guides(fill=guide_legend(title="Gene category")) +
  theme_linedraw() +
  theme(aspect.ratio=1,
        plot.title =  element_text(size = 15, face="bold"), 
        axis.title = element_text(size = 14, face="bold"), 
        axis.text = element_text(size = 13),
        panel.grid.major=element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 13), 
        legend.title=element_text(size=14, face="bold"),
        text = element_text(family = "Helvetica"))

plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_change_disease_genes.png", sep=""), sep="/")
ggsave(
  plot_file,
  plot = plot_flow,
  dpi = 1200,
  width = 9000,
  height = 6000,
  units = c("px")
)


#---- Plot change of category of relevant genes ----#

if(!(is.na(nodes_to_follow_file))){
  # Read nodes to follow
  nodes_to_follow = fread(nodes_to_follow_file, header = F)$V1
  if(!(length(reps) == 1 && reps == "consensus-consensus")){
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
      
      # Plot
      important_df = important_df %>%
        left_join(name2color, by=c("Phi_name"="name"))
      important_node_plot = important_df %>% 
        ggplot(aes(fill=Phi_name, y=n_evidences, x=size)) + 
        geom_bar(position="fill", stat="identity", width = 0.7, col="black") +
        #geom_bar(position = position_dodge(), stat="identity", col="black", width = 0.6) +
        xlab("Number of samples") +
        ylab("Fraction of evidences") +
        scale_fill_manual(values=setNames(important_df$color.codes, important_df$Phi_name)) +
        guides(fill=guide_legend(title="Gene category")) +
        theme_linedraw() +
        theme(aspect.ratio=1,
              plot.title =  element_text(size = 15, face="bold"), 
              axis.title = element_text(size = 14, face="bold"), 
              axis.text = element_text(size = 13),
              panel.grid.major=element_blank(), 
              panel.grid.minor = element_blank(),
              legend.text = element_text(size = 13), 
              legend.title=element_text(size=14, face="bold"),
              text = element_text(family = "Helvetica"))

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
  } else {
    for (important_node in nodes_to_follow){
      
      # Skip if node not in analysis
      if (!(important_node %in% unique(nodes_expanded_df$Node))){
        next
      }
      
      important_df = nodes_expanded_df %>% 
        filter(Node == important_node) %>%
        filter(!(is.na(DiseaseName.no.sp.char))) %>%
        filter(Node %in% nodes_to_follow) %>%
        select(Node, Phi_name, size)
      
      # Plot
      important_df = important_df %>%
        left_join(name2color, by=c("Phi_name"="name"))
      important_node_plot = important_df %>% 
        ggplot(aes(x=size, y=Phi_name, col=Phi_name, group=Node)) + 
        geom_line(size=1) +
        xlab("Number of samples") +
        ylab(NULL) +
        guides(col=guide_legend(title="Gene category")) +
        scale_color_manual(values=setNames(important_df$color.codes, important_df$Phi_name)) +
        theme_linedraw() +
        theme(aspect.ratio=1,
              plot.title =  element_text(size = 15, face="bold"), 
              axis.title = element_text(size = 14, face="bold"), 
              axis.text = element_text(size = 13),
              panel.grid.major=element_blank(), 
              panel.grid.minor = element_blank(),
              legend.text = element_text(size = 13), 
              legend.title=element_text(size=14, face="bold"),
              text = element_text(family = "Helvetica"))

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
  }
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


#---- Plot an improved representation of the flow between gene categories across sample size ----#

if (max(sizes_list) <= 160) {
  sizes_to_plot = c(20, 40, 60, 80, 100, 120, 140, 160)
} else if ((max(sizes_list) > 160) & (max(sizes_list) <= 300)) {
  sizes_to_plot = c(20, 60, 100, 140, 180, 220, 260, 300)
} else if ((max(sizes_list) > 300) & (max(sizes_list) <= 440)) {
  sizes_to_plot = c(20, 80, 140, 200, 260, 320, 380, 440)
} else {
  sizes_to_plot = seq(20, max(sizes_list), 100)
}

# Create last column of plot (size 100)
change_of_category_size100_df = change_of_category_df %>%
  filter(size.curr == 100) %>%
  dplyr::select(Node, Phi_name.curr, size.curr) %>%
  dplyr::rename("Phi_name.prev"="Phi_name.curr", "size.prev"="size.curr") %>%
  mutate(Phi_name.curr = NA) %>%
  mutate(size.curr = NA) %>%
  mutate(transition_type = NA) %>%
  mutate(transition_name = NA) %>%
  mutate(transition_name_simple = Phi_name.prev) %>%
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "common", "common (constant)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "disease-specific", "disease-specific (constant)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "normal-specific", "normal-specific (constant)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "undefined", "undefined (constant)")) %>% 
  left_join(name2color, by=c("Phi_name.prev"="name"))

# Count transitions / merge with colors / merge with last column
change_of_category_col_df = change_of_category_df %>%
  mutate(transition_name_simple = transition_name) %>%
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "common / common", "common (constant)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple %in% c("common / disease-specific", "common / normal-specific", "common / undefined"), "common (change)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "disease-specific / disease-specific", "disease-specific (constant)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple %in% c("disease-specific / common", "disease-specific / normal-specific", "disease-specific / undefined"), "disease-specific (change)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "normal-specific / normal-specific", "normal-specific (constant)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple %in% c("normal-specific / common", "normal-specific / disease-specific", "normal-specific / undefined"), "normal-specific (change)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple == "undefined / undefined", "undefined (constant)")) %>% 
  mutate(transition_name_simple = replace(transition_name_simple, transition_name_simple %in% c("undefined / disease-specific", "undefined / normal-specific", "undefined / common"), "undefined (change)")) %>% 
  left_join(name2color, by=c("transition_name_simple"="name")) %>%
  full_join(change_of_category_size100_df)

# Define factors
change_of_category_col_df$transition_name_simple <- factor(change_of_category_col_df$transition_name_simple, levels=c("common (constant)", "common (change)", "disease-specific (constant)",  "disease-specific (change)", "normal-specific (constant)", "normal-specific (change)", "undefined (constant)", "undefined (change)"))
change_of_category_col_df$size.prev = as.character(change_of_category_col_df$size.prev)
change_of_category_col_df$size.prev <- factor(change_of_category_col_df$size.prev, levels=as.character(sort(unique(as.integer(change_of_category_col_df$size.prev)))))

# Plot
plot_flow = change_of_category_col_df %>% 
  filter(size.prev %in% sizes_to_plot) %>%
  ggplot(aes(x = size.prev, stratum = transition_name_simple, alluvium = Node, 
             fill = transition_name_simple, label = transition_name_simple)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum() +
  scale_fill_manual(values=setNames(change_of_category_col_df$color.codes, change_of_category_col_df$transition_name_simple)) + # If I want to plot only specific elements in the legend, use the argument breaks with a list of the elements that I want to show
  scale_color_manual(values=setNames(change_of_category_col_df$color.codes, change_of_category_col_df$transition_name_simple)) +
  xlab("Number of samples") +
  ylab("Number of disease genes") +
  guides(fill=guide_legend(title="Gene category")) +
  theme_linedraw() +
  theme(aspect.ratio=1,
        plot.title =  element_text(size = 15, face="bold"), 
        axis.title = element_text(size = 14, face="bold"), 
        axis.text = element_text(size = 13),
        panel.grid.major=element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 13), 
        legend.title=element_text(size=14, face="bold"),
        text = element_text(family = "Helvetica"))

plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_change_disease_genes_improved.png", sep=""), sep="/")
ggsave(
  plot_file,
  plot = plot_flow,
  dpi = 1200,
  width = 10000,
  height = 6000,
  units = c("px")
)


#---- Plot stability in gene categories ----#

# Calculate total of stable genes
total_stable_change_genes_df = change_of_category_df %>%
  group_by(size.curr) %>% 
  mutate(n_total = n()) %>%
  ungroup() %>%
  select(Node, Phi_name.curr, size.curr, transition_type, n_total) %>% 
  group_by(size.curr, transition_type, n_total) %>% 
  summarize(n_genes = n()) %>% 
  ungroup() %>%
  mutate(frac_genes = n_genes / n_total)
total_stable_change_genes_df$Phi_name.curr = "all"

# Calculate stable genes by category
categories_stable_change_genes_df = change_of_category_df %>%
  group_by(size.curr) %>% 
  mutate(n_total = n()) %>%
  ungroup() %>%
  filter(transition_type == "stable") %>%
  select(Node, Phi_name.curr, size.curr, transition_type, n_total) %>% 
  group_by(Phi_name.curr, size.curr, transition_type, n_total) %>% 
  summarize(n_genes = n()) %>% 
  ungroup() %>%
  mutate(frac_genes = n_genes / n_total)

# Merge two tables
stable_genes_df = rbind(total_stable_change_genes_df, categories_stable_change_genes_df) %>%
  filter(transition_type == "stable") %>%
  left_join(name2color, by=c("Phi_name.curr"="name"))

# Plot
plot_stable_genes = stable_genes_df %>%
  ggplot(aes(x = size.curr, y = frac_genes, col = Phi_name.curr)) + 
  geom_line(aes(size = Phi_name.curr)) +
  xlab("Number of samples") +
  ylab("Fraction of stable genes") +
  scale_size_manual(values = c("all" = 1, "common" = 0.5, "disease-specific" = 0.5, "normal-specific" = 0.5, "undefined" = 0.5, "different" = 0.5)) +
  scale_color_manual(values=setNames(stable_genes_df$color.codes, stable_genes_df$Phi_name.curr)) +
  guides(col=guide_legend(title="Gene category"), size="none") +
  theme_linedraw() +
  theme(aspect.ratio=1,
        plot.title =  element_text(size = 15, face="bold"), 
        axis.title = element_text(size = 14, face="bold"), 
        axis.text = element_text(size = 13),
        panel.grid.major=element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 13), 
        legend.title=element_text(size=14, face="bold"),
        text = element_text(family = "Helvetica"))

# Save plot
plot_file = paste(plots_dir, paste("codina_", name_D, "_", name_N, "_stable_genes.png", sep=""), sep="/")
ggsave(
  plot_file,
  plot = plot_stable_genes,
  dpi = 1200,
  width = 9000,
  height = 6000,
  units = c("px")
)


#---- Plot stability + flow patchwork ----#

(plot_stable_genes /
  plot_flow) + 
  plot_annotation(tag_levels = 'A') & # Include letters 
  theme(plot.tag = element_text(face = 'bold')) # that are bold


plot_file = paste(plots_dir, "/codina_stability_results.png", sep="")
ggsave(
  plot_file,
  dpi = 1200,
  width = 9000,
  height = 12000,
  units = c("px")
)


#---- Read PPI file and drug targets file ----#

ppi_df = fread(ppi_file) %>% dplyr::select(HGNC_Symbol.1, HGNC_Symbol.2) %>% dplyr::rename("Node.1"="HGNC_Symbol.1", "Node.2"="HGNC_Symbol.2")
ppi_net = graph_from_data_frame(ppi_df, directed=F) %>% simplify()
if ((is.na(ppi_distances_file)) | (is.null(ppi_distances_file))){
  ppi_distances = distances(ppi_net) %>% as.data.frame()
} else {
  ppi_distances = fread(ppi_distances_file) %>% as.data.frame()
}

# Write background genes of the co-expression network
codina_ppi_background_nodes_file = paste(tables_nodes_dir, "/codina_ppi_background_nodes.txt", sep="")
codina_ppi_background_nodes = codina_coexpr_background_nodes[codina_coexpr_background_nodes %in% V(ppi_net)$name]
data.frame(Node=codina_ppi_background_nodes) %>% fwrite(codina_ppi_background_nodes_file)

# Read drug targets and filter by targets in the PPI network
drug_targets_df = fread(drug_targets_file) %>% 
  filter(Target %in% V(ppi_net)$name) %>% 
  select(ID, Target) %>% 
  unique() 
possible_targets = drug_targets_df$Target %>% unique()

#---- Calculate PPI module of the different categories of genes ----#

cols = c("size", "type_genes", "Node", "Phi_name")
codina_ppi_lcc_nodes_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))

cols = c("size", "type_genes", "Phi_name", "num_nodes", "num_edges", "num_components", "num_nodes_lcc", "num_edges_lcc", "lcc_significance_z", "lcc_significance_pval")
codina_ppi_analysis_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))

# cols = c("size", "type_genes", "Phi_name", "GO.ID", "Term", "Annotated", "Significant", "Expected", "Rank in classic", "classic", "KS", "weight")
# codina_ppi_go_enrichment_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))

# cols = c("size", "type_genes", "Proximity", "Targets", "Drug", "Disease", "p")
# codina_ppi_lcc_proximity_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))

#type_genes_list = c("all", "disease-gene")
type_genes_list = c("all")

category_list = sort(c("disease-gene", as.vector(unique(nodes_filtered_df$Phi_name))))

for(size_selected in sizes_list){
  
  for (type_genes_selected in type_genes_list){

    for (category_selected in category_list){
      
      if ((type_genes_selected == "all") & (category_selected == "disease-gene") & ("disease-gene" %in% type_genes_list)){
        next # Skip because it is the same as type_genes_selected == "disease-gene" and category_selected == "disease-gene"
      }
      
      # Filter by type of analysis
      nodes_filtered_type_analysis_df = nodes_filtered_df %>% filter(size == size_selected)
      if ((type_genes_selected == "disease-gene") | (category_selected == "disease-gene")){
        nodes_filtered_type_analysis_df = nodes_filtered_df %>% filter(!(is.na(DiseaseName.no.sp.char)))
      }
      
      # Select category genes
      if (!(category_selected == "disease-gene")){
        nodes_filtered_type_analysis_df = nodes_filtered_type_analysis_df %>% filter(Phi_name == category_selected)
      }
      
      # Save all genes in co-expression network from selected category in a file
      codina_coexpr_nodes_category_file = paste(tables_nodes_dir, "/codina_coexpr_nodes_type_", type_genes_selected, "_cat_", category_selected, "_size_", size_selected, ".txt", sep="")
      data.frame(Node=unique(nodes_filtered_type_analysis_df$Node)) %>% fwrite(codina_coexpr_nodes_category_file)
      
      # Select genes in PPI
      category_ppi_genes_in_analysis = unique(nodes_filtered_type_analysis_df$Node[nodes_filtered_type_analysis_df$Node %in% V(ppi_net)$name])
      
      if (length(category_ppi_genes_in_analysis) > 0){
        
        # Calculate subgraph of the category
        category_subgraph = induced.subgraph(ppi_net, vids=category_ppi_genes_in_analysis)
        category_subgraph_df = as.data.frame(get.edgelist(category_subgraph))
        
        # Calculate LCC and significance
        category_components = igraph::components(category_subgraph)
        category_lcc = igraph::induced.subgraph(category_subgraph, vids = V(category_subgraph)[category_components$membership == which.max(category_components$csize)] )
        if (length(V(category_lcc)$name) > 0){
          category_lcc_sig = LCC_Significance(N = 1000 
                                              , Targets = V(category_lcc)$name
                                              , G = ppi_net
                                              #, bins=1 #With bin=1, the degree distribution is not maintained and the hypothesis changes
          )
          codina_ppi_lcc_nodes_df = rbind(codina_ppi_lcc_nodes_df, data.frame(size=size_selected, type_genes=type_genes_selected, Node=V(category_lcc)$name, Phi_name=category_selected))
          codina_ppi_lcc_nodes_category_file = paste(tables_nodes_dir, "/codina_ppi_lcc_nodes_type_", type_genes_selected, "_cat_", category_selected, "_size_", size_selected, ".txt", sep="")
          data.frame(Node=unique(V(category_lcc)$name)) %>% fwrite(codina_ppi_lcc_nodes_category_file)
        } else {
          category_lcc_sig = list(Z=NA, emp_p=NA) # If empty network, leave result as NA
        }
        codina_ppi_analysis_df = rbind(codina_ppi_analysis_df, data.frame(size=size_selected, type_genes=type_genes_selected, Phi_name=category_selected, 
                                                                          num_nodes=gorder(category_subgraph), num_edges=gsize(category_subgraph),
                                                                          num_components=category_components$no, num_nodes_lcc=gorder(category_lcc), num_edges_lcc=gsize(category_lcc),
                                                                          lcc_significance_z=category_lcc_sig$Z, lcc_significance_pval=category_lcc_sig$emp_p))
        
        # Calculate GO enrichment of PPI LCC nodes of the selected category
        # GO_enrichment = GO(ID_type = "symbol", 
        #                    g = V(category_lcc)$name,
        #                    ONTO = "BP", 
        #                    bg = V(ppi_net)$name)
        # codina_ppi_go_enrichment_df = rbind(codina_ppi_go_enrichment_df, cbind(data.frame(size=size_selected, type_genes=type_genes_selected, Phi_name=category_selected), GO_enrichment$Sign))
        
        # Calculate proximity to drug targets
        # codina_proximity_df = avr_proximity_multiple_target_sets(
        #   set = unique(drug_targets_df$ID),
        #   G = ppi_net,
        #   ST = drug_targets_df,
        #   source = V(category_lcc)$name,
        #   N = 1000
        # )
        # codina_proximity_df = proximity_sign(N = 1000, 
        #                                    disease_genes = V(category_lcc)$name, # Disease genes
        #                                    Drugs = drug_targets_df, # Dataframe of drugs (ID) and targets (Target)
        #                                    distance = ppi_distances, # PPI network distances matrix
        #                                    disease_name = category_selected,  # Disease name
        #                                    possible_targets = possible_targets) # Targets
        # codina_proximity_df = codina_proximity_df %>% dplyr::rename("Phi_name" = "Disease")
        
        # codina_ppi_lcc_proximity_df = rbind(codina_ppi_lcc_proximity_df, cbind(data.frame(data.frame(size=size_selected, type_genes=type_genes_selected)), codina_proximity_df))

      } else {
        codina_ppi_analysis_df = rbind(codina_ppi_analysis_df, data.frame(size=size_selected, type_genes=type_genes_selected, Phi_name=category_selected, num_nodes=0, num_edges=0, num_components=0, num_nodes_lcc=0, num_edges_lcc=0, lcc_significance_z=NA, lcc_significance_pval=NA))
      }
    }
        
  }

}

codina_ppi_lcc_nodes_file = paste(tables_dir, "codina_ppi_lcc_nodes.txt", sep="/")
codina_ppi_lcc_nodes_df %>% fwrite(codina_ppi_lcc_nodes_file)
codina_ppi_analysis_file = paste(tables_dir, "codina_ppi_analysis.txt", sep="/")
codina_ppi_analysis_df %>% fwrite(codina_ppi_analysis_file)
# codina_ppi_go_enrichment_file = paste(tables_dir, "codina_ppi_go_enrichment.txt", sep="/")
# codina_ppi_go_enrichment_df %>% fwrite(codina_ppi_go_enrichment_file)
# codina_ppi_lcc_proximity_file = paste(tables_dir, "codina_ppi_lcc_proximity.txt", sep="/")
# codina_ppi_lcc_proximity_df %>% fwrite(codina_ppi_lcc_proximity_file)


#---- Calculate separation between PPI module of different categories of genes ----#

# cols = c("size", "type_genes", "Phi_name.1", "Phi_name.2", "Jaccard.Index", "Sab", "pvalue_lt")
# codina_ppi_lcc_separation_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
for(size_selected in sizes_list){
  
  for (type_genes_selected in type_genes_list){
    
    All = codina_ppi_lcc_nodes_df %>% filter((size == size_selected) & (type_genes == type_genes_selected)) %>% select(Phi_name, Node) %>% unique()
    
    # Write input file for separation
    codina_ppi_lcc_nodes_separation_category_file = paste(tables_nodes_dir, "/codina_ppi_lcc_nodes_separation_type_", type_genes_selected, "_size_", size_selected, ".txt", sep="")
    All %>% fwrite(codina_ppi_lcc_nodes_separation_category_file)
    
    # # Calculate Jaccard
    # Jac = Jaccard(All)
    # 
    # # Calculate separation
    # S = separation_Significance(G = ppi_net,
    #                            ST = All,
    #                            correct_by_target = FALSE,
    #                            Threads = 1)
    # 
    # # Merge results
    # S_result = Jac %>% full_join(S, by=c("Node.1"="x", "Node.2"="y")) %>% dplyr::rename("Phi_name.1"="Node.1", "Phi_name.2"="Node.2")
    # codina_ppi_lcc_separation_df = rbind(codina_ppi_lcc_separation_df, cbind(data.frame(size=size_selected, type_genes=type_genes_selected), S_result))

  }
  
}

# codina_ppi_lcc_separation_file = paste(tables_dir, "codina_ppi_lcc_separation.txt", sep="/")
# codina_ppi_lcc_separation_df %>% fwrite(codina_ppi_lcc_separation_file)


#---- Read edge files and calculate ranking of nodes based on their number of links in the predicted category ----#

cols = c("Node", "max", "common", "disease-specific", "normal-specific", "different", "DiseaseName.no.sp.char", "Phi_tilde", "pval", "size", "rep_D", "rep_N", "Phi_name", "rank")
ranking_nodes_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
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
    individual_edges_df = fread(paste(input_dir, input_file, sep="/"))
    individual_nodes_df = nodes_expanded_df %>% 
      filter((size == !!size) & (rep_D == r_D) & (rep_N == r_N))
    if (nrow(individual_edges_df) > 0){
      # Define Phi_name
      individual_edges_df$Phi_name = ifelse(individual_edges_df$Phi_tilde == "a", "common", ifelse(individual_edges_df$Phi_tilde == "b.D", "different", ifelse(individual_edges_df$Phi_tilde == "g.D", "disease-specific", ifelse(individual_edges_df$Phi_tilde == "g.N", "normal-specific", "undefined"))))

      # Calculate the ranking of nodes
      ranking_nodes_individual_df = individual_edges_df %>%
        dplyr::select(Node.1, Node.2, Phi_name) %>%
        tidyr::pivot_longer(-Phi_name) %>%
        dplyr::group_by(value, Phi_name) %>%
        dplyr::summarise(n = n()) %>% # Get the number of significant links in the network associated to each category for each gene
        dplyr::ungroup() %>%
        dplyr::group_by(value) %>%
        dplyr::mutate(max = max(n)) %>% # Create new variable max where we store the maximum number of links among the different categories associated to each gene
        tidyr::pivot_wider(names_from = Phi_name, # Create new columns associated to the number of links for each biomarker category
                    values_from = n, 
                    values_fill = 0) %>%
        dplyr::inner_join(., individual_nodes_df, by = c("value" = "Node")) %>% # Join the nodes table
        dplyr::rename("Node" = "value") %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Phi_tilde) %>%
        dplyr::mutate(rank = rank(-max, ties.method = "first")) %>% # Rank nodes by each biomarker category
        dplyr::ungroup() %>% 
        as.data.frame() %>%
        dplyr::arrange(rank)
      
      ranking_nodes_df = rbind(ranking_nodes_df, ranking_nodes_individual_df)
    } else {
      missing_edge_files = c(missing_edge_files, input_file)
    }
    rm(individual_edges_df)
    rm(ranking_nodes_individual_df)
  }
}

# Write output file
ranking_nodes_file = paste(tables_dir, paste("codina_", name_D, "_", name_N, "_nodes_ranked.txt", sep=""), sep="/")
ranking_nodes_df %>% fwrite(ranking_nodes_file)



