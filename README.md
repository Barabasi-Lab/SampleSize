# Enhancing the accuracy of Network Medicine through understanding the impact of sample size in gene co-expression networks

## Description

This repository contains the scripts used to analyze the impact of sample size on gene co-expression networks. We used the following datasets:

* GTEx
* TCGA
* Rheumatoid Arthritis dataset

We created the gene co-expression networks using **Pearson correlation** and employed **bootstrapping** to generate gene co-expression networks of different sample sizes. We also calculated the **consensus networks** between the 5 replicate genes co-expression networks from each sample size and condition.

We made the following gene co-expression network analyses:

* **Link discovery rate analysis**: We analyzed the link discovery rate of gene expression datasets across gene co-expression networks of different sample sizes.
* **Power-law model fitting**: We fitted the gene co-expression networks of each dataset to a power-law model that describes the rate of discovery of statistically significant links across sample size.
* **Co-expression variation analysis**: We analyzed the variation in co-expression values across different network repetitions and sample sizes.
* **Differential co-expression analysis**: We calculated the differential co-expression between disease and normal gene co-expression networks.
* **Protein-protein interaction co-expression analysis**: We analyzed the influence of sample size on protein-protein interaction co-expression.

We also created summary tables and figures to illustrate the results.


## Table of contents

1. [Code](#code)
    1. [Extract the gene expression data](#1-extract-the-gene-expression-data)
    2. [Process the gene expression datasets](#2-process-the-gene-expression-datasets)
    3. [Bootstrapping](#3-bootstrapping)
    4. [Create the gene co-expression networks](#4-create-the-gene-co-expression-networks)
    5. [Analyze the gene co-expression networks](#5-analyze-the-gene-co-expression-networks)
    6. [Parse the results of the analyses and create summary tables](#6-parse-the-results)
    7. [Differential co-expression analysis](#7-differential-co-expression-analysis)
    8. [Analyze the gene co-expression variation in individual links](#8-analyze-the-gene-co-expression-variation-in-individual-links)
    9. [Analyze the gene co-expression of protein-protein interactions](#9-analyze-the-gene-co-expression-of-protein-protein-interactions)
    10. [Create the figures](#10-create-the-figures)


## Code

### 1. Extract the gene expression data

=> folder: `scripts/extract_data`

This is a detailed explanation on how to extract the datasets used in the analysis.

#### 1.1. GTEx

The data is downloaded from the GTEx portal (https://gtexportal.org/home/).

We downloaded the datasets from v8 release:

* **Gene expression reads**: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct
* **Sample information**: GTEx/v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
* **Subject information**: GTEx/v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

We executed the following notebook to analyze the data and prepare it for the analysis:

```
scripts/extract_data/GTEx_analysis_for_coexpr.Rmd
```

#### 1.3. TCGA

##### 1.3.1. Download the GDC-Client App

We downloaded the GDC-CLIENT app (Ubuntu - Client), necessary to download the dataset, at the following website: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool

This video was also useful to get to know how to download the data: https://www.youtube.com/watch?v=GDxj8DrkZok 

##### 1.3.2. Use the GenomicDataCommons R package to download the metadata and manifest

To download the TCGA dataset, we need to define the the parameters that we want to download and create a manifest file that will be used by the GDC-CLIENT app to download the dataset.

To analyze the data from TCGA, decide which parameters we wanted to use and create the manifest file, we used the GenomicDataCommons R package.

The whole process is detailed in the following Rmarkdown script: 

```
create_TCGA_manifest.Rmd
```

The script will generate the following files:
* Manifest file: `gdc_manifest.2022-11-18.txt`
* Metadata file: `metadata.txt`

##### 1.3.3. Download the dataset

To download the dataset, we executed the following command:

```
gdc-client download -m /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/raw/additional/gdc_manifest.2022-11-18.txt -d /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/raw/data
```

The process of downloading the dataset requires a lot of computational memory, so we used a computational cluster to run the command.

##### 1.3.4. Compile the dataset

To compile the scattered TCGA expression files into a unique file we use the following script:

```
python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/compile_tcga.py -d /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/raw/data -o /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/out/reads/TCGA_reads.csv -t unstranded
```

We can use a computational cluster to run the script, because it requires a lot of memory.

##### 1.3.5. Analyze and prepare the dataset

We run the following notebook to analyze and prepare the dataset:

```
TCGA_preparation.Rmd
```

### 2. Process the gene expression datasets

==> folder: process_gene_expression

Notebooks to run:

* For GTEx:

```
GTEx_preprocessing_for_coexpr.Rmd
```

* For GEO GSE193677:

```
GEO_preprocessing_for_coexpr_GSE193677.Rmd
```

* For Rheumatoid Arthritis dataset:

```
scipher_preprocessing.Rmd
```

* For TCGA:

```
TCGA_preprocessing_for_coexpr.Rmd
```

### 3. Bootstrapping

==> folder: bootstrapping

We create groups of samples of different sizes using random sampling with replacement, starting from size 20 and increasing the size by 20 until reaching the maximum number of samples in the dataset.

For this, we have to run a R markdown notebook for each dataset.

* For GTEx:

```
GTEx_subsampling.Rmd
```

* For GEO GSE193677

```
GEO_subsampling_GSE193677.Rmd
```

* For Rheumatoid Arthritis dataset:

```
scipher_subsampling.Rmd
```

* For TCGA:

```
TCGA_subsampling.Rmd
```


### 4. Create the gene co-expression networks

==> folder: gene_coexpression_networks

#### Individual co-expression networks

Execution:

```
Rscript create_gene_coexpression_network_from_samples_list.R -s <samples_file> -f <rnaseq_file> -o <output_file> -m <metric> -n <wto_n> -d <wto_delta> -p 6 -t signed -e pearson -a 0 -c bonferroni
```

This is done for all datasets.

To facilitate the calculations, we ran them in a computational cluster.

#### Consensus co-expression networks

Then, we calculated the consensus networks between the 5 replicate genes co-expression networks from each sample size and condition.

Execution:

```
Rscript create_consensus_gene_coexpression_network.R -l <list_coexpression_networks_file> -n <consensus_network_file> -t <threshold> -m <method>
```

This is done for all datasets.

To facilitate the calculations, we ran them in a computational cluster.

### 5. Analyze the gene co-expression networks

==> folder: gene_coexpression_networks

Execution:

```
Rscript analyze_coexpression_network_by_significant_edges.R -c <coexpression_network_file> -o <output_results_dir> -s <output_subgraphs_dir> -f <file_name> -t <threshold> -p <ppi_file> -d <disease_genes_file> -e <essential_genes_file> -g <genes_dataset_file>
```

This is done for all datasets.

To facilitate the calculations, we ran them in a computational cluster.


### 6. Parse the results of the analyses and create summary tables

==> folder: parse_results

1. Run: "parse_coexpression_analysis_results.Rmd"

2. Run: "server_example.R"

3. Run: "calculate_rnaseq_datasets_variation.Rmd"

4???? Is necessary?. Run: "calculate_convergence_correlation_types.Rmd"


### 7. Differential co-expression analysis

==> folder: differential_coexpression

#### Calculate differential co-expression analysis

Requires a lot of computational memory: 100000-120000 MB (150000 for a couple of them).

Execution:

```
calculate_differentially_coexpressed_genes.R -d <coexpression_network_file_D> -n <coexpression_network_file_N> -l <output_edges_file> -v <output_nodes_file> -t <threshold> -s <stretch_normalization> -f <filter_by_common_nodes>
```

#### Analyze the results from the differential co-expression analysis

```
Rscript analyze_differentially_coexpressed_genes.R -a <input_dir> -b <name_disease> -c <name_normal> -d <disease_gene_associations_file> -e <disease_name_in_associations_file> -f <ppi_file> -g <ppi_distances_file> -i <plots_dir> -j <tables_dir> -k <pval_threshold> -l <pval_correction> -m <nodes_to_follow_file> -n <drug_targets_file>
```

### 8. Analyze the gene co-expression variation in individual links

==> folder: analyze_coexpression_variation

#### Analysis of gene co-expression variation across five links

We illustrate the variation in correlation values of five randomly selected links with different correlation strengths from GTEx whole blood networks across five replicates of the same sample size (Figure 3B).

Execution:

```
Rscript analysis_specific_edges.R -n <networks_dir> -o <output_file> -m <method> -l <num_links_selected>
```

#### Analysis of standard deviation of links across repetitions of each sample size

First, we merge the consensus networks from different sample sizes into a unique file. This operation requires a lot of computational memory:

```
Rscript merge_networks.R -n <networks_dir> -m <method> -o <output_dataframe> -s <step> -r <max_rep>
```

Then, we read the file and plot the standard deviation of the repetitions across  sample size using a half violin plot (Figure 3C):

```
Rscript analysis_sd_coexpression_weight.R -n <networks_file> -t <output_sd_table_file> -p <output_sd_plot_file>
```

### 9. Analyze the gene co-expression of protein-protein interactions

==> folder: analyze_coexpression_ppi

We analyze the influence of sample size on protein-protein interaction co-expression.

The protein-protein interactions network was extracted from Gysi et al., PNAS (2023) (https://doi.org/10.1073/pnas.2301342120).

We have to run the following notebook:

```
compare_coexpression_to_ppi.Rmd
```

###  10. Create the figures

==> folder: general_analysis

We have to run the following notebook:

```
create_figures.Rmd
```
