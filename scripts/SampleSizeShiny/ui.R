library(markdown)

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Enable the use of shinyjs
  shinyjs::useShinyjs(),

  # Application title
  titlePanel("Sample size value in Network Medicine"),
  
  
  # Define navigation bar
  navbarPage("Categories",
             tabPanel("Topology",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_topology", "Dataset:",
                                       c('RA non-responder patients (Scipher)' = 'scipher_nonresponders',
                                         'RA responder patients (Scipher)' = 'scipher_responders',
                                         'GTEx' = 'gtex'
                                       )),
                          radioButtons("method", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       )),
                          radioButtons("parameter_topology", "Parameter:",
                                       c('Nodes' = 'nodes',
                                         'Edges' = 'edges',
                                         'Av. degree' = 'av_degree',
                                         'Av. path length' = 'av_path_length',
                                         'Av. clust. coef.' = 'av_clustering_coef',
                                         'Num. of components' = 'num_components',
                                         'Size of the LCC' = 'size_lcc',
                                         'Lost nodes' = 'lost_nodes',
                                         'Lost edges' = 'lost_edges',
                                         'Gained nodes' = 'gained_nodes',
                                         'Gained edges' = 'gained_edges')),
                          radioButtons("pval_topology", "P-value:",
                                       c('0.01' = 0.01,
                                         '0.05' = 0.05,
                                         'Group' = 'group_topology_pval')),
                          radioButtons("top_topology", "% of top scoring edges:",
                                           c('0.1%' = 0.1,
                                             '0.5%' = 0.5,
                                             '1%' = 1,
                                             '100%' = 100,
                                             'Group' = 'group_topology_top')),
                          verbatimTextOutput(outputId = "gtex_tissue"),
                          shinyWidgets::pickerInput(
                            inputId = "tissue",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'spleen',
                            multiple = TRUE),
                          verbatimTextOutput(outputId = "gtex_sex"),
                          shinyWidgets::pickerInput(
                            inputId = "sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'female',
                            multiple = TRUE)
                        ),
                        
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("topologyBoxPlot")
                        )
                      )
             ),
             
             
             tabPanel("Co-expressed PPIs",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_coexpressed_ppis", "Dataset:",
                                       c('GTEx' = 'gtex')),
                          radioButtons("method_coexpressed_ppis", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       )),
                          radioButtons("metric_coexpressed_ppis", "Evaluation metric:",
                                       c('Accuracy' = 'accuracy',
                                         'F1' = 'F1',
                                         'MCC' = 'MCC',
                                         'Correlation' = 'corr',
                                         'TPR' = 'TPR',
                                         'FPR' = 'FPR',
                                         'TP' = 'TP',
                                         'FP' = 'FP',
                                         'TN' = 'TN',
                                         'FN' = 'FN')),
                          radioButtons("parameter_coexpressed_ppis", "Parameter to predict high co-expression:",
                                       c('Short distance' = 'balanced',
                                         'Short distance & high clustering' = 'clustering',
                                         'Short distance & high degree' = 'degree',
                                         'Short distance & high betweenness' = 'betweenness')),
                          shinyWidgets::pickerInput(
                            inputId = "coexpressed_ppis_gtex_tissues",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'spleen',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "coexpressed_ppis_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "coexpressed_ppis_gtex_sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'female',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "coexpressed_ppis_gtex_sex")
                        ),
                          
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("coexpressedPPIsBoxPlot")
                        )
                      )
             ),


             tabPanel("Stability",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_stability", "Dataset:",
                                       c('GTEx' = 'gtex')),
                          radioButtons("method_stability", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       )),
                          radioButtons("type_interaction_stability", "Types of protein pairs (from the PPI network) to analyze:",
                                       c('All protein pairs' = 'all',
                                         'Only direct interactions' = 'direct',
                                         'Only pairs not interacting directly' = 'indirect')),
                          radioButtons("metric_stability", "Evaluation metric:",
                                       c('Difference' = 'difference.mean',
                                         'Number of edges' = 'number.of.edges',
                                         'Co-expression (network from subset of samples)' = 'score.subset.mean',
                                         'Co-expression (network from all samples)' = 'score.all.samples.mean',
                                         'Distance' = 'distance.mean',
                                         'Local clustering coeff.' = 'local_clustering_coefficient.mean',
                                         'Degree' = 'degree_centrality.mean',
                                         'Betweenness centrality' = 'betweenness_centrality.mean'
                                         )),
                          radioButtons("difference_range_stability", "Co-expression score difference:",
                                       c('> 0' = '0-Inf',
                                         '0 - 0.01' = '0-0.01',
                                         '0.01 - 0.05' = '0.01-0.05',
                                         '0.05 - 0.1' = '0.05-0.1',
                                         '> 0.1' = '0.1-Inf'
                                         )),
                          radioButtons("stability_group_by", "Group by:",
                                       c('Type of interaction' = 'group_type_interaction',
                                         'Range of co-expression score difference' = 'group_difference_range',
                                         'Not group' = 'no_group')),
                          shinyWidgets::pickerInput(
                            inputId = "stability_gtex_tissues",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'spleen',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "stability_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "stability_gtex_sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'female',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "stability_gtex_sex")
                        ),
                        
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("stabilityBoxPlot")
                        )
                      )
             ),

             navbarMenu("Co-expression scores",
                        
               tabPanel("By ranges",
                        
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(
                            radioButtons("dataset_scores", "Dataset:",
                                         c('GTEx' = 'gtex')),
                            radioButtons("method_scores", "Method:",
                                         c('wTO' = 'wto',
                                           'WGCNA' = 'wgcna',
                                           'ARACNE' = 'aracne',
                                           'Pearson' = 'pearson',
                                           'Spearman' = 'spearman'
                                         )),
                            radioButtons("type_interaction_scores", "Types of protein pairs (from the PPI network) to analyze:",
                                         c('All protein pairs' = 'all',
                                           'Only direct interactions' = 'direct',
                                           'Only pairs not interacting directly' = 'indirect',
                                           'Group by type of interaction' = 'group')),
                            radioButtons("metric_scores", "Evaluation metric:",
                                         c('Difference' = 'difference.mean',
                                           'Number of edges' = 'number.of.edges',
                                           'Co-expression (network from subset of samples)' = 'score.subset.mean',
                                           'Co-expression (network from all samples)' = 'score.all.samples.mean',
                                           'Distance' = 'distance.mean',
                                           'Local clustering coeff.' = 'local_clustering_coefficient.mean',
                                           'Degree' = 'degree_centrality.mean',
                                           'Betweenness centrality' = 'betweenness_centrality.mean'
                                         )),
                            checkboxGroupInput("score_range", label = "Co-expression difference:", 
                                               choices = list("All ranges grouped" = "group", 
                                                              "-1 to -0.9" = "-1--0.9", 
                                                              "-0.9 to -0.8" = "-0.9--0.8",
                                                              "0 to 0.1" = "0-0.1",
                                                              "0.1 to 0.2" = "0.1-0.2",
                                                              "0.2 to 0.3" = "0.2-0.3",
                                                              "0.3 to 0.4" = "0.3-0.4",
                                                              "0.4 to 0.5" = "0.4-0.5",
                                                              "0.5 to 0.6" = "0.5-0.6",
                                                              "0.6 to 0.7" = "0.6-0.7",
                                                              "0.7 to 0.8" = "0.7-0.8",
                                                              "0.8 to 0.9" = "0.8-0.9",
                                                              "0.9 to 1" = "0.9-1"
                                               ),
                                               selected = "group"),
                            shinyWidgets::pickerInput(
                              inputId = "scores_gtex_tissues",
                              label = "Tissue (if GTEx dataset)",
                              choices = c("spleen", "whole blood"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'spleen',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_gtex_tissues"),
                            shinyWidgets::pickerInput(
                              inputId = "scores_gtex_sex",
                              label = "Sex (if GTEx dataset)",
                              choices = c("female", "male"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'female',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_gtex_sex")
                          ),
                          
                          # Show boxplot of topological parameters
                          mainPanel(
                            plotOutput("scoresBoxPlot")
                          )
                        )
               ),
               
               tabPanel("By thresholds",
                        
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(
                            radioButtons("dataset_scores_thresholds", "Dataset:",
                                         c('GTEx' = 'gtex')),
                            radioButtons("method_scores_thresholds", "Method:",
                                         c('wTO' = 'wto',
                                           'WGCNA' = 'wgcna',
                                           'ARACNE' = 'aracne',
                                           'Pearson' = 'pearson',
                                           'Spearman' = 'spearman'
                                         )),
                            radioButtons("type_interaction_scores_thresholds", "Types of protein pairs (from the PPI network) to analyze:",
                                         c('All protein pairs' = 'all',
                                           'Only direct interactions' = 'direct',
                                           'Only pairs not interacting directly' = 'indirect')),
                            radioButtons("metric_scores_thresholds", "Evaluation metric:",
                                         c('Number of edges (network from subset of samples)' = 'number.of.edges.subset',
                                           'Number of edges (network from all samples)' = 'number.of.edges.all.samples',
                                           'Number of edges overlapped' = 'number.of.edges.overlapped',
                                           'Number of edges lost (from the network of all samples)' = 'number.of.edges.lost',
                                           'Number of edges gained (new in the subset network)' = 'number.of.edges.gained'
                                         )),
                            checkboxGroupInput("score_thresholds", label = "Score threshold:", 
                                               choices = list("All ranges grouped" = "group", 
                                                              "0.1" = "0.1",
                                                              "0.2" = "0.2",
                                                              "0.3" = "0.3",
                                                              "0.4" = "0.4",
                                                              "0.5" = "0.5",
                                                              "0.6" = "0.6",
                                                              "0.7" = "0.7",
                                                              "0.8" = "0.8",
                                                              "0.9" = "0.9"
                                               ),
                                               selected = "group"),
                            shinyWidgets::pickerInput(
                              inputId = "scores_thresholds_gtex_tissues",
                              label = "Tissue (if GTEx dataset)",
                              choices = c("spleen", "whole blood"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'spleen',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_thresholds_gtex_tissues"),
                            shinyWidgets::pickerInput(
                              inputId = "scores_thresholds_gtex_sex",
                              label = "Sex (if GTEx dataset)",
                              choices = c("female", "male"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'female',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_thresholds_gtex_sex")
                          ),
                          
                          # Show boxplot of topological parameters
                          mainPanel(
                            plotOutput("scoresThresholdsBoxPlot")
                          )
                        )
               )
             ),
             
             
             navbarMenu("Disease genes",
                        
               tabPanel("By score",
                        
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(
                            radioButtons("dataset_disease_genes_by_score", "Dataset:",
                                         c('GTEx' = 'gtex')),
                            radioButtons("method_disease_genes_by_score", "Method:",
                                         c('wTO' = 'wto',
                                           'WGCNA' = 'wgcna',
                                           'ARACNE' = 'aracne',
                                           'Pearson' = 'pearson',
                                           'Spearman' = 'spearman'
                                         )),
                            radioButtons("metric_disease_genes_by_score", "Evaluation metric:",
                                         c('Difference' = 'difference',
                                           'Number of edges' = 'num.edges',
                                           'Co-expression (network from subset of samples)' = 'score.subset.mean',
                                           'Co-expression (network from all samples)' = 'score.all.samples.mean',
                                           'Distance' = 'distances',
                                           'Local clustering coeff.' = 'local_clustering_coefficient.mean',
                                           'Degree' = 'degree_centrality.mean',
                                           'Betweenness centrality' = 'betweenness_centrality.mean'
                                         )),
                            radioButtons("disease_genes_by_score_group_by", "Group by:",
                                         c('Disease' = 'group_disease',
                                           'Disease class' = 'group_disease_class',
                                           'P-value' = 'group_disease_pval',
                                           '% of top scoring edges' = 'group_disease_top',
                                           'Not group' = 'no_group')),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_immune",
                              label = "Immune system diseases",
                              choices = c("immune system diseases", "arthritis, rheumatoid", "asthma", "graves disease", "multiple sclerosis"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'arthritis, rheumatoid',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_immune"),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_cardiovascular",
                              label = "Cardiovascular diseases",
                              choices = c("cardiovascular diseases", "arrhythmias, cardiac", "cardiomyopathies", "heart arrest", "myocardial ischemia"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'arrhythmias, cardiac',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_cardiovascular"),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_nervous",
                              label = "Nervous system diseases",
                              choices = c("nervous system diseases", "alzheimer disease", "epilepsy", "parkinson disease"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'alzheimer disease',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_nervous"),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_neoplasms",
                              label = "Neoplasms",
                              choices = c("neoplasms", "breast neoplasms", "colorectal neoplasm", "kidney neoplasms", "leukemia, myeloid", "lung neoplasms"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'breast neoplasms',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_neoplasms"),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_digestive",
                              label = "Digestive system diseases",
                              choices = c("digestive system diseases", "crohn disease", "cholestasis", "liver cirrhosis", "gastroenteritis"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              #selected = 'liver cirrhosis',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_digestive"),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_nutritional",
                              label = "Nutritional and metabolic diseases",
                              choices = c("nutritional and metabolic diseases", "amyotrophic lateral sclerosis", "celiac disease", "diabetes mellitus, type 1", "diabetes mellitus, type 2", "obesity"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              #selected = 'diabetes mellitus, type 2',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_nutritional"),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_endocrine",
                              label = "Endocrine system diseases",
                              choices = c("endocrine system diseases", "diabetes mellitus, type 2", "hyperthyroidism", "goiter"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              #selected = 'hyperthyroidism',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_endocrine"),
                            
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_gtex_tissues",
                              label = "Tissue (if GTEx dataset)",
                              choices = c("spleen", "whole blood"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'spleen',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_gtex_tissues"),
                            shinyWidgets::pickerInput(
                              inputId = "disease_genes_by_score_gtex_sex",
                              label = "Sex (if GTEx dataset)",
                              choices = c("female", "male"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'female',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "disease_genes_by_score_gtex_sex")
                          ),
                          
                          # Show boxplot of topological parameters
                          mainPanel(
                            plotOutput("diseaseByScoreBoxPlot")
                          )
                        )
               ),
               
               
               tabPanel("By thresholds",
  
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_disease_genes", "Dataset:",
                                       c('RA non-responder patients (Scipher)' = 'scipher_nonresponders',
                                         'RA responder patients (Scipher)' = 'scipher_responders')),
                          radioButtons("method_disease_genes", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       )),
                          radioButtons("parameter_disease_genes", "Parameter:",
                                       c('Num. disease genes in network' = 'disease_genes_in_network',
                                         '% disease genes in network' = 'percent_disease_genes_in_network',
                                         'Num. disease genes forming a LCC' = 'disease_genes_in_lcc',
                                         '% disease genes forming a LCC' = 'percent_disease_genes_in_lcc')),
                          radioButtons("pval_disease_genes", "P-value:",
                                           c('0.01' = 0.01,
                                             '0.05' = 0.05)),
                          radioButtons("top_disease_genes", "% of top scoring edges:",
                                       c('0.1%' = 0.1,
                                         '0.5%' = 0.5,
                                         '1%' = 1,
                                         '100%' = 100)),
                          radioButtons("group_by", "Group by:",
                                       c('Disease' = 'group_disease',
                                         'Disease class' = 'group_disease_class',
                                         'P-value' = 'group_disease_pval',
                                         '% of top scoring edges' = 'group_disease_top',
                                         'Not group' = 'no_group')),
                          shinyWidgets::pickerInput(
                            inputId = "immune",
                            label = "Immune system diseases",
                            choices = c("immune system diseases", "arthritis, rheumatoid", "asthma", "graves disease", "multiple sclerosis"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'arthritis, rheumatoid',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "diseases_immune"),
                          shinyWidgets::pickerInput(
                            inputId = "cardiovascular",
                            label = "Cardiovascular diseases",
                            choices = c("cardiovascular diseases", "arrhythmias, cardiac", "cardiomyopathies", "heart arrest", "myocardial ischemia"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'heart arrest',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "diseases_cardiovascular"),
                          shinyWidgets::pickerInput(
                            inputId = "nervous",
                            label = "Nervous system diseases",
                            choices = c("nervous system diseases", "alzheimer disease", "epilepsy", "parkinson disease"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'alzheimer disease',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "diseases_nervous"),
                          shinyWidgets::pickerInput(
                            inputId = "neoplasms",
                            label = "Neoplasms",
                            choices = c("neoplasms", "breast neoplasms", "colorectal neoplasm", "kidney neoplasms", "leukemia, myeloid", "lung neoplasms"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'breast neoplasms',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "diseases_neoplasms"),
                          shinyWidgets::pickerInput(
                            inputId = "digestive",
                            label = "Digestive system diseases",
                            choices = c("digestive system diseases", "crohn disease", "cholestasis", "liver cirrhosis", "gastroenteritis"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            #selected = 'liver cirrhosis',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "diseases_digestive"),
                          shinyWidgets::pickerInput(
                            inputId = "nutritional",
                            label = "Nutritional and metabolic diseases",
                            choices = c("nutritional and metabolic diseases", "amyotrophic lateral sclerosis", "celiac disease", "diabetes mellitus, type 1", "diabetes mellitus, type 2", "obesity"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            #selected = 'diabetes mellitus, type 2',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "diseases_nutritional"),
                          shinyWidgets::pickerInput(
                            inputId = "endocrine",
                            label = "Endocrine system diseases",
                            choices = c("endocrine system diseases", "diabetes mellitus, type 2", "hyperthyroidism", "goiter"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            #selected = 'hyperthyroidism',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "diseases_endocrine"),
                        ),
  
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("diseaseBoxPlot")
                        )
                      )
               )
            )
  )
  
)
