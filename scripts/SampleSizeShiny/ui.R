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
                                         'RA responder patients (Scipher)' = 'scipher_responders')),
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
                          ),
                        
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("topologyBoxPlot")
                        )
                      )
                      
                      
             ),
             tabPanel("Disease genes",

                    # Sidebar with a slider input for number of bins 
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("dataset_disease_genes", "Dataset:",
                                     c('RA non-responder patients (Scipher)' = 'scipher_nonresponders',
                                       'RA responder patients (Scipher)' = 'scipher_responders')),
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
