library(shiny)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridlayout)
library(bslib)
library(dplyr)
library(tidyr)
library(DT)
library(data.table)

#shiny::runApp('/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeShiny/app.R')
#shiny::runApp('/Users/j.aguirreplans/WORK/Postdoc/Projects/Scipher/SampleSize/scripts/SampleSizeShiny/app.R')


ui <- navbarPage(
  title = "MTCARS",
  selected = "Signature Visualization",
  collapsible = TRUE,
  theme = bslib::bs_theme(),
  tabPanel(
    title = "Signature Visualization",
    grid_container(
      layout = c(
        "sig_settings sig_viz"
      ),
      row_sizes = c(
        "1fr"
      ),
      col_sizes = c(
        "250px",
        "1fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "sig_settings",
        card_header("Settings"),
        card_body_fill(
          div(style = "overflow-y: auto;",
              verbatimTextOutput(outputId = "input_cars"),
              shinyWidgets::pickerInput(
                inputId = "input_cars",
                label = "Select cars:",
                choices = c("Mazda RX4", "Mazda RX4 Wag", "Datsun 710"),
                options = list(
                  `actions-box` = TRUE,
                  size = 10,
                  `selected-text-format` = "count > 1",
                  container = "body"
                ),
                selected = "Mazda RX4",
                multiple = TRUE
              ),
              shinyWidgets::materialSwitch(inputId = "input_filter",
                                           label = "Filter signatures:",
                                           status = "primary"),
              conditionalPanel(
                condition = "input.input_filter == true",
                style = "margin-bottom: 20px;", # Add margin-bottom style here
                verbatimTextOutput(outputId = "input_carb"),
                shinyWidgets::pickerInput(inputId = "input_carb",
                                          label = "Filter by carb:",
                                          choices = unique(mtcars$carb),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          selected = unique(mtcars$carb),
                                          multiple = TRUE),
                verbatimTextOutput(outputId = "input_gear"),
                shinyWidgets::pickerInput(inputId = "input_gear",
                                          label = "Filter by gear:",
                                          choices = unique(mtcars$gear),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          selected = unique(mtcars$gear),
                                          multiple = TRUE),
                verbatimTextOutput(outputId = "pert_itime"),
                shinyWidgets::pickerInput(inputId = "pert_itime",
                                          label = "Filter by exposure time:",
                                          choices = list(
                                            "6 h" = "6 h",
                                            "24 h" = "24 h"
                                          ),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          multiple = TRUE),
                sliderInput(inputId = "input_mpg_cutoff",
                            label = "Filter by mpg value:",
                            value = 20,
                            min = 10,
                            max = 34),
                sliderInput(inputId = "ctl_sim_cutoff",
                            label = "Filter by control similarity:",
                            value = 0.51,
                            min = 0,
                            max = 1),
                actionButton(inputId = "filter_action",
                             label = "Filter"),
              ),
              shinyWidgets::materialSwitch(inputId = "input_merge",
                                           label = "Merge signatures:",
                                           status = "primary"),
              conditionalPanel(
                condition = "input.input_merge == true",
                verbatimTextOutput(outputId = "merge_cond"),
                shinyWidgets::pickerInput(inputId = "merge_cond",
                                          label = "Merge by same:",
                                          choices = c("carb", "gear"),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          selected = c("carb", "gear"),
                                          multiple = TRUE),
                selectInput(inputId = "merge_method",
                            label = "Method:", 
                            choices = list("DGES" = "dges",
                                           "DGES & TAS" = "dges.tas",
                                           "Spearman" = "corr"),
                            selected = "dges"
                ),
                actionButton(inputId = "merge_action",
                             label = "Merge"),
              ),
          ),
        )
      ),
      grid_card(
        area = "sig_viz",
        card_body_fill(
          tabsetPanel(
            id = "sig_viz_type",
            tabPanel(
              title = "Information",
              DTOutput(outputId = "sigTable", width = "100%")
            ),
            tabPanel(
              title = "Plot",
              grid_container(
                layout = c(
                  "sig_heatmap_area"
                ),
                row_sizes = c(
                  "1fr"
                ),
                col_sizes = c(
                  "1fr"
                ),
                gap_size = "10px",
                grid_card_plot(
                  area = "sig_heatmap_area",
                  outputId = "sigHeatmap"
                )
              )
            )
          )
        )
      )
    )
  )
)


#-----------#
# Read data #
#-----------#

# Read numbers and sd data
data_dir <- "data"
numbers_file <- paste(data_dir, "dataset_numbers_complete_graph.txt", sep = "/")
numbers_df <- fread(numbers_file)
numbers_df$dataset <- tolower(numbers_df$dataset)
meta_file <- paste(data_dir, "metadata_summary.txt", sep = "/")
meta_df <- fread(meta_file)
sd_genes_file <- paste(data_dir, "variation_by_gene.txt", sep = "/") # from => calculate_rnaseq_datasets_variation.Rmd
sd_genes_df <- fread(sd_genes_file)

# Read unnormalized data
input_dir <- file.path(data_dir, "example_pearson_pval_0.05")
topology_results_file <- paste(input_dir, "topology_results_pearson_pval_0.05.txt", sep = "/")
results_selected_df <- fread(topology_results_file)
topology_results_by_size_file <- paste(input_dir, "topology_results_mean_pearson_pval_0.05.txt", sep = "/")
topology_results_selected_by_size_df <- fread(topology_results_by_size_file)
predictions_file <- paste(input_dir, "predictions_pearson_pval_0.05.txt", sep='/')
predicted_results_df <- fread(predictions_file)
analytical_results_file <- paste(input_dir, "analytical_model_results_pearson_pval_0.05.txt", sep = "/")
topology_results_selected_analytical_df <- fread(analytical_results_file)
analytical_summary_file <- paste(input_dir, 'analytical_model_summary_pearson_pval_0.05.txt', sep = "/")
analytical_model_summary_df <- fread(analytical_summary_file)
analytical_regression_results_file <- paste(input_dir, 'analytical_model_regression_results_pearson_pval_0.05.txt', sep = "/")
stretched_exponential_regression_df <- fread(analytical_regression_results_file)
theoretical_sample_size_file <- paste(input_dir, 'theoretical_sample_size_for_correlations_of_datasets.txt', sep = "/") # from => calculate_convergence_correlation_types.Rmd
sample_size_correlation_df <- fread(theoretical_sample_size_file)

# Read normalized data
topology_type_normalization <- "divide.L"
topology_results_file <- paste(input_dir, "/", "topology_results_norm_", topology_type_normalization, "_pearson_pval_0.05.txt", sep = "")
results_selected_norm_df <- fread(topology_results_file)
topology_results_by_size_file <- paste(input_dir, "/", 'topology_results_mean_norm_', topology_type_normalization, '_pearson_pval_0.05.txt', sep = "")
topology_results_selected_by_size_norm_df <- fread(topology_results_by_size_file)
predictions_file <- paste(input_dir, "/", 'predictions_', topology_type_normalization, '_pearson_pval_0.05.txt', sep = "")
predicted_results_norm_df <- fread(predictions_file)
analytical_results_file <- paste(input_dir, "/", 'analytical_model_results_', topology_type_normalization, '_pearson_pval_0.05.txt', sep = "")
topology_results_selected_analytical_norm_df <- fread(analytical_results_file)

# Create dictionary of dataset names
name_dict <- c(
  "scipher:scipher.sample.per.patient.baseline" = "R. arthritis",
  "gtex:whole.blood" = "GTEx: Whole blood",
  "gtex:lung" = "GTEx: Lung",
  "gtex:breast.mammary.tissue_female" = "GTEx: Breast",
  "tcga:tcga-luad" = "TCGA: Lung cancer",
  "tcga:tcga-brca_female" = "TCGA: Breast cancer"
)

# Create a dictionary name to color
name2color <- c(
  "TCGA: Lung cancer" = "#56B4E9",
  "tcga:tcga-luad" = "#56B4E9",
  "tcga:tcga-luad-tumor" = "#56B4E9", # #D55E00
  "TCGA: Breast cancer" = "#0072B2",
  "TCGA: Breast c." = "#0072B2",
  "tcga:tcga-brca_female" = "#0072B2",
  "tcga:tcga-brca-tumor_female" = "#0072B2",
  "breast neoplasms"="#0072B2", 
  "GTEx: Whole blood" = "#65F3BF",
  "gtex:whole.blood" = "#65F3BF",
  "GTEx: Lung" = "#24D157",
  "gtex:lung" = "#24D157",
  "GTEx: Breast" = "#00BA37", #44AA99
  "gtex:breast.mammary.tissue" = "#00BA37",
  "gtex:breast.mammary.tissue_female" = "#00BA37",
  "R. arthritis" = "#D55E00", #E69F00
  "scipher:scipher.sample.per.patient.baseline" = "#D55E00",
  "arthritis rheumatoid"="#D55E00", 
  "asthma"="#d9f365", 
  "thyroid neoplasms"="#0072B2",
  "tcga" = "#0072B2",
  "TCGA" = "#0072B2",
  "gtex" = "#00BA37",
  "GTEx" = "#00BA37",
  "scipher" = "#D55E00",
  #"Analytical model"="#e3f705",
  "Analytical model"="#cc68f7",
  "Power law"="#cc68f7",
  "Logarithm"="#09b863",
  "Exponential decay"="#09b863",
  "\u2265 0.2" = "#F8766D",
  "\u2265 0.4" = "#7CAE00",
  "\u2265 0.6" = "#00BFC4",
  "\u2265 0.8" = "#C77CFF",
  "Very Weak" = "#ff7b00",
  "Weak" = "#F8766D",
  "Moderate" = "#7CAE00",
  "Strong" = "#00BFC4",
  "Very Strong" = "#C77CFF",
  "Convergence point" = "red",
  "Model prediction" = "red",
  "P-value < 0.05" = "red",
  "common" = "#619CFF",
  "disease-gene" = "#cc68f7",
  "disease-specific" = "#F8766D",
  "normal-specific" = "#00BA38",
  "undefined" = "#808080",
  "different" = "#C77CFF",
  "all" = "#CD9600",
  "common (constant)" = "#619CFF",
  "common (change)" = "#b5d1ff",
  "disease-specific (constant)" = "#F8766D",
  "disease-specific (change)" = "#ffcbc7",
  "normal-specific (constant)" = "#00BA38",
  "normal-specific (change)" = "#a3d1b1",
  "undefined (constant)" = "#808080",
  "undefined (change)" = "#c2c2c2"
)

# Define selected parameters
datasets_selected <- unique(topology_results_selected_by_size_df$type_dataset)
model_selected <- "Stretched exponential (by linear fit)"

# Calculate data from scaling relationship
cols <- c("size", "S", "type_dataset", "S_lag1", "Fn", "slope", "intercept", "adj.r.squared")

scaling_relation_df <- data.frame(matrix(ncol = length(cols),
                                        nrow = 0,
                                        dimnames = list(NULL, cols)))

for (dataset_selected in datasets_selected){
  scaling_relation_filt <- topology_results_selected_by_size_df %>% 
    filter(type_dataset == dataset_selected) %>%
    select(size, mean, type_dataset) %>%
    unique() %>%
    rename("S"="mean") %>%
    arrange(size) %>%
    mutate(S_lag1 = if_else(size == min(size), 0, lag(S))) %>%
    mutate(Fn = ((S - S_lag1)/S_lag1)) %>%
    filter((!(is.infinite(Fn))) & (Fn > 0))
  lm_summary <- summary(lm(log(scaling_relation_filt$Fn) ~ log(scaling_relation_filt$size)))
  slope <- coef(lm_summary)[2]
  intercept <- coef(lm_summary)[1]
  adj.r.squared <- lm_summary$adj.r.squared
  scaling_relation_df <- rbind(scaling_relation_df, 
                              cbind(scaling_relation_filt, 
                                    data.frame(slope = slope, 
                                               intercept = intercept, 
                                               adj.r.squared = adj.r.squared)))
}

# Create scaling relation summary table
scaling_relation_summary_df <- scaling_relation_df %>% 
  dplyr::select(type_dataset, adj.r.squared) %>% 
  unique() %>% 
  arrange(factor(type_dataset, levels = datasets_selected)) %>%
  mutate(type_dataset = dplyr::recode(type_dataset, !!!name_dict)) %>%
  rename("R**2" = "adj.r.squared") %>%
  mutate_if(is.numeric, ~sprintf("%.2f",.))

# Show all information for model selected
analytical_model_summary_df <- analytical_model_summary_df %>% 
  inner_join((results_selected_norm_df %>%
                select("type_dataset", "subclassification")), 
             by = c("type_dataset")) %>%
  unique()
analytical_model_summary_df$L_vs_total <- abs((analytical_model_summary_df$unnorm_L - analytical_model_summary_df$total_num_edges)) / analytical_model_summary_df$total_num_edges
analytical_model_summary_df$dataset_name <- analytical_model_summary_df$type_dataset
analytical_model_summary_df <- 
  analytical_model_summary_df %>% 
  tidyr::separate(col = "dataset_name", 
           into = c("dataset_name", "sex"), 
           sep = "_", 
           remove = TRUE) %>%
  tidyr::separate(col = "dataset_name", 
           into = c("dataset", NULL), 
           sep = ":", 
           remove = FALSE) %>%
  mutate(sex = if_else(is.na(sex), "", sex))

analytical_model_summary_df$dataset_name = 
  ifelse(!(analytical_model_summary_df$subclassification == ""),
         paste(analytical_model_summary_df$dataset_name, 
               analytical_model_summary_df$subclassification, 
               sep = "-"), 
         analytical_model_summary_df$dataset_name)
analytical_model_summary_df$dataset_name = 
  ifelse(!(analytical_model_summary_df$sex == ""),
         paste(analytical_model_summary_df$dataset_name, 
               analytical_model_summary_df$sex, 
               sep = "_"), 
         analytical_model_summary_df$dataset_name)

analytical_model_summary_df %>% 
  dplyr::select(type_dataset, model, a, b, L, L_vs_total, adj.r.squared, relative.error.mean) %>% 
  filter((type_dataset %in% datasets_selected) & (model == model_selected)) %>% unique() %>%
  arrange(factor(type_dataset, levels = datasets_selected, datasets_selected))

# Show all information for logarithmic
log_model_summary_df <- analytical_model_summary_df %>% 
  dplyr::select(type_dataset, model, a, b, L, L_vs_total, adj.r.squared, relative.error.mean) %>% 
  filter((type_dataset %in% datasets_selected) & (model == "Logarithmic")) %>% unique() %>% 
  dplyr::select(type_dataset, adj.r.squared, relative.error.mean) %>% 
  rename("R**2" = "adj.r.squared", "epsilon" = "relative.error.mean") %>%
  arrange(factor(type_dataset, levels = datasets_selected)) %>%
  mutate_if(is.numeric, ~sprintf("%.2f",.))
print(log_model_summary_df)

# Create table for figure
power_law_summary_df <- analytical_model_summary_df %>% 
  dplyr::select(type_dataset, model, a, adj.r.squared, relative.error.mean) %>% 
  filter((type_dataset %in% datasets_selected) & (model == model_selected)) %>%
  dplyr::select(-model) %>%
  arrange(factor(type_dataset, levels = datasets_selected)) %>%
  mutate(type_dataset = dplyr::recode(type_dataset, !!!name_dict)) %>%
  rename("alpha"="a", "R**2" = "adj.r.squared", "epsilon" = "relative.error.mean") %>%
  unique() %>% 
  mutate_if(is.numeric, ~sprintf("%.2f",.))
print(power_law_summary_df)

# Bind power law with logarithm results
figure_summary_table <- cbind(scaling_relation_summary_df, (power_law_summary_df %>% select(-type_dataset)), (log_model_summary_df %>% select(-type_dataset)))



# Read mtcars data
mtcars_mod <- mtcars %>% 
  mutate(car = row.names(mtcars)) %>%
  dplyr::select(car, carb, gear, mpg)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  #freezeReactiveValue(input, "input_cars")
  shinyWidgets::updatePickerInput(session,
                                  inputId = "input_cars",
                                  label = "Select cars:",
                                  choices = unique(mtcars_mod$car),
                                  selected = unique(mtcars_mod$car),
                                  options = list(
                                    `actions-box` = TRUE,
                                    size = 10,
                                    `selected-text-format` = "count > 1"
                                  ))
  
  process_inputs <- reactive({
    mtcars_mod_filt <- mtcars_mod %>%
      filter(car %in% input$input_cars)
    if (input$input_filter) {
      input$filter_action
      isolate({
        mtcars_mod_filt <- mtcars_mod_filt %>%
          filter((carb %in% input$input_carb) &
                   (gear %in% input$input_gear) &
                   (mpg > input$input_mpg_cutoff)
          )
      })
    }
    if (input$input_merge) {
      input$merge_action
      isolate({
        mtcars_mod_filt <- mtcars_mod_filt %>%
          group_by_at(input$merge_cond) %>%
          summarize(mpg = median(mpg)) %>%
          unite(car, all_of(input$merge_cond), sep = "|")
      })
    }
    return(mtcars_mod_filt)
  })
  
  # Observe changes in input_cars and update input_gear and input_carb choices
  observeEvent(input$input_cars, {
    selected_cars <- input$input_cars
    filtered_data <- mtcars_mod %>% filter(car %in% selected_cars)
    gear_selected <- input$input_gear
    carb_selected <- input$input_carb
    if ((is.null(gear_selected)) | (!(isTRUE(input$input_filter)))) {
      gear_selected <- unique(filtered_data$gear)
    }
    if ((is.null(carb_selected)) | (!(isTRUE(input$input_filter)))) {
      carb_selected <- unique(filtered_data$carb)
    }
    # Update input_gear choices and retain previous selections
    shinyWidgets::updatePickerInput(
      session,
      inputId = "input_gear",
      label = "Filter by gear:",
      choices = unique(filtered_data$gear),
      selected = intersect(gear_selected, unique(filtered_data$gear)),
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = "count > 1"
      )
    )
    
    # Update input_carb choices and retain previous selections
    shinyWidgets::updatePickerInput(
      session,
      inputId = "input_carb",
      label = "Filter by carb:",
      choices = unique(filtered_data$carb),
      selected = intersect(carb_selected, unique(filtered_data$carb)),
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = "count > 1"
      )
    )
  })
  
  output$sigTable <- DT::renderDataTable({
    mtcars_mod_filt <- process_inputs()
    mtcars_mod_filt
  })

  output$sigHeatmap <- renderPlot({
    mtcars_mod_filt <- process_inputs()
    mtcars_mod_filt %>%
      ggplot(aes(x = mpg)) +
      geom_histogram(aes(fill = car))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
