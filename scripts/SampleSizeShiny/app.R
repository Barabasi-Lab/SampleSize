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
  title = "SAMPLE SIZE",
  selected = "Power-law scaling relationship",
  collapsible = TRUE,
  theme = bslib::bs_theme(),
  tabPanel(
    title = "Power-law scaling relationship",
    grid_container(
      layout = c(
        "power_law_settings power_law_viz"
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
        area = "power_law_settings",
        card_header("Settings"),
        card_body_fill(
          div(
            style = "overflow-y: auto;",
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
            shinyWidgets::materialSwitch(
              inputId = "input_filter",
              label = "Filter signatures:",
              status = "primary"
            ),
            conditionalPanel(
              condition = "input.input_filter == true",
              style = "margin-bottom: 20px;", # Add margin-bottom style here
              verbatimTextOutput(outputId = "input_carb"),
              shinyWidgets::pickerInput(
                inputId = "input_carb",
                label = "Filter by carb:",
                choices = unique(mtcars$carb),
                options = list(
                  `actions-box` = TRUE,
                  size = 10,
                  `selected-text-format` = "count > 1"
                ),
                selected = unique(mtcars$carb),
                multiple = TRUE
              ),
              verbatimTextOutput(outputId = "input_gear"),
              shinyWidgets::pickerInput(
                inputId = "input_gear",
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


            checkboxGroupInput(
              "dataset_input",
              label = "Dataset",
              c("GTEx" = "gtex",
                "TCGA" = "tcga",
                "R. arthritis" = "scipher"
              ),
              selected = "gtex"
            ),
            conditionalPanel(
              condition = "input.dataset_input.includes('gtex')",
              verbatimTextOutput(outputId = "type_dataset_gtex"),
              shinyWidgets::pickerInput(
                inputId = "type_dataset_gtex",
                label = "GTEx tissue:",
                choices = list(
                  "Breast" = "gtex:breast.mammary.tissue",
                  "Breast (female)" = "gtex:breast.mammary.tissue_female",
                  "Lung" = "gtex:lung",
                  "Skin" = "gtex:skin.sun.exposed.lower.leg",
                  "Thyroid" = "gtex:thyroid",
                  "Whole blood" = "gtex:whole.blood"
                ),
                options = list(
                  `actions-box` = TRUE,
                  size = 10,
                  `selected-text-format` = "count > 1"
                ),
                selected = "gtex:breast.mammary.tissue_female",
                multiple = TRUE
              ),
              # verbatimTextOutput(outputId = "sex_gtex"),
              # shinyWidgets::pickerInput(
              #   inputId = "sex_gtex",
              #   label = "GTEx sex:",
              #   choices = c("female", "male", "both"),
              #   options = list(
              #     `actions-box` = TRUE,
              #     size = 10,
              #     `selected-text-format` = "count > 1"
              #   ),
              #   selected = "both",
              #   multiple = TRUE
              # ),
            ),
            conditionalPanel(
              condition = "input.dataset_input.includes('tcga')",
              verbatimTextOutput(outputId = "type_dataset_tcga"),
              shinyWidgets::pickerInput(
                inputId = "type_dataset_tcga",
                label = "TCGA subset:",
                choices = list(
                  "BRCA (Breast cancer)" = "tcga:tcga-brca_female",
                  "BRCA.LumA (Breast Luminal A cancer)" = "tcga:brca.luma",
                  "BRCA.LumB (Breast Luminal B cancer)" = "tcga:brca.lumb",
                  "KIRC (Kidney renal clear cell cancer)" = "tcga:tcga-kirc",
                  "KIRP (Kidney renal papillary cell cancer)" = "tcga:tcga-kirp",
                  "LUAD (Lung adenocarcinoma)" = "tcga:tcga-luad",
                  "LUSC (Lung squamous cell carcinoma)" = "tcga:tcga-lusc",
                  "THCA (Thyroid carcinoma)" = "tcga:tcga-thca",
                  "Breast tissue" = "tcga:breast_female",
                  "Kidney tissue" = "tcga:kidney",
                  "Lung cancer" = "tcga:tcga-luad"
                ),
                options = list(
                  `actions-box` = TRUE,
                  size = 10,
                  `selected-text-format` = "count > 1"
                ),
                selected = "tcga:tcga-brca_female",
                multiple = TRUE
              ),
              # verbatimTextOutput(outputId = "type_tissue_tcga"),
              # shinyWidgets::pickerInput(
              #   inputId = "type_tissue_tcga",
              #   label = "TCGA tissue type:",
              #   choices = list("Tumor" = "tumor", "Tissue" = "normal"),
              #   options = list(
              #     `actions-box` = TRUE,
              #     size = 10,
              #     `selected-text-format` = "count > 1"
              #   ),
              #   selected = "tumor",
              #   multiple = TRUE),
            ),
            conditionalPanel(
              condition = "input.dataset_input.includes('scipher')",
              verbatimTextOutput(outputId = "type_dataset_scipher"),
              shinyWidgets::pickerInput(
                inputId = "type_dataset_scipher",
                label = "R. Arthritis subset:",
                choices = list(
                  "1 sample/patient at baseline" = "scipher:scipher.sample.per.patient.baseline",
                  "Complete dataset" = "scipher:scipher.complete.dataset"
                ),
                options = list(
                  `actions-box` = TRUE,
                  size = 10,
                  `selected-text-format` = "count > 1"
                ),
                selected = "scipher:scipher.sample.per.patient.baseline",
                multiple = TRUE
              ),
            ),
            # Add select input for type of plot
            selectInput(
              inputId = "plot_type",
              label = "Plot type:",
              choices = c(
                "Corr. vs. size",
                "Scaling rel.",
                "Model prediction"
              ),
              selected = "Corr. vs. size"
            )


          ),
        )
      ),
      grid_card(
        area = "power_law_viz",
        card_body_fill(
          tabsetPanel(
            id = "power_law_tabset",
            tabPanel(
              title = "Summary table",
              DTOutput(outputId = "summaryTable", width = "100%")
            ),
            tabPanel(
              title = "Corr. vs  size plots",
              grid_container(
                layout = c(
                  "corr_vs_size_area"
                ),
                row_sizes = c(
                  "1fr"
                ),
                col_sizes = c(
                  "1fr"
                ),
                gap_size = "10px",
                grid_card_plot(
                  area = "corr_vs_size_area",
                  outputId = "scalingPlot"
                )
              )
            ),
            tabPanel(
              title = "Discovery rate plots",
              grid_container(
                layout = c(
                  "discovery_rate_area"
                ),
                row_sizes = c(
                  "1fr"
                ),
                col_sizes = c(
                  "1fr"
                ),
                gap_size = "10px",
                grid_card_plot(
                  area = "discovery_rate_area",
                  outputId = "discoveryRatePlot"
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
dataset_to_type_df <- results_selected_df %>% dplyr::select(dataset, type_dataset) %>% unique()
gtex_datasets <- results_selected_df %>% filter(dataset == "gtex") %>% pull(type_dataset) %>% unique() # Get GTEx datasets
tcga_datasets <- results_selected_df %>% filter(dataset == "tcga") %>% pull(type_dataset) %>% unique() # Get TCGA datasets
scipher_datasets <- results_selected_df %>% filter(dataset == "scipher") %>% pull(type_dataset) %>% unique() # Get scipher datasets
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
  "gtex" = "GTEx",
  "tcga" = "TCGA",
  "scipher" = "R. arthritis",
  "gtex:whole.blood" = "GTEx: Whole blood",
  "gtex:lung" = "GTEx: Lung",
  "gtex:breast.mammary.tissue" = "GTEx: Breast",
  "gtex:breast.mammary.tissue_female" = "GTEx: Breast (Female)",
  "tcga:tcga-brca.luma" = "TCGA: Breast cancer (Luminal A)",
  "tcga:tcga-brca.lumb" = "TCGA: Breast cancer (Luminal B)",
  "tcga:breast_female" = "TCGA: Breast tissue (Female)",
  "tcga:kidney" = "TCGA: Kidney tissue",
  "tcga:tcga-luad" = "TCGA: Lung cancer",
  "tcga:tcga-brca_female" = "TCGA: Breast cancer (Female)",
  "scipher:scipher.complete.dataset" = "R. arthritis (all samples)",
  "scipher:scipher.sample.per.patient.baseline" = "R. arthritis (sample/patient)"
)

# Create a dictionary name to color
color_dict <- c(
  "TCGA: Lung cancer" = "#56B4E9",
  "tcga:tcga-luad" = "#56B4E9",
  "TCGA: Breast cancer" = "#0072B2",
  "TCGA: Breast cancer (Female)" = "#0072B2",
  "tcga:tcga-brca_female" = "#0072B2",
  "breast neoplasms" = "#0072B2",
  "tcga:tcga-brca.luma" = "#00b2a0",
  "TCGA: Breast cancer (Luminal A)" = "#00b2a0",
  "tcga:tcga-brca.lumb" = "#74f6ff",
  "TCGA: Breast cancer (Luminal B)" = "#74f6ff",
  "TCGA: Breast tissue" = "#0072B2",
  "TCGA: Breast tissue (Female)" = "#0072B2",
  "TCGA: Kidney tissue" = "#009E73",
  "tcga:kidney" = "#009E73",
  "GTEx: Whole blood" = "#65F3BF",
  "gtex:whole.blood" = "#65F3BF",
  "GTEx: Lung" = "#24D157",
  "gtex:lung" = "#24D157",
  "GTEx: Breast" = "#00BA37",
  "GTEx: Breast (Female)" = "#00BA37",
  "gtex:breast.mammary.tissue" = "#00BA37",
  "gtex:breast.mammary.tissue_female" = "#00BA37",
  "R. arthritis" = "#D55E00", #E69F00
  "R. arthritis (all samples)" = "#D55E00", #E69F00
  "R. arthritis (sample/patient)" = "#D55E00", #E69F00
  "scipher:scipher.sample.per.patient.baseline" = "#D55E00",
  "arthritis rheumatoid" = "#D55E00",
  "asthma" = "#d9f365",
  "thyroid neoplasms" = "#0072B2",
  "tcga" = "#0072B2",
  "TCGA" = "#0072B2",
  "gtex" = "#00BA37",
  "GTEx" = "#00BA37",
  "scipher" = "#D55E00",
  "Analytical model" = "#cc68f7",
  "Power law" = "#cc68f7",
  "Logarithm" = "#09b863",
  "Exponential decay" = "#09b863",
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

analytical_model_summary_df$dataset_name <-
  ifelse(!(analytical_model_summary_df$subclassification == ""),
         paste(analytical_model_summary_df$dataset_name,
               analytical_model_summary_df$subclassification,
               sep = "-"),
         analytical_model_summary_df$dataset_name)
analytical_model_summary_df$dataset_name <-
  ifelse(!(analytical_model_summary_df$sex == ""),
         paste(analytical_model_summary_df$dataset_name,
               analytical_model_summary_df$sex,
               sep = "_"), 
         analytical_model_summary_df$dataset_name)

# Show all information for logarithmic
log_model_summary_df <- analytical_model_summary_df %>%
  dplyr::select(dataset, type_dataset, model, a, b, L, L_vs_total, adj.r.squared, relative.error.mean) %>% 
  filter((type_dataset %in% datasets_selected) & (model == "Logarithmic")) %>% unique() %>% 
  dplyr::select(type_dataset, adj.r.squared, relative.error.mean) %>% 
  rename("R**2 (Log)" = "adj.r.squared", "epsilon (Log)" = "relative.error.mean") %>%
  arrange(factor(type_dataset, levels = datasets_selected)) %>%
  mutate_if(is.numeric, ~sprintf("%.2f",.))

# Create table for figure
power_law_summary_df <- analytical_model_summary_df %>%
  dplyr::select(dataset, type_dataset, model, a, adj.r.squared, relative.error.mean) %>% 
  filter((type_dataset %in% datasets_selected) & (model == model_selected)) %>%
  dplyr::select(-model) %>%
  arrange(factor(type_dataset, levels = datasets_selected)) %>%
  rename("alpha (model)"="a", "R**2 (model)" = "adj.r.squared", "epsilon (model)" = "relative.error.mean") %>%
  unique() %>% 
  mutate_if(is.numeric, ~sprintf("%.2f",.))

# Join model results with the mean of the empirical results
predicted_results_mean_df <-
  rbind((topology_results_selected_by_size_df %>%
           dplyr::select("type_dataset", "size", "mean") %>%
           unique() %>%
           rename("model_result" = "mean") %>%
           mutate(model = "Mean")),
        (predicted_results_df %>%
           dplyr::select(type_dataset, size, model_result, model) %>%
           unique()
         %>% mutate_at(c("model_result"), as.numeric))) %>%
  filter(type_dataset %in% datasets_selected) %>%
  filter(model == model_selected) %>%
  inner_join((analytical_model_summary_df %>% 
                dplyr::select("type_dataset", "model", "adj.r.squared")),
             by = c("type_dataset", "model")) %>%
  mutate_at(c('adj.r.squared'), ~sprintf("%.2f",.)) %>%
  mutate_at(c('adj.r.squared'), as.character)

# Join predicted results with empirical results
goodnessfit_df <- results_selected_df %>% 
  dplyr::select(dataset, type_dataset, size, rep, num_edges) %>%
  right_join(predicted_results_mean_df, by = c("type_dataset", "size")) %>%
  filter(type_dataset %in% datasets_selected) %>%
  filter(
    size <= max(
      (results_selected_df %>% filter(type_dataset %in% datasets_selected))$size
    )
  ) %>%
  filter(
    size >= min(
      (results_selected_df %>% filter(type_dataset %in% datasets_selected))$size
    )
  ) %>%
  filter(!(is.na(num_edges)))

# Join predicted results with empirical results
model_prediction_df <- results_selected_norm_df %>%
  dplyr::select("type_dataset", "size", "rep", "num_edges") %>%
  # Include predictions from model
  dplyr::right_join(
    (predicted_results_norm_df %>%
      dplyr::select(type_dataset, model, model_result, size) %>%
      mutate_at(c('model_result'), as.numeric) %>%
      unique()
    ), by=c("type_dataset", "size")) %>%
  dplyr::filter(
    (type_dataset %in% datasets_selected) &
    (model == model_selected)
  ) %>%
  # Add column dataset
  dplyr::left_join(dataset_to_type_df, by = c("type_dataset")) %>%
  dplyr::select(dataset, type_dataset, size, rep, num_edges, model_result)

# Calculate scaling relationship
cols <- c("size", "S", "type_dataset", "S_lag1", "Fn", "slope", "intercept",
          "adj.r.squared")
scaling_relation_df <- data.frame(
  matrix(ncol = length(cols), nrow = 0, dimnames = list(NULL, cols))
)
for (dataset_selected in datasets_selected) {
  scaling_relation_filt <- topology_results_selected_by_size_df %>%
    filter(type_dataset == dataset_selected) %>%
    dplyr::select(size, mean, type_dataset) %>%
    unique() %>%
    dplyr::rename("S" = "mean") %>%
    arrange(size) %>%
    mutate(S_lag1 = if_else(size == min(size), 0, lag(S))) %>%
    mutate(Fn = ((S - S_lag1) / S_lag1)) %>%
    filter((!(is.infinite(Fn))) & (Fn > 0))
  lm_summary <- summary(lm(log(scaling_relation_filt$Fn)~log(scaling_relation_filt$size)))
  slope <- coef(lm_summary)[2]
  intercept <- coef(lm_summary)[1]
  adj_r_squared <- lm_summary$adj.r.squared
  scaling_relation_df <- rbind(
    scaling_relation_df,
    cbind(
      scaling_relation_filt,
      data.frame(
        slope = slope,
        intercept = intercept,
        adj.r.squared = adj_r_squared
      )
    )
  )
}
scaling_relation_df <- scaling_relation_df %>%
  left_join(dataset_to_type_df, by = c("type_dataset")) %>%
  filter(Fn > 0) %>%
  mutate_at(c('adj.r.squared'), ~sprintf("%.2f",.)) %>%
  mutate_at(c('adj.r.squared'), as.character)

# Create scaling relation summary table
scaling_relation_summary_df <- scaling_relation_df %>%
  dplyr::select(type_dataset, adj.r.squared) %>%
  unique() %>%
  arrange(factor(type_dataset, levels = datasets_selected)) %>%
  rename("R**2 (scaling)" = "adj.r.squared") %>%
  mutate_if(is.numeric, ~sprintf("%.2f",.))

# Bind power law with logarithm results
figure_summary_table <- cbind(
  scaling_relation_summary_df,
  (power_law_summary_df %>% dplyr::select(-type_dataset)),
  (log_model_summary_df %>% dplyr::select(-type_dataset))
)
figure_summary_table <- scaling_relation_summary_df %>%
  full_join(power_law_summary_df, by = "type_dataset") %>%
  full_join(log_model_summary_df, by = "type_dataset") %>%
  dplyr::select("dataset", "type_dataset", "R**2 (scaling)", "alpha (model)", "R**2 (model)", "epsilon (model)", "R**2 (Log)", "epsilon (Log)") %>%
  arrange(dataset, type_dataset)

# Read discovery rate results
a_vs_fraction_corr_file <- paste(input_dir, "a_vs_fraction_sig_correlations_pearson_pval_0.05.txt", sep = "/")
sample_size_vs_a_df <- fread(a_vs_fraction_corr_file)

# Read mtcars data
mtcars_mod <- mtcars %>% 
  mutate(car = row.names(mtcars)) %>%
  dplyr::select(car, carb, gear, mpg)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # shinyWidgets::updatePickerInput(session,
  #                                 inputId = "input_cars",
  #                                 label = "Select cars:",
  #                                 choices = unique(mtcars_mod$car),
  #                                 selected = unique(mtcars_mod$car),
  #                                 options = list(
  #                                   `actions-box` = TRUE,
  #                                   size = 10,
  #                                   `selected-text-format` = "count > 1"
  #                                 ))

  process_inputs <- reactive({

    # Get selected datasets
    selected_dataset_types <- c(
      input$type_dataset_gtex,
      input$type_dataset_tcga,
      input$type_dataset_scipher
    )
    print(selected_dataset_types)

    # Filter summary table by dataset
    figure_summary_table_filt <- figure_summary_table %>%
      filter(
        (dataset %in% input$dataset_input) &
        (type_dataset %in% selected_dataset_types)
      ) %>%
      # Rename dataset and type_dataset to human readable names
      mutate(type_dataset = dplyr::recode(type_dataset, !!!name_dict),
             dataset = dplyr::recode(dataset, !!!name_dict))

    # Filter goodnessfit table by dataset
    goodnessfit_filt <- goodnessfit_df %>%
      filter(
        (dataset %in% input$dataset_input) &
        (type_dataset %in% selected_dataset_types)
      ) %>%
      # Include color codes for each dataset
      left_join(
        (color_dict %>%
          as.data.frame() %>%
          dplyr::rename("rgb_col" = ".") %>%
          tibble::rownames_to_column("type_dataset")
        ), by = c("type_dataset")
      ) %>%
      # Rename dataset and type_dataset to human readable names
      mutate(type_dataset = dplyr::recode(type_dataset, !!!name_dict),
             dataset = dplyr::recode(dataset, !!!name_dict),
             # Create new column that pastes type_dataset and adj.r.squared
             r2labels = paste(type_dataset, ": ", adj.r.squared, sep = ""))

    # Filter scaling relation table by dataset
    scaling_relation_filt <- scaling_relation_df %>%
      filter(
        (dataset %in% input$dataset_input) &
        (type_dataset %in% selected_dataset_types)
      ) %>%
      # Include color codes for each dataset
      left_join(
        (color_dict %>%
          as.data.frame() %>%
          dplyr::rename("rgb_col" = ".") %>%
          tibble::rownames_to_column("type_dataset")
        ), by = c("type_dataset")
      ) %>%
      # Rename dataset and type_dataset to human readable names
      mutate(type_dataset = dplyr::recode(type_dataset, !!!name_dict),
             dataset = dplyr::recode(dataset, !!!name_dict),
             # Create new column that pastes type_dataset and adj.r.squared
             r2labels = paste(type_dataset, ": ", adj.r.squared, sep = ""))

    # Filter model prediction table by dataset
    model_prediction_filt <- model_prediction_df %>%
      filter(
        (dataset %in% input$dataset_input) &
        (type_dataset %in% selected_dataset_types)
      ) %>%
      # Include color codes for each dataset
      left_join(
        (color_dict %>%
          as.data.frame() %>%
          dplyr::rename("rgb_col" = ".") %>%
          tibble::rownames_to_column("type_dataset")
        ), by = c("type_dataset")
      ) %>%
      # Rename dataset and type_dataset to human readable names
      mutate(type_dataset = dplyr::recode(type_dataset, !!!name_dict),
             dataset = dplyr::recode(dataset, !!!name_dict))

    # Filter sample size vs. a table by dataset
    type_correlation_selected <- "weak"
    sample_size_vs_a_filt <- sample_size_vs_a_df %>% 
      filter(
        (dataset %in% input$dataset_input) &
        (type_dataset %in% selected_dataset_types)
      ) %>%
      # Filter by correlation type
      filter(type_correlation == type_correlation_selected) %>%
      # Include color codes for each dataset
      left_join(
        (color_dict %>%
          as.data.frame() %>%
          dplyr::rename("rgb_col" = ".") %>%
          tibble::rownames_to_column("type_dataset")
        ), by = c("type_dataset")
      ) %>%
      # Rename dataset and type_dataset to human readable names
      mutate(type_dataset = dplyr::recode(type_dataset, !!!name_dict),
             dataset = dplyr::recode(dataset, !!!name_dict))

    # Calculate correlation and regression between sample size and a
    ss_vs_a_cor <- cor(
      sample_size_vs_a_filt$a,
      sample_size_vs_a_filt$num_edges_from_statistical_corrected_norm
    )
    ss_vs_a_lm <- summary(
      lm(
        sample_size_vs_a_filt$num_edges_from_statistical_corrected_norm ~
        sample_size_vs_a_filt$a
      )
    )

    ss_vs_a_slope <- coef(ss_vs_a_lm)[2]
    ss_vs_a_intercept <- coef(ss_vs_a_lm)[1]
    ss_vs_a_adj_r_squared <- ss_vs_a_lm$adj.r.squared
    ss_vs_a_cor_lm <- paste(
      "R² =",
      round(ss_vs_a_adj_r_squared, 3),
      "\ncorr. =",
      round(ss_vs_a_cor, 3)
    )

    # Initialize the plot variable as NULL
    scaling_plot <- NULL
    discovery_rate_plot <- NULL

    if (input$plot_type == "Corr. vs. size") {

      # Create plot based on goodnessfit table
      scaling_plot <- goodnessfit_filt %>%
        ggplot(aes(x = size, y = num_edges, col = r2labels)) +
        geom_point(alpha = 0.6, size = 3) +
        geom_line(aes(x = size, y = model_result, col = r2labels), linewidth = 1) +
        scale_color_manual(
          guide = "legend",
          values = setNames(goodnessfit_filt$rgb_col, goodnessfit_filt$r2labels)
        ) +
        theme_linedraw() +
        xlab("Num. samples") +
        ylab("Num. significant correlations") +
        guides(col = guide_legend(title = bquote(bold(.("Analytical model") ~ R^2)))) +
        theme(
          aspect.ratio = 1,
          plot.title =  element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 17, face = "bold"),
          axis.text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(family = "Helvetica"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16, face = "bold"),
          legend.background = element_rect(
            fill = "transparent",
            color = "transparent"
          )
        )

    } else if (input$plot_type == "Scaling rel.") {

      scaling_plot <- scaling_relation_filt %>%
        ggplot(aes(x = log(size), y = log(Fn), col = r2labels)) +
        geom_point(alpha = 0.9, size = 3) +
        geom_abline(aes(intercept = intercept, slope = slope, col = r2labels), size = 1, key_glyph = "smooth") +
        theme_linedraw() +
        xlab("log(n)") +
        ylab(bquote(bold(log((S[n]-S[n-1])/S[n-1])))) +
        scale_color_manual(
          guide = "legend",
          values = setNames(scaling_relation_filt$rgb_col, scaling_relation_filt$r2labels)
        ) +
        guides(col = guide_legend(title = bquote(bold(.("Power law") ~ R^2)))) +
        theme(
          aspect.ratio = 1,
          plot.title =  element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 17, face = "bold"),
          axis.text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(family = "Helvetica"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16, face = "bold"),
          legend.background = element_rect(
            fill = "transparent",
            color = "transparent"
          )
        )

    } else if (input$plot_type == "Model prediction") {

      # Create plot based on model prediction table
      scaling_plot <- model_prediction_filt %>%
        ggplot(aes(x = size, y = num_edges, col = type_dataset)) +
        geom_point(alpha = 0.6, size = 3) +
        geom_line(aes(x = size, y = model_result, col = type_dataset), linewidth = 1) +
        scale_color_manual(
          guide = "legend",
          values = setNames(model_prediction_filt$rgb_col, model_prediction_filt$type_dataset)
        ) +
        theme_linedraw() +
        xlab("Num. samples") +
        ylab("Frac. significant correlations") +
        guides(col = guide_legend(title = "Dataset")) +
        theme(
          aspect.ratio = 1,
          plot.title =  element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 17, face = "bold"),
          axis.text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(family = "Helvetica"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16, face = "bold"),
          legend.background = element_rect(
            fill = "transparent",
            color = "transparent"
          )
        )

    } else if (input$plot_type == "Corr. vs. rate") {

      discovery_rate_plot <- sample_size_vs_a_filt %>%
        ggplot(aes(x = a, y = num_edges_from_statistical_corrected_norm)) +
          geom_point(
            aes(col = type_dataset),
            alpha = 0.6,
            size = 3,
            show.legend = TRUE
          ) +
          scale_color_manual(
            guide = "legend",
            values = setNames(sample_size_vs_a_filt$rgb_col, sample_size_vs_a_filt$type_dataset)
          ) +
          geom_abline(aes(slope = ss_vs_a_slope,
                          intercept = ss_vs_a_intercept),
                      linetype = "dashed", 
                      col = "gray50",
                      show.legend = FALSE) +
        annotate("text", x = 2.00, y = 0.45, label = ss_vs_a_cor_lm, size = 5) +
        xlab(expression(alpha)) +
        ylab("Frac. correlations above 0.2") +
        theme_linedraw() +
        theme(
          aspect.ratio = 1,
          plot.title =  element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 17, face = "bold"),
          axis.text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(family = "Helvetica"),
          legend.text = element_text(size = 16),
          legend.title = element_blank(),
          legend.background = element_rect(
            fill = "transparent",
            color = "transparent"
          )
        )

    }

    # mtcars_mod_filt <- mtcars_mod %>%
    #   filter(car %in% input$input_cars)
    # if (input$input_filter) {
    #   input$filter_action
    #   isolate({
    #     mtcars_mod_filt <- mtcars_mod_filt %>%
    #       filter((carb %in% input$input_carb) &
    #                (gear %in% input$input_gear) &
    #                (mpg > input$input_mpg_cutoff)
    #       )
    #   })
    # }
    # if (input$input_merge) {
    #   input$merge_action
    #   isolate({
    #     mtcars_mod_filt <- mtcars_mod_filt %>%
    #       group_by_at(input$merge_cond) %>%
    #       summarize(mpg = median(mpg)) %>%
    #       unite(car, all_of(input$merge_cond), sep = "|")
    #   })
    # }
    return(
      list(
        "summary_table" = figure_summary_table_filt,
        "scaling_plot" = scaling_plot,
        "discovery_rate_plot" = discovery_rate_plot
      )
    )
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

  # Initialize a reactive value to store the current tab
  current_tab <- reactiveVal()

  # Update 'current_tab' whenever 'power_law_tabset' changes
  observeEvent(input$power_law_tabset, {
    current_tab(input$power_law_tabset)
  })

  # Initialize reactive values to store the last selected 'plot_type' for each tab
  last_plot_type <- reactiveValues(
    discovery = "Corr. vs. rate",
    scaling = "Corr. vs. size"
  )

  # Update 'last_plot_type' whenever 'plot_type' changes
  observeEvent(input$plot_type, {
    if (input$power_law_tabset == "Discovery rate plots") {
      last_plot_type$discovery <- input$plot_type
    } else {
      last_plot_type$scaling <- input$plot_type
    }
  })

  # Update 'plot_type' when switching tabs
  observeEvent(input$power_law_tabset, {
    if (input$power_law_tabset == "Discovery rate plots") {
      updateSelectInput(
        session,
        "plot_type",
        choices = c("Corr. vs. rate", "CV vs. rate", "Rate vs. size"),
        selected = last_plot_type$discovery
      )
    } else {
      updateSelectInput(
        session,
        "plot_type",
        choices = c("Corr. vs. size", "Scaling rel.", "Model prediction"),
        selected = last_plot_type$scaling
      )
    }
  })

  output$summaryTable <- DT::renderDataTable({
    results <- process_inputs()
    results$summary_table
  })

  output$scalingPlot <- renderCachedPlot({
    results <- process_inputs()
    results$scaling_plot
  }, cacheKeyExpr = { list(current_tab(), input$plot_type, input$dataset_input, input$type_dataset_gtex, input$type_dataset_tcga, input$type_dataset_scipher) })  # Include all relevant inputs

  output$discoveryRatePlot  <- renderCachedPlot({
    results <- process_inputs()
    results$discovery_rate_plot
  }, cacheKeyExpr = { list(current_tab(), input$plot_type, input$dataset_input, input$type_dataset_gtex, input$type_dataset_tcga, input$type_dataset_scipher) })  # Include all relevant inputs

}

# Run the application
shinyApp(ui = ui, server = server)
