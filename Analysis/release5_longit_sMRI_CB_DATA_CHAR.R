#### LOAD PACKAGES AND ABCD 5.1 DATA, DEFINE DIRECTORIES, AND SOURCE EXTERNAL FUNCTIONS ####

# packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, rlang, lme4, plotrix, ggrepel, gt, ggridges, stringr, reshape2, plotly, shiny, shinybg, data.table)

# directories
base_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/'
image_dir <- paste0(base_dir, 'Images/')
shiny_dir <- paste0(base_dir, 'plots/Shinyapps_plotly_plots/')
table1_dir <- paste0(shiny_dir, 'Counts.Table')
table2_dir <- paste0(shiny_dir, 'ROC.Summary.Table')
gwas_dir <- paste0(base_dir, 'GCTA_GWAS/')
anc_pc_dir <- paste0(gwas_dir, 'ANCESTRY_PCS/')
pheno_dir <- paste0(gwas_dir, 'PHENOTYPE_DATA')

# ABCD 5.1 data and external functions
# load(paste0(base_dir,'ABCD.Release.5.1.C4.Merged.Data.Rda'))
abcdData.R5.1 <- readRDS(paste0(base_dir,'ABCD.Release.5.1.Merged.Data.ICV.rds'))
source(paste0(base_dir, 'release5_external_functions.R'))

#### RECODE SEX AND TIMEPOINT COLUMNS AND DEFINE SAMPLES AND TIMEPOINT FIELDS ####

# recode sex and timepoint columns
abcdData.R5.1 <- abcdData.R5.1 %>%
  mutate(
    timepoint = recode.timepoint(eventname),
    sex = recode.sex(sex)
  )

# define samples and timepoint fields
samples <- unique(abcdData.R5.1$src_subject_id)
id.timepoint.fields <- c('src_subject_id', 'eventname', 'timepoint')

#### SMRI FILTERING AND PREPARATION ####

# filter for smri t1w qc metric and other relevant columns
smri.R5.1 <- abcdData.R5.1 %>%
  filter(imgincl_t1w_include == 1) %>%
  dplyr::select(matches('smri|src_subject_id|eventname|timepoint|ethnicity|pc|ses|interview_age|sex|site_id_l|rel_family_id|mri_info_deviceserialnumber|mrisdp_453|pps')) %>%
  filter(rowSums(is.na(.)) != ncol(.))

# average across brain hemispheres for each applicable scs region
smri.R5.1 <- average_hemispheres(smri.R5.1)

# extract samples with smri data at just baseline and year 2 timepoints
smri.R5.1.baseline.y2 <- smri.R5.1 %>%
  filter(timepoint %in% c(0, 2)) %>%
  group_by(src_subject_id) %>%
  filter(all(c(0, 2) %in% timepoint)) %>%
  ungroup()

# extract samples with smri data at all timepoints
smri.R5.1.all <- smri.R5.1 %>%
  filter(timepoint %in% c(0, 2, 4)) %>%
  group_by(src_subject_id) %>%
  filter(all(c(0, 2, 4) %in% timepoint)) %>%
  ungroup()

#### DATA CHARACTERIZATION  ####

# Define scs columns for volume measurements
roc_volumes <- grep("smri_vol_scs_", names(smri.R5.1.baseline.y2), value = TRUE)

# Calculate the percent rates of change for each scs ROI volume in smri.R5.1.baseline.y2 and smri.R5.1.all
# remove sex == NA and only return values for sex, mri_info_deviceserialnumber, and interview_age from latest timepoint
smri.R5.1.baseline.y2.ROC <- baseline_y2_roc(smri.R5.1.baseline.y2, roc_volumes)
smri.R5.1.all.ROC <- all_timepoints_roc(smri.R5.1.all, roc_volumes)

# convert smri.R5.1.all.ROC and smri.R5.1.all to ggplot2 long format
smri.R5.1.all.ROC.long <- pivot_roc_to_long_format(smri.R5.1.all.ROC, is_baseline_y2 = FALSE)
smri.R5.1.all.long <- pivot_original_to_long_format(smri.R5.1.all, roc_volumes)

# convert smri.R5.1.baseline.y2.ROC and smri.R5.1.baseline.y2 to ggplot2 long format
smri.R5.1.baseline.y2.ROC.long <- pivot_roc_to_long_format(smri.R5.1.baseline.y2.ROC, is_baseline_y2 = TRUE)
smri.R5.1.baseline.y2.long <- pivot_original_to_long_format(smri.R5.1.baseline.y2, roc_volumes)

#### PLOTTING ####

# view sample size tables for smri.R5.1.all and smri.R5.1.baseline.y2 and run in background
port <- httpuv::randomPort(); runBackgroundApp(appDir = table1_dir, port = port); Sys.sleep(1); view_app(port = port)
# kill after viewing
# shinybg::kill_all_apps()

# view ROC summary tables for smri.R5.1.all.ROC.long and smri.R5.1.baseline.y2.ROC.long and run in background
port <- httpuv::randomPort(); runBackgroundApp(appDir = table2_dir, port = port); Sys.sleep(1); view_app(port = port)
# kill after viewing
# shinybg::kill_all_apps()

# create violin/box plots for annual rate of change for volume (toggle all timepoints or baseline y2)
# violin plot of annual rate of change by region
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
            body { text-align: center; margin: auto; }
            .shiny-input-container { margin: 20px auto; width: 50%; } /* Center and adjust width of dropdown */
            .shiny-plot-output { margin: auto; width: 90%; } /* Adjust plot width and centering */
            .js-plotly-plot .plotly .legend { font-size: 16px; } /* Increase legend font size */
        "))
  ),
  titlePanel("Interactive Violin Plots of SCS ROI Rate of Change", windowTitle = "Rate of Change Analysis"),
  selectInput("dataset", "Select Dataset:", choices = c("All Timepoints" = "all", "Baseline + Y2 Only" = "baseline")),
  selectInput("roi_select", "Select ROI:", choices = NULL),  # Choices will be updated based on dataset
  plotlyOutput("roi_plot", height = "600px")
)
# Define server logic
server <- function(input, output, session) {
  # Load data
  all_data <- smri.R5.1.all.ROC.long
  baseline_data <- smri.R5.1.baseline.y2.ROC.long
  # Observe changes in dataset and update ROI selection choices
  observe({
    data <- if (input$dataset == "all") {
      all_data
    } else {
      baseline_data
    }
    updateSelectInput(session, "roi_select", choices = unique(data$roi))
  })
  # Generate plot based on selected dataset and ROI
  output$roi_plot <- renderPlotly({
    data <- if (input$dataset == "all") {
      all_data
    } else {
      baseline_data
    }
    # Ensure ROI is selected
    req(input$roi_select)
    filtered_data <- data[data$roi == input$roi_select, ]
    plot_ly(data = filtered_data, type = "violin",
            x = ~roi, y = ~Value,
            color = ~roi, hoverinfo = 'text',
            text = ~paste('ROI:', roi, '<br>Value:', Value)) %>%
      layout(
        title = paste("Distribution of Rate of Change for", input$roi_select),
        yaxis = list(title = "Rate of Change (mm³/yr)", titlefont = list(size = 16), tickfont = list(size = 12), standoff = 20),
        xaxis = list(title = "ROI", tickangle = 0, tickfont = list(size = 12)),
        legend = list(title = "Time Comparison", orientation = "v", x = 1.05, y = 0.5, xanchor = 'left', yanchor = 'middle', font = list(size = 16)),
        margin = list(l = 40, r = 40, t = 50, b = 50)
      )
  })
}
shinyApp(ui = ui, server = server)

# box plot of region volumes

# Calculate the correlation matrix for scs ROI volume ROCs
# Generate melted correlatoin matrix
melted_correlation <- reshape2::melt(cor(smri.R5.1.baseline.y2.ROC[, roc_volumes], use = "pairwise.complete.obs"))
# Plot correlation matrix heatmap
ggplot(melted_correlation, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title = element_text(size = 14, face = "bold")) +
  labs(x = element_blank(), y = element_blank(), 
       title = "Subcortical ROC Correlation Matrix Heatmap") +
  coord_fixed()
ggsave(filename = paste0(image_dir, 'Subcortical_ROC_Correlation_Matrix_Heatmap.pdf'), width = 8, height = 6)

# Generate pooled and sex-separated histograms for each scs ROI volume ROC
# Loop through each volume to generate and save histograms
for (roc_volume in roc_volumes) {
  # Generate pooled histogram
  Pooled.ROC <- ggplot(smri.R5.1.baseline.y2.ROC, aes_string(x = roc_volume)) +
    geom_histogram(bins = 50, fill = "cornflowerblue", color = "black") +
    ggtitle(paste0('Rate of Change of ', roc_volume, ' (Pooled)')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(Pooled.ROC)
  ggsave(filename = paste0(image_dir, roc_volume, '_Pooled.pdf'), plot = Pooled.ROC, width = 8, height = 6)
  # Generate sex-separated histograms
  Sex.ROC <- ggplot(smri.R5.1.baseline.y2.ROC, aes_string(x = roc_volume)) +
    geom_histogram(bins = 50, fill = "salmon", color = "black") +
    facet_wrap(~sex) +
    ggtitle(paste0('Rate of Change of ', roc_volume, ' (By Sex)')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(image_dir, roc_volume, '_BySex.pdf'), plot = Sex.ROC, width = 8, height = 6)
}

# Check the normality of these distributions
for (roc_volume in roc_volumes) {
  print(paste0("Shapiro-Wilk Test for ", roc_volume, " (Pooled):"))
  print(shapiro.test(smri.R5.1.baseline.y2.ROC[[roc_volume]]))
  
  print(paste0("Shapiro-Wilk Test for ", roc_volume, " (Male):"))
  print(shapiro.test(smri.R5.1.baseline.y2.ROC %>% filter(sex == "M")[[roc_volume]]))
  
  print(paste0("Shapiro-Wilk Test for ", roc_volume, " (Female):"))
  print(shapiro.test(smri.R5.1.baseline.y2.ROC %>% filter(sex == "F")[[roc_volume]]))
}

