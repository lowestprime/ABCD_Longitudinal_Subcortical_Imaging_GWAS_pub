#### PREP ####
smri.R5.1.baseline.y2.ROC <- readRDS("smri.R5.1.baseline.y2.ROC.rds")
p_load(genpwr)

#### SAMPLE SIZE CALCS ####
# Summary Statistics for smri.R5.1.baseline.y2.ROC
baseline.y2.ROC.stats <- smri.R5.1.baseline.y2.ROC %>%
  summarise(across(starts_with("smri_vol"), list(
    mean = ~ mean(.),
    sd = ~ sd(.)
  ), .names = "{col}_{fn}"))

# Step 1: Rename columns ending with std_dev to std.dev
colnames(baseline.y2.ROC.stats) <- gsub("std_dev$", "std.dev", colnames(baseline.y2.ROC.stats))

# Step 2: Pivot the table
pivoted_df <- baseline.y2.ROC.stats %>%
  pivot_longer(
    cols = everything(),
    names_to = c("phenotype", "stat"),
    names_pattern = "(smri_vol_scs_.*)_(.*)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  select(phenotype, mean, sd)

# Function to calculate sample size using genpwr package
calculate_sample_size2 <- function(sd = sds, maf = maf_seq, es = 0.65, alpha = 5e-8, power = 0.8) {
  # Calculate sample size using ss.calc.linear
  result <- ss.calc.linear(
    power = power,
    MAF = maf,
    ES = es,
    sd_y = sd,  # Using the standard deviation from the data
    Alpha = alpha,
    True.Model = "Additive",
    Test.Model = "Additive"
  )
  return(result)
}

# define variables
maf_seq = seq(0.01, 0.5, 0.01)
sds = pivoted_df$sd

# sample_size2 function test
sample_sizes <- calculate_sample_size2()

# Merge the data frames on 'sd' and 'SD_Y'
merged_df <- left_join(pivoted_df, sample_sizes, by = c("sd" = "SD_Y"))

#### PLOTTING ####
# Create the initial plot with ggplot2
p <- ggplot(merged_df, aes(x = MAF, y = `N_total_at_Alpha_5e-08`, color = phenotype, group = phenotype)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Minimum Sample Size by MAF and Phenotype",
       x = "Minor Allele Frequency (MAF)",
       y = "Minimum Sample Size")

# Convert to plotly object
p_plotly <- ggplotly(p)

# Add checkboxes for selecting phenotypes
# hplot <- highlight(
#   p_plotly,
#   on = "plotly_click",
#   off = "plotly_doubleclick",
#   opacityDim = 0.2,
#   selected = attrs_selected(showlegend = FALSE),
#   dynamic = TRUE,
#   persistent = FALSE
# )

# Customize layout for compactness, normal y-axis ticks, and autoscaling
hplot <- layout(
  p_plotly,
  title = list(text = "Minimum Sample Size by MAF and Phenotype", x = 0.5),
  xaxis = list(title = "Minor Allele Frequency (MAF)"),
  yaxis = list(
    title = "Minimum Sample Size",
    tickformat = ",d",  # Ensure full number format without scientific notation
    ticks = "outside",
    tickmode = "auto",
    nticks = 5,  # Ensure there are at least 4-5 ticks on the y-axis
    rangemode = "auto"
  )
)

# Display the interactive plot
hplot

#### ORIGINAL CODE ####
# Function to calculate sample size using genpwr package
calculate_sample_size <- function(mean, sd, maf = 0.1, alpha = 5e-8, power = 0.8) {
  if (is.numeric(mean) & is.numeric(sd) & sd != 0) {
    effect_size <- mean / sd
    # Calculate sample size using ss.calc.linear
    result <- tryCatch({
      res <- ss.calc.linear(
        power = power,
        MAF = maf,
        ES = effect_size,
        sd_y = sd,  # Using the standard deviation from the data
        Alpha = alpha,
        True.Model = "Additive",
        Test.Model = "Additive"
      )
      res$N_total_at_Alpha_5e-08
    }, error = function(e) {
      return(NA)
    })
    return(result)
  } else {
    return(NA)
  }
}

# Actual sample sizes
actual_sample_sizes <- list(AFR = 1173, AMR = 1445, EUR = 4198)

# Step 3: Calculate the sample size for each phenotype and store it in a new column
sample_sizes <- pivoted_df %>%
  rowwise() %>%
  mutate(sample_size = calculate_sample_size(mean, sd)) %>%
  ungroup()

# Step 4: Check for sufficient power
sample_sizes_with_power <- sample_sizes %>%
  mutate(
    sufficient_power_AFR = ifelse(!is.na(sample_size) & sample_size <= actual_sample_sizes$AFR, TRUE, FALSE),
    sufficient_power_AMR = ifelse(!is.na(sample_size) & sample_size <= actual_sample_sizes$AMR, TRUE, FALSE),
    sufficient_power_EUR = ifelse(!is.na(sample_size) & sample_size <= actual_sample_sizes$EUR, TRUE, FALSE)
  )

# Print the sample sizes with power sufficiency
print(sample_sizes_with_power)