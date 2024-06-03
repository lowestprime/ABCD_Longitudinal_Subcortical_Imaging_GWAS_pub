#### EXTERNAL FUNCTIONS FOR ABCD.Release.5.1.C4.Merged.Data ####
### Cooper Beaman, 5/18/24

# recode timepoint data
recode.timepoint <- Vectorize(function(eventname) {
  parts <- unlist(strsplit(eventname, split = '_'))
  t <- parts[1]
  y <- parts[2]

  if (t == 'baseline') {
    return(0)
  } else if (y == 'year') {
    return(as.numeric(t))
  } else {
    return(as.numeric(t) / 12)
  }
})

# recode sex data
recode.sex <- Vectorize(function(kbi_sex_assigned_at_birth) {
  if (is.na(kbi_sex_assigned_at_birth)) {
    return('NA')
  } else if (kbi_sex_assigned_at_birth == 1) {
    return('M')
  } else if (kbi_sex_assigned_at_birth == 2) {
    return('F')
  } else {
    return('NA')
  }
})

# average hemispheres
average_hemispheres <- function(df) {
  regions <- df %>%
    names() %>%
    str_extract("^smri_vol_scs_.*(?=(lh|rh|l|r)$)") %>%
    na.omit() %>%
    unique()
  avg_data <- sapply(setNames(
    lapply(regions, function(x) {
      c(paste0(x, ifelse(x == "smri_vol_scs_aa", "l", "lh")),
        paste0(x, ifelse(x == "smri_vol_scs_aa", "r", "rh")))
    }),
    regions
  ), function(region) {
    if (all(region %in% names(df)))
      rowMeans(df[, region, drop = FALSE], na.rm = TRUE)
    else
      rep(NA_real_, nrow(df))
  })
  result <- cbind(df, avg_data) %>%
    select(-matches("^smri_vol_scs_.*(lh|rh|aal|aar)$"))
  
  return(result)
}

# Helper function to get latest timepoint information
get_latest_timepoint_info <- function(df, timepoint) {
  df %>%
    filter(timepoint == !!timepoint & sex != 'NA') %>%
    select(src_subject_id, rel_family_id, sex, mri_info_deviceserialnumber, interview_age)
}

# baseline_y2_roc
baseline_y2_roc <- function(df, roc_volumes) {
  # Compute ROCs
  roc_df <- df %>%
    group_by(src_subject_id, rel_family_id) %>%
    reframe(across(
      .cols = all_of(roc_volumes),
      .fns = list(
        ROC0_2 = ~ ((.x[timepoint == 2] - .x[timepoint == 0]) / .x[timepoint == 0]) * (100 / 2)
      ),
      .names = "{.col}_{.fn}"
    ))
  
  # Get latest timepoint information for timepoint 2
  latest_info <- get_latest_timepoint_info(df, 2)
  
  # Join latest timepoint information
  roc_df <- roc_df %>%
    left_join(latest_info, by = c("src_subject_id", "rel_family_id"))
  
  # Filter out rows with 'NA' in the sex column
  roc_df <- roc_df %>%
    filter(sex != 'NA')
  
  return(roc_df)
}

# all_timepoints_roc
all_timepoints_roc <- function(df, roc_volumes) {
  # Compute ROCs
  roc_df <- df %>%
    group_by(src_subject_id, rel_family_id) %>%
    reframe(across(
      .cols = all_of(roc_volumes),
      .fns = list(
        ROC0_2 = ~ ((.x[timepoint == 2] - .x[timepoint == 0]) / .x[timepoint == 0]) * (100 / 2),
        ROC0_4 = ~ ((.x[timepoint == 4] - .x[timepoint == 0]) / .x[timepoint == 0]) * (100 / 4),
        ROC2_4 = ~ ((.x[timepoint == 4] - .x[timepoint == 2]) / .x[timepoint == 2]) * (100 / 2)
      ),
      .names = "{.col}_{.fn}"
    ))
  
  # Get latest timepoint information for timepoint 4
  latest_info <- get_latest_timepoint_info(df, 4)
  
  # Join latest timepoint information
  roc_df <- roc_df %>%
    left_join(latest_info, by = c("src_subject_id", "rel_family_id"))
  
  # Filter out rows with 'NA' in the sex column
  roc_df <- roc_df %>%
    filter(sex != 'NA')
  
  return(roc_df)
}

# pivot_roc_to_long_format.R
pivot_roc_to_long_format <- function(df, is_baseline_y2 = FALSE) {
  df %>%
    pivot_longer(
      cols = -all_of(c("src_subject_id", "rel_family_id")),
      names_to = if (is_baseline_y2) "roi" else c("roi", "Time_Comparison"),
      names_pattern = if (is_baseline_y2) NULL else "(.*)_(ROC.*)",
      values_to = "Value"
    )
}

# pivot_original_to_long_format.R
pivot_original_to_long_format <- function(df, roc_volumes) {
  df %>%
    pivot_longer(
      cols = all_of(roc_volumes),
      names_to = "volume_type",
      values_to = "volume"
    ) %>%
    select(all_of(c("src_subject_id", "rel_family_id", "sex", "interview_age", "eventname", "timepoint", "ethnicity", "volume_type", "volume")))
}
