### external functions for release5_longit_mri_C4.R
### Jack Dodson, 3/20/24

# recode timepoint data
# input = eventname e.g "2_year_follow_up_y_arm_1"
# output = numeric timepoint in years
recode.timepoint <- function(eventname) {
  t <- unlist(strsplit(eventname, split = '_'))[1];
  y <- unlist(strsplit(eventname, split = '_'))[2];
  if (t == 'baseline') {
    return(0);
  }
  else if (y == 'year') {
    return(as.numeric(t));
  }
  else {
    return(as.numeric(t) / 12);
  }
}

# recode sex data
# input = kbi_sex_assigned_at_birth [1,2,777,999]
# output = ['M', 'F', 'NA']
recode.sex <- function(kbi_sex_assigned_at_birth) {
  if (!is.na(kbi_sex_assigned_at_birth) && kbi_sex_assigned_at_birth == 1) {
    return('M');
  }
  else if (!is.na(kbi_sex_assigned_at_birth) && kbi_sex_assigned_at_birth == 2) {
    return('F');
  }
  else{
    return('NA');
  }
}

# Create PheWAS plots
plot.pheWAS <- function(df) {
  plot <- ggplot(df, aes(x=Category, y=-log10(P_C4A), shape=Direction)) + 
    geom_point(aes(col=Category, size=abs(beta_C4A), fill=Category), position = position_jitter(seed = 1)) + 
    scale_shape_manual(values=c(25,24)) +
    theme_classic() + 
    theme(text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks=element_blank()) + 
    labs(color="Category", size="Abs. Beta", x="Phenotypes", y="-log10(P)") + 
    geom_label_repel(aes(label = ifelse(-log10(P_C4A) > 2, ID, "")), #label only significant points
                     size=4,
                     box.padding = 0.5,
                     force = 1,
                     force_pull = 1,
                     min.segment.length = 0,
                     point.padding = 0,
                     direction = "both",
                     segment.linetype = 6,
                     segment.color = "black",
                     segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc")),
                     max.overlaps = Inf,
                     position = position_jitter(seed = 1)
    ) + 
    geom_hline(yintercept=-log10(0.05), color="red", linewidth=1, alpha=0.3, linetype="dashed") +
    ggtitle(paste("C4A PheWAS", df$data[1]))
  return(plot);
}

create.severity.number.dataframe <- function(timepoint.data, pheno) {
  male.timepoint <- filter(timepoint.data, sex == "M");
  female.timepoint <- filter(timepoint.data, sex == "F");
  
  pheno.male <- male.timepoint[[pheno]];
  pheno.female <- female.timepoint[[pheno]];
  
  # Male
  quantiles.male <- quantile(
    x = pheno.male,
    na.rm = TRUE
    );
  quantiles.female <- quantile(
    x = pheno.female,
    na.rm = TRUE
    );
  
  BigBrain <- subset(male.timepoint, subset = get(pheno) >= quantiles.male[4]);
  BigBrain$Size <- "Top 25%"
  
  SmallBrain <- subset(male.timepoint, subset = get(pheno) <= quantiles.male[2])
  SmallBrain$Size <- "Bottom 25%"
  
  mean <- mean(BigBrain$pps_y_ss_number, na.rm=TRUE)
  sd <- sd(BigBrain$pps_y_ss_number, na.rm=TRUE)
  se <- plotrix::std.error(BigBrain$pps_y_ss_number, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_number, SmallBrain$pps_y_ss_number, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  wilcox$conf.int
  big.df.male.number <- data.frame(mean, se, p)
  big.df.male.number$Sex <- "Male"
  big.df.male.number$Size <- "Top 25%"
  big.df.male.number$Variable <- "Number of Events"
  
  mean <- mean(SmallBrain$pps_y_ss_number, na.rm=TRUE)
  sd <- sd(SmallBrain$pps_y_ss_number, na.rm=TRUE)
  se <- plotrix::std.error(SmallBrain$pps_y_ss_number, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_number, SmallBrain$pps_y_ss_number, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  wilcox$conf.int
  small.df.male.number <- data.frame(mean, se, p)
  small.df.male.number$Sex <- "Male"
  small.df.male.number$Size <- "Bottom 25%"
  small.df.male.number$Variable <- "Number of Events"
  
  mean <- mean(BigBrain$pps_y_ss_severity_score, na.rm=TRUE)
  sd <- sd(BigBrain$pps_y_ss_severity_score, na.rm=TRUE)
  se <- plotrix::std.error(BigBrain$pps_y_ss_severity_score, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_severity_score, SmallBrain$pps_y_ss_severity_score, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  wilcox$conf.int
  big.df.male.severity <- data.frame(mean, se, p)
  big.df.male.severity$Sex <- "Male"
  big.df.male.severity$Size <- "Top 25%"
  big.df.male.severity$Variable <- "Severity of Events"
  
  mean <- mean(SmallBrain$pps_y_ss_severity_score, na.rm=TRUE)
  sd <- sd(SmallBrain$pps_y_ss_severity_score, na.rm=TRUE)
  se <- plotrix::std.error(SmallBrain$pps_y_ss_severity_score, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_severity_score, SmallBrain$pps_y_ss_severity_score, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  small.df.male.severity <- data.frame(mean, se, p)
  small.df.male.severity$Sex <- "Male"
  small.df.male.severity$Size <- "Bottom 25%"
  small.df.male.severity$Variable <- "Severity of Events"
  
  df.male.number <- rbind(big.df.male.number, small.df.male.number)
  df.male.severity <- rbind(big.df.male.severity, small.df.male.severity)
  df.male.big.small <- rbind (BigBrain, SmallBrain)
  
  male.iqr.number <- df.male.big.small %>%
    group_by(Size) %>%
    get_summary_stats(pps_y_ss_number, type = "median_iqr") 

  male.test.number <- df.male.big.small %>% 
    wilcox.test(pps_y_ss_number ~ Size, data = .) %>%
    rstatix::add_significance()
  
  male.iqr.severity <- df.male.big.small %>%
    group_by(Size) %>%
    get_summary_stats(pps_y_ss_severity_score, type = "median_iqr")
  
  male.test.severity <- df.male.big.small %>% 
    wilcox.test(pps_y_ss_severity_score ~ Size, data = .) %>%
    rstatix::add_significance()
  
  # Female
  BigBrain <- subset(female.timepoint, subset = get(pheno) >= quantiles.female[4]);
  BigBrain$Size <- "Top 25%"
  
  SmallBrain <- subset(female.timepoint, subset = get(pheno) <= quantiles.female[2])
  SmallBrain$Size <- "Bottom 25%"
  
  mean <- mean(BigBrain$pps_y_ss_number, na.rm=TRUE)
  sd <- sd(BigBrain$pps_y_ss_number, na.rm=TRUE)
  se <- plotrix::std.error(BigBrain$pps_y_ss_number, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_number, SmallBrain$pps_y_ss_number, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  big.df.female.number <- data.frame(mean, se, p)
  big.df.female.number$Sex <- "Female"
  big.df.female.number$Size <- "Top 25%"
  big.df.female.number$Variable <- "Number of Events"
  
  mean <- mean(SmallBrain$pps_y_ss_number, na.rm=TRUE)
  sd <- sd(SmallBrain$pps_y_ss_number, na.rm=TRUE)
  se <- plotrix::std.error(SmallBrain$pps_y_ss_number, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_number, SmallBrain$pps_y_ss_number, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  small.df.female.number <- data.frame(mean, se, p)
  small.df.female.number$Sex <- "Female"
  small.df.female.number$Size <- "Bottom 25%"
  small.df.female.number$Variable <- "Number of Events"
  
  mean <- mean(BigBrain$pps_y_ss_severity_score, na.rm=TRUE)
  sd <- sd(BigBrain$pps_y_ss_severity_score, na.rm=TRUE)
  se <- plotrix::std.error(BigBrain$pps_y_ss_severity_score, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_severity_score, SmallBrain$pps_y_ss_severity_score, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  big.df.female.severity <- data.frame(mean, se, p)
  big.df.female.severity$Sex <- "Female"
  big.df.female.severity$Size <- "Top 25%"
  big.df.female.severity$Variable <- "Severity of Events"
  
  mean <- mean(SmallBrain$pps_y_ss_severity_score, na.rm=TRUE)
  sd <- sd(SmallBrain$pps_y_ss_severity_score, na.rm=TRUE)
  se <- plotrix::std.error(SmallBrain$pps_y_ss_severity_score, na.rm=TRUE)
  wilcox <- wilcox.test(BigBrain$pps_y_ss_severity_score, SmallBrain$pps_y_ss_severity_score, conf.int = TRUE, conf.level = 0.95)
  p <- wilcox$p.value
  small.df.female.severity <- data.frame(mean, se, p)
  small.df.female.severity$Sex <- "Female"
  small.df.female.severity$Size <- "Bottom 25%"
  small.df.female.severity$Variable <- "Severity of Events"
  
  df.female.number <- rbind(big.df.female.number, small.df.female.number)
  df.female.severity <- rbind(big.df.female.severity, small.df.female.severity)
  df.female.big.small <- rbind (BigBrain, SmallBrain)
  
  female.iqr.number <- df.female.big.small %>%
    group_by(Size) %>%
    get_summary_stats(pps_y_ss_number, type = "median_iqr") 

  female.test.number <- df.female.big.small %>% 
    wilcox.test(pps_y_ss_number ~ Size, data = .) %>%
    add_significance()

  female.iqr.severity <- df.female.big.small %>%
    group_by(Size) %>%
    get_summary_stats(pps_y_ss_severity_score, type = "median_iqr")

  female.test.severity <- df.female.big.small %>% 
    wilcox.test(pps_y_ss_severity_score ~ Size, data = .) %>%
    add_significance()

  # Combine
  df.number <- rbind(df.male.number, df.female.number)
  df.severity <- rbind(df.male.severity, df.female.severity)
  df.year1 <- rbind(df.number, df.severity)
  return(df.year1);
}
