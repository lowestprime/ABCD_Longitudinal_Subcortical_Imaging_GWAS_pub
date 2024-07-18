# **ABCD Subcortical Volume Rate of Change GWAS**

## Motivation and Background
  - [Relevant Literature+](https://docs.google.com/spreadsheets/d/1daRx5JcFafNdxd7xn3jf4QrojfYBzgE11LgFUIyP4KY)

## Investigation Overview

### **Stage 1 - Data characterization and preparation**
#### _Summary Tables_
  - [Overview Table and Barplots of Sample Sizes by Timepoint Ethnicity and Sex](https://lowestprime.shinyapps.io/Ethnicity_and_Sex_Counts_by_Timepoint/)
  - [Overview Table of Subcortical Volume ROCs](https://lowestprime.shinyapps.io/ROC_Summary_Table/)

#### _Phenotype Distributions_
  - [Violin Plots of Subcortical Volume ROCs](https://lowestprime.shinyapps.io/Interactive_SCS_ROI_ROC_Violin_Plots_y0_2/)
  - [Box Plots of Subcortical Volumes](https://lowestprime.shinyapps.io/Interactive_SCS_ROI_Volume_Box_Plots_y0_2/)

### _Tasks_
- [ ] Find out if relatedness needs to be accounted for in GCTA setup (crossref with Sruthi and Emma phenoprep scripts)
- [ ] Check if GRM files and file names are appropriate for this analysis and that contents are compatible with job script and split txts
- [ ] Optional Phenotype Splitting Functionality
  - [ ] Add histogram/transformed phenotype data QC as toggleable arg in save_split_data function to be called within split txt prep script function.
- [ ] Finalize job script
- [ ] Optional Additional Normalization QC
  - [ ] Add Correlation Analysis
  - [ ] Pairwise Comparisons
  - [ ] Linear Models and Covariance Checks
  - [ ] Analysis of Covariance (ANCOVA)
  - [ ] Phenotypes-Covariate Relationship Visualizations
  - [ ] Phenotypes-Coefficient Association Analsyses and Rates of Change Comparisons

### **Stage 2 - Longitudinal Subcortical Volume GWAS**
  - [MLMA](https://yanglab.westlake.edu.cn/software/gcta/#MLMA)
  - [SAIGE](https://saigegit.github.io/SAIGE-doc/docs/single.html)

### **Stage 3 - Adult Neuropsychiatric GWAS Integration**
  - Condition existing neurodevelopmental disorder PRS on ROCs using [PleioPGS](https://www.biologicalpsychiatryjournal.com/article/S0006-3223(21)01865-5)
  - Perform [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) to investigate joint-genetic architectures and ROC mediation

