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
- [ ] 1.	Plot histograms to ensure proper rank inverse log normalize AFTER sex and ethnicity split
- [ ] 2.	Add to txt naming info to have sruthi/emma/leanna info for each in pheno prep script (date/number of samples/phenonames)
- [ ] 3.	Clean and simplify directory structure to account for the use of master covar and just the two qcovar txt options
- [ ] 4.	Make the separate qcovar file for wholebrain so one with ICV and one without
- [ ] 5.	Check grm files to ensure they have all required samples in the pheno and other files
- [ ] 6.	Check the txt dummy batch and MRI serial variable col names
- [ ] 7.	Check txts are in format that gcta mlma can read
- [ ] 8.	Improve/finalize job script
- [ ] 9.	Finalize split txt formatting, file names and directory structure
- [ ] 10.	Switch from split to master covar and qcovar files and adjust job script accordingly
- [ ] 11.	Check dummy variable formatting for batch and MRI serial cols in txt splitting script
- [ ] 12.	Check input txt save format same as Sruthi/Emma and compatible with GCTA MLMA
- [ ] 13.	Check contents of each split txt to confirm successful splitting and expected contents
- [ ] 14. Make sure you can add covars and qcovars as single file instead of split

### **Stage 2 - Longitudinal Subcortical Volume GWAS**
  - [MLMA](https://yanglab.westlake.edu.cn/software/gcta/#MLMA)
  - [SAIGE](https://saigegit.github.io/SAIGE-doc/docs/single.html)

### **Stage 3 - Adult Neuropsychiatric GWAS Integration**
  - Condition existing neurodevelopmental disorder PRS on ROCs using [PleioPGS](https://www.biologicalpsychiatryjournal.com/article/S0006-3223(21)01865-5)
  - Perform [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) to investigate joint-genetic architectures and ROC mediation

