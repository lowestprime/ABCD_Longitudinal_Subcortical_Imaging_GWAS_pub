#### Get sex info to add sex column in the genotype file ####
# will be used to update sex column in the plink fam file
rm(list=ls())

# library
library(tidyverse)

#### directories and files ####
dir <- "/u/project/gandalm"
# dir <- "/Users/Emma/Hoffman"
genodir <- paste(dir, "shared/GenomicDatasets/ABCD_Release_4/genomics_sample03/genotype_QCed", sep = "/")
phenodir <- paste(dir, "shared/GenomicDatasets/ABCD_Release_4/phenotypes", sep = "/")
qcdir <- paste(dir, "emmamk/ASD_GWAS/ABCDr4/02.QCplink", sep = "/")


#### Read in longitudinal tracking file ####
setwd(phenodir)
Longit_Tracking <- read.delim("abcd_lt01.txt", header = TRUE, na.strings = c("", "NA"))
Longit_Tracking <- Longit_Tracking[-c(1),]
Longit_Tracking <- Longit_Tracking %>% dplyr::select(subjectkey, eventname, sex, interview_date, interview_age, site_id_l)
nrow(Longit_Tracking)
head(Longit_Tracking)
table(Longit_Tracking$sex)
summary(factor(Longit_Tracking$sex))

# remove NAs in sex column
keep <- !is.na(Longit_Tracking$sex)
table(keep)
# All have sex info

# get baseline year 1
Longit_Tracking_base <- filter(Longit_Tracking, eventname == "baseline_year_1_arm_1")
nrow(Longit_Tracking_base)

# check inconsistent subjectkey found in other files
Longit_Tracking_base$subjectkey[grep("INVF3FYXH1G|INVPWLFYWTX", Longit_Tracking_base$subjectkey)]
# no problems: continue


#### Read in .fam file ####
setwd(genodir)
dt <- read.table("ABCD_release_3.0_QCed.fam", na.strings = c("", "NA"))
head(dt)

#. correct inconsistent subjectkey format
dt[grep("INVF3FYXH1G|INVPWLFYWTX", dt$V2),2]
# '`NDAR_INVF3FYXH1G''NDARINVPWLFYWTX'

dt[grep("INVF3FYXH1G", dt$V2),] <- str_replace(dt[grep("INVF3FYXH1G", dt$V2),], "`NDAR", "NDAR")
dt[grep("INVPWLFYWTX", dt$V2),] <- str_replace(dt[grep("INVPWLFYWTX", dt$V2),], "NDARINVPWLFYWTX", "NDAR_INVPWLFYWTX")
dt[grep("INVF3FYXH1G|INVPWLFYWTX", dt$V2),2]



#### Merge with Longit_Tracking_base ####
dt <- merge(dt, Longit_Tracking_base, by.x = "V2", by.y = "subjectkey", all.x = TRUE)
nrow(dt)
head(dt)

summary(factor(dt$sex))
# F 5214
# M 5849
# NA's 36

#. check if the 36 samples have sex info somewhere else ----
setwd(phenodir)

# abcd_screen01.txt
dt_scrn <- read.delim("abcd_screen01.txt", header = TRUE, na.strings = c("", "NA"))
dt_scrn <- dt_scrn[-c(1),]
names(dt_scrn)

dt_scrn <- dt_scrn %>% dplyr::select(subjectkey, eventname, sex, interview_date, interview_age)
nrow(dt_scrn)

table(dt$V2[is.na(dt$sex)] %in% dt_scrn$subjectkey)
# FALSE 36

# abcd_screen02.txt
dt_scrn <- read.delim("abcd_screen02.txt", header = TRUE, na.strings = c("", "NA"))
dt_scrn <- dt_scrn[-c(1),]
names(dt_scrn)

dt_scrn <- dt_scrn %>% dplyr::select(subjectkey, eventname, sex, interview_date, interview_age)
nrow(dt_scrn)

table(dt$V2[is.na(dt$sex)] %in% dt_scrn$subjectkey)

# conclusion:
# seems these 36 samples does not have sex info
# use Longit_Tracking_base data for sex info


#### Create gender info file for plink --update-sex ####
dt$fid <- paste(dt$V1, dt$V2, sep = "_")
dt$iid <- paste(dt$V1, dt$V2, sep = "_")
dt_plinkformat <- dt %>% dplyr::select(fid, iid, sex)
head(dt_plinkformat)

# replace NA's with 0
dt_plinkformat$sex[is.na(dt_plinkformat$sex)] <- 0

# save file
setwd(qcdir)
filename <- "sex_info.txt"
write.table(dt_plinkformat, filename, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

nrow(dt_plinkformat)
table(dt_plinkformat$sex)
#  0    F    M 
# 36 5214 5849 

