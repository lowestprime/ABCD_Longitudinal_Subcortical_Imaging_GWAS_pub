#!/bin/sh
# Perform Genome-Wide association analysis using GCTA - Females

#$ -wd /u/project/lhernand/sganesh/gwas_srs/association_analysis/related/mlma/bigsnpr/pc20/females
#$ -l highp,h_rt=70:00:00,h_data=8G
#$ -pe shared 16
#$ -o Females_MLMA_all_sor.log
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-5:1

date="05032024"

pops=("EUR" "AMR" "AFR" "EAS" "SAS")
pop=${pops[$SGE_TASK_ID-1]}
phenotypes="sor"

# Directories
indir=/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/splitted_by_ancestry_groups/females
grmDir=/u/project/lhernand/sganesh/gwas_srs/grm/grm_females
phenodataDir=/u/project/lhernand/sganesh/gwas_srs/phenotypes
outdir=/u/project/lhernand/sganesh/gwas_srs/association_analysis/related/mlma/bigsnpr/pc20/females
outFile="$pop.females.$phenotypes.$date.pc20"
covarFile="covar_AllSubj_batch_no_gender_noNAs_baseline_11665.txt"
qcovarFile="qcovar_ABCD5_${pop}_PC20_SOR_related_bigsnpr.txt"

infile="$indir/${pop}.females.genotype"

# Software
gcta=/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1

# Phenotype and covariates
phenotypeFile="Phen_AllSubj_abcd_pssrs01.txt_1yr_followup_ssrs_42_p_within_ancestry_group_noNAs_11023_w_fid.txt"
grmFile="${pop}.females_GRM"

# Run the association analysis using MLMA on pre-computed GRM
$gcta \
    --mlma \
    --bfile $infile \
    --grm $grmDir/$grmFile \
    --pheno $phenodataDir/$phenotypeFile \
    --qcovar $phenodataDir/$qcovarFile \
    --covar $phenodataDir/$covarFile \
    --thread-num 10 \
    --out $outdir/$outFile


