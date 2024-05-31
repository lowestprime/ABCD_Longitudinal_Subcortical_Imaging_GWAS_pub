#!/bin/sh
# Perform Genome-Wide association analysis using GCTA

#$ -wd /u/project/lhernand/sganesh/gwas_srs/association_analysis/related/mlma/pc20
#$ -l highp,h_rt=200:00:00,h_data=8G
#$ -pe shared 16
#$ -o MLMA_allsamples_pc20.log 
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1

date="05032024"
cohort="ABCDr5"
ref="MikeDB_db155_rm"
filebase_1="${cohort}_TopMed_chr1.22_rsID"
filebase=$filebase_1.$ref

pops=("EUR" "AMR" "AFR" "EAS" "SAS")
pop=${pops[$SGE_TASK_ID-1]}

# Directories
indir=/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/splitted_by_ancestry_groups/related_individuals

ancestry_dir=/u/project/lhernand/sganesh/gwas_srs/ancestry_estimation_results/splitted_by_ancestry_groups
grmDir=/u/project/lhernand/sganesh/gwas_srs/grm/grm_all
phenodataDir=/u/project/lhernand/sganesh/gwas_srs/phenotypes
covarFile="covar_AllSubj_batch_gender_noNAs_baseline_11665.txt"
qcovarFile="qcovar_ABCD5_${pop}_PC20_SOR_related_bigsnpr.txt"

individuals="AllSubj"
group="ancestry_group"
timepoint="1yr_followup"

# Software
gcta=/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1

filepath=`ls -1v $ancestry_dir/*${pop}*.rmsexmismatch*.fam | grep -v "allproblem"`
num_no_sex_mismatch=`wc -l $filepath | cut -f1 -d" "`

infile="$indir/$filebase_1.$ref.sexmismatch.${pop}.${num_no_sex_mismatch}"

# Phenotype File
phenotypeFile="Phen_AllSubj_abcd_pssrs01.txt_1yr_followup_ssrs_42_p_within_ancestry_group_noNAs_11023_w_fid.txt"
pheno_short="sor"

grmFile="ABCD_202209.updated.nodups.curated.cleaned_indivs.qc.basic_withsexinfo_RSid_NoDuplicates_RSidOnly_${pop}.${num_no_sex_mismatch}.rmsexmismatch.0.chr1.22_indep_pruned_qc_GRM"

outdir=/u/project/lhernand/sganesh/gwas_srs/association_analysis/related/mlma/bigsnpr/pc20/all
outFile="$pop.$pheno_short.$date.pc20"
mkdir -p $outdir/log

# Run the association analysis using MLMA on pre-computed GRM

$gcta \
    --mlma \
    --bfile $infile \
    --grm $grmDir/$grmFile \
    --pheno $phenodataDir/$phenotypeFile \
    --covar $phenodataDir/$covarFile \
    --qcovar $phenodataDir/$qcovarFile \
    --thread-num 10 \
    --out $outdir/$outFile


