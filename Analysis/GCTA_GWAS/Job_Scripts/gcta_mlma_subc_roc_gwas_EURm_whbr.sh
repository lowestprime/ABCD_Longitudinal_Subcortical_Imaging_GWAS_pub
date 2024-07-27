#!/bin/sh
# Perform Genome-Wide association analysis using GCTA

# Set working directory
#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data

# Request resources
#$ -l highp,h_rt=70:00:00,h_data=5G
#$ -pe shared 32
#$ -l arch=intel-gold*  # Requesting nodes with Intel Gold CPUs

# Output and error handling
#$ -o GCTA_GWAS_EURM_WhBr${JOB_ID}.log
#$ -j y # join std error and std output streams, yes

# Email notifications
#$ -M cobeaman@g.ucla.edu
#$ -m bea # email when job begins, ends, and if it aborts

# Directories
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
grmDir="/u/project/lhernand/cobeaman/grm/grm_all"
results_dir="${base_dir}/Results/test_run"
mkdir -p "${results_dir}"

# Software path
gcta="/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1"

# Date
date=$(date +%Y%m%d)

# Population, Sex, and Phenotype
pop="EUR"
sex="M"
pheno_file="${base_dir}/Phenotypes/EUR/M/988282_pheno_07152024_EUR_M_2268_smri_vol_scs_wholeb_ROC0_2.txt"

# Covariate and Quantitative Covariate Files
covar_file="${base_dir}/Covariates/Discrete/EUR/M/covar_988282_07152024_EUR_M_2268.txt"
qcovar_file="${base_dir}/Covariates/Quantitative/EUR/M/qcovar_noICV_988282_07152024_EUR_M_2268.txt"

# GRM and bfile
num_no_sex_mismatch=$(ls -1v $grmDir/*${pop}*.rmsexmismatch*.fam | grep -v "allproblem" | wc -l)
grm_file_prefix="${grmDir}/ABCD_202209.updated.nodups.curated.cleaned_indivs.qc.basic_withsexinfo_RSid_NoDuplicates_RSidOnly_${pop}.${num_no_sex_mismatch}.rmsexmismatch.0.chr1.22_indep_pruned"

# Output file
pheno_name=$(basename "${pheno_file}" .txt)
num_samples=$(wc -l < "${pheno_file}")
out_file="${results_dir}/${pheno_name}_${date}_n${num_samples}"

# Check if all necessary files exist
if [ ! -f "${pheno_file}" ]; then
  echo "Phenotype file ${pheno_file} does not exist. Exiting."
  exit 1
fi

if [ ! -f "${covar_file}" ]; then
  echo "Covariate file ${covar_file} does not exist. Exiting."
  exit 1
fi

if [ ! -f "${qcovar_file}" ]; then
  echo "Quantitative covariate file ${qcovar_file} does not exist. Exiting."
  exit 1
fi

if [ ! -f "${grm_file_prefix}.bed" ] || [ ! -f "${grm_file_prefix}.bim" ] || [ ! -f "${grm_file_prefix}.fam" ]; then
  echo "One or more of the GRM binary files (${grm_file_prefix}.bed, .bim, .fam) do not exist. Exiting."
  exit 1
fi

if [ ! -f "${grm_file_prefix}_qc_GRM.grm.bin" ] || [ ! -f "${grm_file_prefix}_qc_GRM.grm.id" ] || [ ! -f "${grm_file_prefix}_qc_GRM.grm.N.bin" ]; then
  echo "One or more of the GRM files (${grm_file_prefix}_qc_GRM.grm.bin, .grm.id, .grm.N.bin) do not exist. Exiting."
  exit 1
fi

# Run GCTA MLMA
$gcta --mlma \
      --bfile "${grm_file_prefix}" \
      --grm "${grm_file_prefix}_qc_GRM" \
      --pheno "${pheno_file}" \
      --covar "${covar_file}" \
      --qcovar "${qcovar_file}" \
      --thread-num 32 \
      --out "${out_file}"

if [ $? -ne 0 ]; then
  echo "Error running GCTA MLMA for ${pheno_name}" >&2
  exit 1
fi

echo "Completed GCTA MLMA for ${pheno_name} with ${num_samples} samples"
