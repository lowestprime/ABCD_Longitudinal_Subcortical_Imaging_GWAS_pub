#!/bin/sh
# Perform Genome-Wide association analysis using GCTA

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=200:00:00,h_data=16G  # Memory per process
#$ -pe shared 32  # Using 32 cores
#$ -o GCTA_GWAS_${JOB_ID}.out
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-3:1

# Software path
gcta=/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1

date=$(date +%Y%m%d)
cohort="ABCDr5.1"
ref="TOPMed"

# Define directories
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
pheno_dir="${base_dir}/Phenotypes"
covar_dir="${base_dir}/Covariates"
grm_dir="/u/project/lhernand/sganesh/gwas_srs/grm" # Adjust if different
results_dir="${base_dir}/results"

# Define populations and sexes
pops=("EUR" "AMR" "AFR")
pop=${pops[$SGE_TASK_ID-1]}
sexes=("F" "M")

# Master covariate files
master_covar_file="${covar_dir}/covar_discrete_master.txt"
master_qcovar_file="${covar_dir}/covar_quant_master.txt"
master_qcovar_file_no_icv="${covar_dir}/covar_quant_master_no_icv.txt"

# Function to run GCTA MLMA
run_gcta_mlma() {
  local pheno_file=$1
  local pheno_name=$(basename "${pheno_file}" .txt)
  local grm_file=$2
  local out_dir=$3
  local date=$4
  local qcovar_file=$5

  local out_file="${out_dir}/${pheno_name}_${date}"
  
  echo "Starting GCTA MLMA for ${pheno_name}..."
  $gcta --mlma \
        --bfile "${grm_file}" \
        --grm-sparse "${grm_file}" \
        --pheno "${pheno_file}" \
        --covar "${master_covar_file}" \
        --qcovar "${qcovar_file}" \
        --thread-num 32 \
        --out "${out_file}"

  gzip "${out_file}.fastGWA"
  echo "Completed GCTA MLMA for ${pheno_name}"
}

export -f run_gcta_mlma

total_tasks=0
completed_tasks=0

# Loop through populations and sexes
for sex in "${sexes[@]}"; do
  echo "Preparing GCTA MLMA tasks for: Population - ${pop}, Sex - ${sex}"

  pheno_files=$(ls "${pheno_dir}/${pop}/${sex}/smri_vol_*.txt")
  grm_file=$(ls "${grm_dir}/*${pop}*_GRMsparse.grm.id" | sed 's/.grm.id//g')

  out_dir="${results_dir}/${sex}/${pop}"
  mkdir -p "${out_dir}/log"

  task_list=()
  for pheno_file in ${pheno_files}; do
    if [[ "${pheno_file}" == *"smri_vol_scs_wholeb_ROC0_2.txt" ]]; then
      task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${master_qcovar_file_no_icv}")
    else
      task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${master_qcovar_file}")
    fi
    total_tasks=$((total_tasks + 1))
  done

  # Run tasks in parallel
  parallel --jobs 32 --progress run_gcta_mlma ::: "${task_list[@]}"
  
  completed_tasks=$((completed_tasks + ${#task_list[@]}))
  echo "GCTA MLMA completed for: Population - ${pop}, Sex - ${sex}"
  echo "Completed tasks: ${completed_tasks}/${total_tasks}"
done

echo "All GCTA MLMA tasks completed."
