#!/bin/sh
# Perform Genome-Wide association analysis using GCTA

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=200:00:00,h_data=16G  # Memory per process
#$ -pe shared 32  # Using 32 cores
#$ -o GCTA_GWAS_${JOB_ID}.out
#$ -j y
#$ -M cobeaman@g.ucla.edu
#$ -m bea
#$ -t 1-3:1

# Software path
gcta=/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1

date=$(date +%Y%m%d)
cohort="ABCDr5.1"
ref="MikeDB_db155_rm"
filebase_1="${cohort}_TopMed_chr1.22_rsID"
filebase=$filebase_1.$ref

# Define directories
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
pheno_dir="${base_dir}/Phenotypes"
covar_dir="${base_dir}/Covariates"
grmDir="/u/project/lhernand/sganesh/gwas_srs/grm/grm_all"
results_dir="${base_dir}/results"

# Define populations and sexes
pops=("EUR" "AMR" "AFR")
pop=${pops[$SGE_TASK_ID-1]}
sexes=("F" "M")

# Function to run GCTA MLMA
run_gcta_mlma() {
  local pheno_file=$1
  local pheno_name=$(basename "${pheno_file}" .txt)
  local grm_file=$2
  local out_dir=$3
  local date=$4
  local qcovar_file=$5
  local num_samples=$6

  local out_file="${out_dir}/${pheno_name}_${date}_n${num_samples}"
  
  echo "Starting GCTA MLMA for ${pheno_name} with ${num_samples} samples..."
  $gcta --mlma \
        --bfile "${grm_file}" \
        --grm "${grm_file}" \
        --pheno "${pheno_file}" \
        --covar "${covar_dir}/Discrete/${pop}/${sex}/covar_discrete.txt" \
        --qcovar "${qcovar_file}" \
        --thread-num 32 \
        --out "${out_file}"

  if [ $? -ne 0 ]; then
    echo "Error running GCTA MLMA for ${pheno_name}" >&2
    exit 1
  fi

  echo "Completed GCTA MLMA for ${pheno_name} with ${num_samples} samples"
}

export -f run_gcta_mlma

total_tasks=0
completed_tasks=0

# Loop through populations and sexes
for sex in "${sexes[@]}"; do
  echo "Preparing GCTA MLMA tasks for: Population - ${pop}, Sex - ${sex}"

  pheno_files=$(ls "${pheno_dir}/${pop}/${sex}/smri_vol_*.txt" | grep -v "smri_vol_scs_intracranialv_ROC0_2.txt")
  filepath=$(ls -1v $grmDir/*${pop}*.rmsexmismatch*.fam | grep -v "allproblem")
  num_no_sex_mismatch=$(wc -l $filepath | cut -f1 -d" ")
  grm_file="${grmDir}/ABCD_202209.updated.nodups.curated.cleaned_indivs.qc.basic_withsexinfo_RSid_NoDuplicates_RSidOnly_${pop}.${num_no_sex_mismatch}.rmsexmismatch.0.chr1.22_indep_pruned_qc_GRM"

  out_dir="${results_dir}/${sex}/${pop}"
  mkdir -p "${out_dir}/log"

  task_list=()
  for pheno_file in ${pheno_files}; do
    num_samples=$(wc -l < "${pheno_file}")
    if [[ "${pheno_file}" == *"smri_vol_scs_wholeb_ROC0_2.txt" ]]; then
      task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${covar_dir}/Quantitative/${pop}/${sex}/qcovar_noICV.txt ${num_samples}")
    else
      task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${covar_dir}/Quantitative/${pop}/${sex}/covar_quant.txt ${num_samples}")
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