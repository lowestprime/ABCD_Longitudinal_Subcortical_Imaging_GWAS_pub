#!/bin/sh
# Perform Genome-Wide association analysis using GCTA

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=200:00:00,h_data=5G,arch=intel-gold*  # Memory per process
#$ -pe shared 36
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/GCTA_GWAS_${JOB_ID}.out
#$ -j y
#$ -M cobeaman@g.ucla.edu
#$ -m bea
#$ -t 1-3:1

# Software path
gcta=/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1

# Get current date
date=$(date +%m%d%Y)

# Define directories
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
pheno_dir="${base_dir}/Phenotypes"
covar_dir="${base_dir}/Covariates"
results_dir="${base_dir}/Results"

# Define scratch directory
scratch_base="/u/scratch/c/cobeaman/GCTA_GWAS_${JOB_ID}"
mkdir -p "${scratch_base}"

# Define populations and sexes
pops=("EUR" "AMR" "AFR")
pop=${pops[$SGE_TASK_ID-1]}
sexes=("F" "M")

# Specify combinations to skip (format: "pop_sex_phenotype")
skip_combinations=(
  "EUR_M_smri_vol_scs_wholeb_ROC0_2"
  "AMR_F_smri_vol_scs_wholeb_ROC0_2"
  # Add more combinations to skip as needed
)

# Function to check if a combination should be skipped
should_skip() {
  local pop=$1
  local sex=$2
  local pheno=$3
  local combination="${pop}_${sex}_${pheno}"
  for skip in "${skip_combinations[@]}"; do
    if [[ "$combination" == "$skip" ]]; then
      return 0  # Should skip
    fi
  done
  return 1  # Should not skip
}

# Function to check if all required files exist
check_files() {
  local missing_files=()

  for sex in "${sexes[@]}"; do
    local indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex,,}s"
    local grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex,,}s"
    local grm_file="${grmDir}/${pop}.${sex,,}s_GRM"

    # Check GRM and bfile
    for ext in ".grm.bin" ".grm.id" ".grm.N.bin" ".bed" ".bim" ".fam"; do
      if [ ! -f "${grm_file}${ext}" ]; then
        missing_files+=("${grm_file}${ext}")
      fi
    done

    # Check phenotype files
    local pheno_files=$(ls "${pheno_dir}/${pop}/${sex}/smri_vol_*.txt" 2>/dev/null | grep -v "smri_vol_scs_intracranialv_ROC0_2.txt")
    if [ -z "$pheno_files" ]; then
      missing_files+=("${pheno_dir}/${pop}/${sex}/smri_vol_*.txt")
    fi

    # Check covariate files
    local covar_file="${covar_dir}/Discrete/${pop}/${sex}/covar_discrete.txt"
    local qcovar_file="${covar_dir}/Quantitative/${pop}/${sex}/covar_quant.txt"
    local qcovar_noICV_file="${covar_dir}/Quantitative/${pop}/${sex}/qcovar_noICV.txt"

    for file in "$covar_file" "$qcovar_file" "$qcovar_noICV_file"; do
      if [ ! -f "$file" ]; then
        missing_files+=("$file")
      fi
    done
  done

  if [ ${#missing_files[@]} -ne 0 ]; then
    echo "Error: The following required files are missing:"
    printf '%s\n' "${missing_files[@]}"
    exit 1
  fi

  echo "All required files exist."
}

# Function to run GCTA MLMA
run_gcta_mlma() {
  local pheno_file=$1
  local pheno_name=$(basename "${pheno_file}" .txt)
  local grm_file=$2
  local out_dir=$3
  local date=$4
  local qcovar_file=$5
  local num_samples=$6
  local scratch_dir=$7
  local pop=$8
  local sex=$9

  # Check if this combination should be skipped
  if should_skip "$pop" "$sex" "$pheno_name"; then
    echo "Skipping $pop $sex $pheno_name as specified"
    return 0
  fi

  local out_file="${out_dir}/${pheno_name}_${pop}_${sex}_n${num_samples}_${date}"
  
  # Copy files to scratch directory
  local scratch_grm_file="${scratch_dir}/$(basename ${grm_file})"
  cp "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" "${scratch_dir}/"
  cp "${grm_file}.bed" "${grm_file}.bim" "${grm_file}.fam" "${scratch_dir}/"
  
  echo "Starting GCTA MLMA for ${pheno_name} with ${num_samples} samples..."
  $gcta --mlma \
        --bfile "${scratch_grm_file}" \
        --grm "${scratch_grm_file}" \
        --pheno "${pheno_file}" \
        --covar "${covar_dir}/Discrete/${pop}/${sex}/covar_discrete.txt" \
        --qcovar "${qcovar_file}" \
        --thread-num 36 \
        --out "${out_file}"

  if [ $? -ne 0 ]; then
    echo "Error running GCTA MLMA for ${pheno_name}" >&2
    exit 1
  fi

  echo "Completed GCTA MLMA for ${pheno_name} with ${num_samples} samples"
  
  # Clean up scratch directory
  rm -f "${scratch_dir}"/*
}

export -f run_gcta_mlma
export -f should_skip

# Check all required files before proceeding
check_files

total_tasks=0
completed_tasks=0

# Loop through populations and sexes
for sex in "${sexes[@]}"; do
  echo "Preparing GCTA MLMA tasks for: Population - ${pop}, Sex - ${sex}"

  # Updated GRM and bfile directories
  indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex,,}s"
  grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex,,}s"

  pheno_files=$(ls "${pheno_dir}/${pop}/${sex}/smri_vol_*.txt" | grep -v "smri_vol_scs_intracranialv_ROC0_2.txt")
  filepath=$(ls -1v $grmDir/*${pop}*.rmsexmismatch*.fam | grep -v "allproblem")
  num_no_sex_mismatch=$(wc -l $filepath | cut -f1 -d" ")
  grm_file="${grmDir}/${pop}.${sex,,}s_GRM"

  out_dir="${results_dir}/${sex}/${pop}"
  # create directories for the output dir and logs
  mkdir -p "${out_dir}/log"

  task_list=()
  for pheno_file in ${pheno_files}; do
    num_samples=$(wc -l < "${pheno_file}")
    scratch_dir="${scratch_base}/${pop}_${sex}_${RANDOM}"
    mkdir -p "${scratch_dir}"
    if [[ "${pheno_file}" == *"smri_vol_scs_wholeb_ROC0_2.txt" ]]; then
      task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${covar_dir}/Quantitative/${pop}/${sex}/qcovar_noICV.txt ${num_samples} ${scratch_dir} ${pop} ${sex}")
    else
      task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${covar_dir}/Quantitative/${pop}/${sex}/covar_quant.txt ${num_samples} ${scratch_dir} ${pop} ${sex}")
    fi
    total_tasks=$((total_tasks + 1))
  done

  # Run tasks in parallel
  parallel --jobs 36 --progress run_gcta_mlma ::: "${task_list[@]}"
  
  completed_tasks=$((completed_tasks + ${#task_list[@]}))
  echo "GCTA MLMA completed for: Population - ${pop}, Sex - ${sex}"
  echo "Completed tasks: ${completed_tasks}/${total_tasks}"
done

echo "All GCTA MLMA tasks completed."

# Final cleanup of scratch directory (optional)
rm -rf "${scratch_base}"
echo "Scratch directory cleaned up."