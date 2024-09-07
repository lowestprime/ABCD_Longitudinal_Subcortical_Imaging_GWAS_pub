#!/bin/sh
# Perform Genome-Wide association analysis using GCTA

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=70:00:00,h_data=8G
#$ -pe shared 16
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Results/Processed_Data/Results/GCTA_GWAS_$JOB_ID.out
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-5:1

# Load the GNU Parallel module
. /u/local/Modules/default/init/modules.sh
module load parallel

# Debug step: Check if 'parallel' is working
if ! command -v parallel &> /dev/null; then
  echo "Error: GNU Parallel is not installed or not functioning correctly."
  exit 1
fi

echo "GNU Parallel is working correctly. Proceeding with the analysis..."

# Software path
gcta=/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1

# Get current date
date=$(date +"%m%d%Y")

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

# Function to run GCTA MLMA
run_gcta_mlma() {
  local pheno_file=$1
  local pheno_name=$(basename "${pheno_file}" .txt)
  local grm_file=$2
  local out_dir=$3
  local date=$4
  local covar_file=$5
  local qcovar_file=$6
  local num_samples=$7
  local scratch_dir=$8
  local pop=$9
  local sex=${10}

  # Skip specified combinations
  if should_skip "$pop" "$sex" "$pheno_name"; then
    echo "Skipping $pop $sex $pheno_name"
    return 0
  fi

  local out_file="${out_dir}/${pheno_name}_${pop}_${sex}_n${num_samples}_${date}_${JOB_ID}"

  # Transfer GRM and bfile to scratch directory
  scratch_infile="${scratch_dir}/$(basename ${grm_file})"
  cp "${grm_file}.bed" "${grm_file}.bim" "${grm_file}.fam" "${scratch_dir}/"
  cp "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" "${scratch_dir}/"

  # Run GCTA
  $gcta --mlma \
        --bfile "${scratch_infile}" \
        --grm "${scratch_infile}" \
        --pheno "${pheno_file}" \
        --covar "${covar_file}" \
        --qcovar "${qcovar_file}" \
        --thread-num 16 \
        --out "${out_file}"

  # Check for errors
  if [ $? -ne 0 ]; then
    echo "Error running GCTA MLMA for ${pheno_name}" >&2
    exit 1
  fi

  # Clean up scratch directory
  rm -f "${scratch_dir}"/*
}

export -f run_gcta_mlma

# Prepare tasks
total_tasks=0
completed_tasks=0
for sex in "${sexes[@]}"; do
  echo "Preparing tasks for Population: $pop, Sex: $sex"
  
  # Check and convert 'M' to 'males' and 'F' to 'females'
  if [ "$sex" = "M" ]; then
    sex_dir="males"
  elif [ "$sex" = "F" ]; then
    sex_dir="females"
  fi

  # Update GRM and bfile directories
  indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}"
  grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}"

  pheno_files=$(find "${pheno_dir}/${pop}/${sex}/" -type f -name "smri_vol_*.txt" -not -path "*/archive/*")
  grm_file="${grmDir}/${pop}.${sex_dir}_GRM"

  out_dir="${results_dir}/${sex}/${pop}"
  mkdir -p "${out_dir}/log"

  task_list=()
  for pheno_file in ${pheno_files}; do
    num_samples=$(wc -l < "${pheno_file}")
    scratch_dir="${scratch_base}/${pop}_${sex}_${RANDOM}"
    mkdir -p "${scratch_dir}"
    
    # Find the covariate files
    covar_file=$(find "${covar_dir}/Discrete/${pop}/${sex}/" -type f -name "covar_*${pop}_${sex}*.txt" -not -path "*/archive/*")
    
    # Find the correct quantitative covariate file
    if [[ "$(basename ${pheno_file} .txt)" == "smri_vol_scs_wholeb_ROC0_2" ]]; then
      qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -type f -name "qcovar_noICV_*${pop}_${sex}*.txt")
    else
      qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -type f -name "qcovar_*${pop}_${sex}*.txt")
    fi

    task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${covar_file} ${qcovar_file} ${num_samples} ${scratch_dir} ${pop} ${sex}")
    total_tasks=$((total_tasks + 1))
  done

  # Run tasks in parallel
  parallel --jobs 16 run_gcta_mlma ::: "${task_list[@]}"
  completed_tasks=$((completed_tasks + ${#task_list[@]}))

  echo "Tasks completed: ${completed_tasks}/${total_tasks}"
done

# Clean up scratch directory
rm -rf "${scratch_base}"
echo "Cleanup complete."
