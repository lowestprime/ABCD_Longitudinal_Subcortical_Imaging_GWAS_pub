#!/bin/sh
# Perform Genome-Wide association analysis using GCTA --mlma for all jobs in parallel

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=70:00:00,h_data=8G
#$ -pe shared 16
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run/GCTA_GWAS_$JOB_ID.out
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-3:1

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

# Function to check if all required files exist
check_files_exist() {
  local pheno_file=$1
  local grm_file=$2
  local covar_file=$3
  local qcovar_file=$4
  local infile=$5

  local missing_files=()

  # Check phenotype file
  if [ ! -f "$pheno_file" ]; then
    missing_files+=("Phenotype file: $pheno_file")
  fi

  # Check GRM files
  for ext in .grm.bin .grm.id .grm.N.bin; do
    if [ ! -f "${grm_file}${ext}" ]; then
      missing_files+=("GRM file: ${grm_file}${ext}")
    fi
  done

  # Check bfile files
  for ext in .bed .bim .fam; do
    if [ ! -f "${infile}${ext}" ]; then
      missing_files+=("Bfile: ${infile}${ext}")
    fi
  done

  # Check covariate file
  if [ ! -f "$covar_file" ]; then
    missing_files+=("Covariate file: $covar_file")
  fi

  # Check quantitative covariate file
  if [ ! -f "$qcovar_file" ]; then
    missing_files+=("Quantitative covariate file: $qcovar_file")
  fi

  # If any files are missing, print them and return 1 (error)
  if [ ${#missing_files[@]} -ne 0 ]; then
    echo "Error: The following required files are missing:"
    printf '%s\n' "${missing_files[@]}"
    return 1
  fi

  return 0  # All files exist
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
  local infile=${11}

  # Skip specified combinations
  if should_skip "$pop" "$sex" "$pheno_name"; then
    echo "Skipping $pop $sex $pheno_name"
    return 0
  fi

  # Check if all required files exist
  if ! check_files_exist "$pheno_file" "$grm_file" "$covar_file" "$qcovar_file" "$infile"; then
    echo "Skipping analysis for $pop $sex $pheno_name due to missing files"
    return 1
  fi

  local out_file="${out_dir}/${pheno_name}_${pop}_${sex}_n${num_samples}_${date}_${JOB_ID}"

  # Transfer GRM and bfile to scratch directory
  scratch_bfile="${scratch_dir}/$(basename ${infile})"
  scratch_grm="${scratch_dir}/$(basename ${grm_file})"
  cp "${infile}.bed" "${infile}.bim" "${infile}.fam" "${scratch_dir}/"
  cp "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" "${scratch_dir}/"

  # Run GCTA
  $gcta --mlma \
        --bfile "${scratch_bfile}" \
        --grm "${scratch_grm}" \
        --pheno "${pheno_file}" \
        --covar "${covar_file}" \
        --qcovar "${qcovar_file}" \
        --thread-num 16 \
        --out "${out_file}"

  # Check for errors
  if [ $? -ne 0 ]; then
    echo "Error running GCTA MLMA for ${pheno_name}" >&2
    return 1
  fi

  # Clean up scratch directory
  rm -f "${scratch_dir}"/*
}

export -f run_gcta_mlma check_files_exist should_skip

# Prepare tasks
total_tasks=0
completed_tasks=0

echo "Current population: $pop"

for sex in "${sexes[@]}"; do
  echo "Preparing tasks for Population: $pop, Sex: $sex"
  
  # Convert 'M' to 'males' and 'F' to 'females' for indir and grmDir
  if [ "$sex" = "M" ]; then
    sex_dir="males"
  elif [ "$sex" = "F" ]; then
    sex_dir="females"
  fi

  # Updated GRM and bfile directories
  indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}"
  grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}"

  # Ensure the population and sex directories are properly set
  if [ -z "${pop}" ] || [ -z "${sex}" ]; then
    echo "Error: Population or sex is not set correctly. Skipping."
    continue
  fi

  # Debug: print the directories being used
  echo "Phenotype directory: ${pheno_dir}/${pop}/${sex}"
  echo "Discrete Covariate directory: ${covar_dir}/Discrete/${pop}/${sex}"
  echo "Quantitative Covariate directory: ${covar_dir}/Quantitative/${pop}/${sex}"
  echo "bfile directory: ${indir}"
  echo "GRM directory: ${grmDir}"

  # Find phenotype files and check for errors
  pheno_files=$(find "${pheno_dir}/${pop}/${sex}/" -type f -name "*_pheno_*_${pop}_${sex}_*_smri_vol_*.txt" -not -path "*/archive/*")
  if [ -z "${pheno_files}" ]; then
    echo "Error: No phenotype files found in ${pheno_dir}/${pop}/${sex}/"
    continue  # Skip to next sex/population combination
  fi

  grm_file="${grmDir}/${pop}.${sex_dir}_GRM"
  infile="${indir}/${pop}.${sex_dir}.genotype"

  out_dir="${results_dir}/${pop}/${sex}"
  mkdir -p "${out_dir}/log"

  task_list=()
  for pheno_file in ${pheno_files}; do
    pheno_basename=$(basename "$pheno_file")
    # Extract the date and sample size from the phenotype filename
    date_sample=$(echo "$pheno_basename" | grep -oP '\d{8}_\w+_\w_\d+')
    
    num_samples=$(echo "$date_sample" | grep -oP '\d+$')
    scratch_dir="${scratch_base}/${pop}_${sex}_${RANDOM}"
    mkdir -p "${scratch_dir}"
    
    # Find the covariate files
    covar_file=$(find "${covar_dir}/Discrete/${pop}/${sex}/" -type f -name "covar_*_${date_sample}.txt" -not -path "*/archive/*")
    
    # Find the correct quantitative covariate file
    if [[ "$pheno_basename" == *"_smri_vol_scs_wholeb_ROC0_2.txt" ]]; then
      qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -type f -name "qcovar_noICV_*_${date_sample}.txt" -not -path "*/archive/*")
    else
      qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -type f -name "qcovar_*_${date_sample}.txt" -not -path "*/archive/*" -not -name "*noICV*")
    fi

    task_list+=("${pheno_file} ${grm_file} ${out_dir} ${date} ${covar_file} ${qcovar_file} ${num_samples} ${scratch_dir} ${pop} ${sex} ${infile}")
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