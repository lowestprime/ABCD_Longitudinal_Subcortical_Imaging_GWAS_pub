#!/bin/bash
# Perform Genome-Wide association analysis using GCTA MLMA

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l h_rt=20:00:00,h_data=15G,highp,arch=intel-gold*
#$ -pe shared 32
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run/MLMA_GWAS_$JOB_ID.$TASK_ID.out
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-50:1

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
sexes=("F" "M")

# Function to check if all required files exist
check_files_exist() {
    local pop=$1
    local sex=$2
    local pheno_name=$3

    local pheno_file=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_${pheno_name}.txt" -not -path "*/archive/*" | head -n 1)

    if [[ -z "$pheno_file" ]]; then
        echo "Error: Phenotype file not found for ${pop} ${sex} ${pheno_name}"
        return 1
    fi

    local sex_dir="$sex"
    local grm_file="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}/${pop}.${sex_dir}_GRM"
    local infile="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}/${pop}.${sex_dir}.genotype"
    local covar_file=$(find "${covar_dir}/Discrete/${pop}/${sex}/" -maxdepth 1 -type f -name "covar_*_${pop}_${sex}_*.txt" -not -path "*/archive/*" | head -n 1)
    local qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -maxdepth 1 -type f -name "$( [[ "$pheno_name" == "smri_vol_scs_wholeb_ROC0_2" ]] && echo "qcovar_noICV_*" || echo "qcovar_*")_${pop}_${sex}_*.txt" -not -path "*/archive/*" | head -n 1)

    local missing_files=()
    for file in "$pheno_file" "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" \
                "${infile}.bed" "${infile}.bim" "${infile}.fam" "$covar_file" "$qcovar_file"; do
        if [ ! -f "$file" ]; then
            missing_files+=("$file")
        fi
    done

    if [ ${#missing_files[@]} -ne 0 ]; then
        echo "Error: The following required files are missing for ${pop} ${sex} ${pheno_name}:"
        printf '%s\n' "${missing_files[@]}"
        return 1
    fi

    return 0
}

# Function to copy files with locking
copy_with_lock() {
    local src=$1
    local dest=$2
    local lock_file="${dest}.lock"

    (
        flock -w 60 9 || exit 1
        if [ ! -f "$dest" ]; then
            cp "$src" "$dest"
        fi
    ) 9>"$lock_file"
}

# Function to run GCTA MLMA
run_gcta_mlma() {
    local pop=$1
    local sex=$2
    local pheno_name=$3
    local scratch_dir=$4

    echo "Running GCTA MLMA for Population: $pop, Sex: $sex, Phenotype: $pheno_name"
    echo "Using scratch directory: $scratch_dir"

    # Check if all required files exist
    if ! check_files_exist "$pop" "$sex" "$pheno_name"; then
        echo "Skipping analysis for $pop $sex $pheno_name due to missing files"
        return 1
    fi

    local pheno_file=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_${pheno_name}.txt" -not -path "*/archive/*" | head -n 1)
    local sex_dir="$sex"
    local grm_file="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}/${pop}.${sex_dir}_GRM"
    local infile="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}/${pop}.${sex_dir}.genotype"
    local covar_file=$(find "${covar_dir}/Discrete/${pop}/${sex}/" -maxdepth 1 -type f -name "covar_*_${pop}_${sex}_*.txt" -not -path "*/archive/*" | head -n 1)
    local qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -maxdepth 1 -type f -name "$( [[ "$pheno_name" == "smri_vol_scs_wholeb_ROC0_2" ]] && echo "qcovar_noICV_*" || echo "qcovar_*")_${pop}_${sex}_*.txt" -not -path "*/archive/*" | head -n 1)

    local num_samples=$(wc -l < "$pheno_file")
    local out_file="${results_dir}/${pop}/${sex}/${pheno_name}_${pop}_${sex}_n${num_samples}_${date}_${JOB_ID}"

    local scratch_bfile="${scratch_dir}/$(basename "${infile}")"
    local scratch_grm="${scratch_dir}/$(basename "${grm_file}")"

    # Copy files to scratch directory with locking
    copy_with_lock "${infile}.bed" "${scratch_bfile}.bed"
    copy_with_lock "${infile}.bim" "${scratch_bfile}.bim"
    copy_with_lock "${infile}.fam" "${scratch_bfile}.fam"
    copy_with_lock "${grm_file}.grm.bin" "${scratch_grm}.grm.bin"
    copy_with_lock "${grm_file}.grm.id" "${scratch_grm}.grm.id"
    copy_with_lock "${grm_file}.grm.N.bin" "${scratch_grm}.grm.N.bin"

    local log_file="${results_dir}/${pop}/${sex}/log/${pheno_name}_${pop}_${sex}.log"

    if ! $gcta --mlma \
          --bfile "${scratch_bfile}" \
          --grm "${scratch_grm}" \
          --pheno "${pheno_file}" \
          --covar "${covar_file}" \
          --qcovar "${qcovar_file}" \
          --thread-num 32 \
          --out "${out_file}" > "${log_file}" 2>&1
    then
        echo "Error running GCTA MLMA for ${pheno_name}. Check ${log_file} for details." >&2
        return 1
    fi

    echo "GCTA MLMA completed successfully for ${pheno_name}"
    return 0
}

# Generate list of jobs, prioritizing EUR and removing duplicates
job_list=()
processed_jobs=()
for pop in "EUR" "${pops[@]}"; do
    for sex in "${sexes[@]}"; do
        pheno_names=($(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_smri_vol_*.txt" -not -path "*/archive/*" | xargs -I{} basename {} | grep -oP 'smri_vol_scs_\w+_ROC0_2'))

        for pheno_name in "${pheno_names[@]}"; do
            job_key="${pop}_${sex}_${pheno_name}"
            if [[ ! " ${processed_jobs[@]} " =~ " ${job_key} " ]]; then
                # Check if .mlma file already exists
                if ! ls "${results_dir}/${pop}/${sex}/${pheno_name}_${pop}_${sex}_"*.mlma 1> /dev/null 2>&1; then
                    job_list+=("$pop $sex $pheno_name")
                fi
                processed_jobs+=("$job_key")
            fi
        done
    done
done

echo "Total jobs to run: ${#job_list[@]}"

# Function to run a job
run_job() {
    local job_index=$1
    IFS=' ' read -r pop sex pheno_name <<< "${job_list[$job_index]}"
    scratch_dir="${scratch_base}/${pop}_${sex}"
    mkdir -p "${scratch_dir}"
    run_gcta_mlma "$pop" "$sex" "$pheno_name" "$scratch_dir"
}

# Run jobs with improved job management
total_jobs=${#job_list[@]}
current_job=0
max_concurrent_jobs=10

while [ $current_job -lt $total_jobs ]; do
    # Count running jobs
    running_jobs=$(jobs -p | wc -l)

    # Start new jobs if below the max concurrent limit
    while [ $running_jobs -lt $max_concurrent_jobs ] && [ $current_job -lt $total_jobs ]; do
        run_job $current_job &
        current_job=$((current_job + 1))
        running_jobs=$((running_jobs + 1))
    done

    # Wait for any job to finish before starting the next batch
    wait -n
done

# Wait for all remaining jobs to finish
wait

# Clean up main scratch directory
rm -rf "${scratch_base}"
echo "Analysis complete."
