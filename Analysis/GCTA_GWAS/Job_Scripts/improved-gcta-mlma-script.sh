#!/bin/bash

# Perform Genome-Wide association analysis using GCTA --mlma for all jobs in parallel

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=70:00:00,h_data=8G
#$ -pe shared 16
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run/GCTA_GWAS_$JOB_ID.out
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-3:1

set -euo pipefail

# Load modules
. /u/local/Modules/default/init/modules.sh
module load parallel

# Check if GNU Parallel is available
if ! command -v parallel &> /dev/null; then
    echo "Error: GNU Parallel is not installed or not functioning correctly."
    exit 1
fi

echo "GNU Parallel is working correctly. Proceeding with the analysis..."

# Define variables
gcta="/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1"
date=$(date +"%m%d%Y")
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
pheno_dir="${base_dir}/Phenotypes"
covar_dir="${base_dir}/Covariates"
results_dir="${base_dir}/Results"
scratch_base="/u/scratch/c/cobeaman/GCTA_GWAS_${JOB_ID}"

mkdir -p "${scratch_base}"

# Define populations and sexes
pops=("EUR" "AMR" "AFR")
pop=${pops[$SGE_TASK_ID-1]}
sexes=("F" "M")

# Specify combinations to skip
declare -A skip_combinations=(
    ["EUR_M_smri_vol_scs_wholeb_ROC0_2"]=1
    ["AMR_F_smri_vol_scs_wholeb_ROC0_2"]=1
)

# Function to check if a combination should be skipped
should_skip() {
    local combination="${1}_${2}_${3}"
    [[ -n "${skip_combinations[$combination]:-}" ]]
}

# Function to check if all required files exist
check_files_exist() {
    local pop=$1 sex=$2 pheno_name=$3

    local pheno_file grm_file infile covar_file qcovar_file
    pheno_file=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_${pop}_${sex}_*_${pheno_name}.txt" -not -path "*/archive/*")
    grm_file="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_females/${pop}.females_GRM"
    infile="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/females/${pop}.females.genotype"
    covar_file="${covar_dir}/Discrete/${pop}/${sex}/covar_574229_09052024_${pop}_${sex}_1916.txt"
    qcovar_file="${covar_dir}/Quantitative/${pop}/${sex}/qcovar_574229_09052024_${pop}_${sex}_1916.txt"

    # Handle noICV case
    [[ "$pheno_name" == "smri_vol_scs_wholeb_ROC0_2" ]] && qcovar_file="${covar_dir}/Quantitative/${pop}/${sex}/qcovar_noICV_574229_09052024_${pop}_${sex}_1916.txt"

    local missing_files=()
    for file in "$pheno_file" "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" \
                "${infile}.bed" "${infile}.bim" "${infile}.fam" "$covar_file" "$qcovar_file"; do
        [[ ! -f "$file" ]] && missing_files+=("$file")
    done

    if (( ${#missing_files[@]} > 0 )); then
        echo "Error: The following required files are missing:"
        printf '%s\n' "${missing_files[@]}"
        return 1
    fi

    return 0
}

# Function to run GCTA MLMA
run_gcta_mlma() {
    local pop=$1 sex=$2 pheno_name=$3 scratch_dir=$4

    if should_skip "$pop" "$sex" "$pheno_name"; then
        echo "Skipping $pop $sex $pheno_name"
        return 0
    fi

    if ! check_files_exist "$pop" "$sex" "$pheno_name"; then
        echo "Skipping analysis for $pop $sex $pheno_name due to missing files"
        return 1
    fi

    local pheno_file grm_file infile covar_file qcovar_file out_file num_samples
    pheno_file=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_${pheno_name}.txt" -not -path "*/archive/*")
    grm_file="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_females/${pop}.females_GRM"
    infile="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/females/${pop}.females.genotype"
    covar_file="${covar_dir}/Discrete/${pop}/${sex}/covar_574229_09052024_${pop}_${sex}_1916.txt"
    qcovar_file="${covar_dir}/Quantitative/${pop}/${sex}/qcovar_574229_09052024_${pop}_${sex}_1916.txt"

    [[ "$pheno_name" == "smri_vol_scs_wholeb_ROC0_2" ]] && qcovar_file="${covar_dir}/Quantitative/${pop}/${sex}/qcovar_noICV_574229_09052024_${pop}_${sex}_1916.txt"

    num_samples=$(wc -l < "$pheno_file")
    out_file="${results_dir}/${pop}/${sex}/${pheno_name}_${pop}_${sex}_n${num_samples}_${date}_${JOB_ID}"

    local scratch_bfile="${scratch_dir}/$(basename "${infile}")"
    local scratch_grm="${scratch_dir}/$(basename "${grm_file}")"
    cp "${infile}.bed" "${infile}.bim" "${infile}.fam" "${scratch_dir}/"
    cp "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" "${scratch_dir}/"

    "$gcta" --mlma \
        --bfile "${scratch_bfile}" \
        --grm "${scratch_grm}" \
        --pheno "${pheno_file}" \
        --covar "${covar_file}" \
        --qcovar "${qcovar_file}" \
        --thread-num 16 \
        --out "${out_file}" > "${results_dir}/${pop}/${sex}/log/${pheno_name}_${pop}_${sex}.log" 2>&1

    if (( $? != 0 )); then
        echo "Error running GCTA MLMA for ${pheno_name}" >&2
        return 1
    fi

    rm -f "${scratch_dir}"/*
}

export -f run_gcta_mlma check_files_exist should_skip

# Main loop
for pop in "${pops[@]}"; do
    echo "Current population: $pop"
    for sex in "${sexes[@]}"; do
        echo "Preparing tasks for Population: $pop, Sex: $sex"
        
        out_dir="${results_dir}/${pop}/${sex}"
        mkdir -p "${out_dir}/log"

        pheno_files=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_smri_vol_*.txt" -not -path "*/archive/*")
        if [[ -z "$pheno_files" ]]; then
            echo "Error: No phenotype files found for ${pop} ${sex}. Skipping..."
            continue
        fi

        task_list=()
        while IFS= read -r pheno_file; do
            pheno_name=$(basename "$pheno_file" | grep -oP 'smri_vol_scs_\w+_ROC0_2')
            scratch_dir="${scratch_base}/${pop}_${sex}_${RANDOM}"
            mkdir -p "${scratch_dir}"
            task_list+=("$pop" "$sex" "$pheno_name" "$scratch_dir")
        done <<< "$pheno_files"

        parallel --jobs 16 run_gcta_mlma ::: "${task_list[@]}"
    done
done

rm -rf "${scratch_base}"
echo "Analysis complete."