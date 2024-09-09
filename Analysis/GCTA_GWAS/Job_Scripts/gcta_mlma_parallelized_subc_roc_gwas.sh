#!/bin/sh
# Perform Genome-Wide association analysis using GCTA MLMA for all jobs in parallel

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=70:00:00,h_data=8G
#$ -pe shared 16
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run/MLMA_GWAS_$JOB_ID.out
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-3:1

# Load the GNU Parallel module
. /u/local/Modules/default/init/modules.sh
module load parallel

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
    local pop=$1
    local sex=$2
    local pheno_name=$3

    # Ensure variables are properly expanded and paths are correct
    local pheno_file=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_${pheno_name}.txt" -not -path "*/archive/*")

    if [[ -z "$pheno_file" ]]; then
        echo "Error: Phenotype file not found for ${pop} ${sex} ${pheno_name}"
        return 1
    fi

    # Define other file paths
    local sex_dir=$([ "$sex" = "F" ] && echo "females" || echo "males")
    local grm_file="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}/${pop}.${sex_dir}_GRM"
    local infile="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}/${pop}.${sex_dir}.genotype"
    local covar_file=$(find "${covar_dir}/Discrete/${pop}/${sex}/" -maxdepth 1 -type f -name "covar_*_${pop}_${sex}_*.txt" -not -path "*/archive/*")

    # Condensed qcovar file logic
    local qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -maxdepth 1 -type f -name "$( [[ "$pheno_name" == "smri_vol_scs_wholeb_ROC0_2" ]] && echo "qcovar_noICV_*" || echo "qcovar_*")_${pop}_${sex}_*.txt" -not -path "*/archive/*")

    # Collect missing files
    local missing_files=()
    for file in "$pheno_file" "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" \
                "${infile}.bed" "${infile}.bim" "${infile}.fam" "$covar_file" "$qcovar_file"; do
        if [ ! -f "$file" ]; then
            missing_files+=("$file")
        fi
    done

    if [ ${#missing_files[@]} -ne 0 ]; then
        echo "Error: The following required files are missing:"
        printf '%s\n' "${missing_files[@]}"
        return 1
    fi

    return 0
}

# Function to run GCTA MLMA
run_gcta_mlma() {
    local pop=$1
    local sex=$2
    local pheno_name=$3
    local scratch_dir=$4

    # Skip specified combinations
    if should_skip "$pop" "$sex" "$pheno_name"; then
        echo "Skipping $pop $sex $pheno_name"
        return 0
    fi

    # Check if all required files exist
    if ! check_files_exist "$pop" "$sex" "$pheno_name"; then
        echo "Skipping analysis for $pop $sex $pheno_name due to missing files"
        return 1
    fi

    # Get phenotype file path
    local pheno_file=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_${pheno_name}.txt" -not -path "*/archive/*")

    # Define other file paths
    local sex_dir=$([ "$sex" = "F" ] && echo "females" || echo "males")
    local grm_file="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}/${pop}.${sex_dir}_GRM"
    local infile="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}/${pop}.${sex_dir}.genotype"
    local covar_file=$(find "${covar_dir}/Discrete/${pop}/${sex}/" -maxdepth 1 -type f -name "covar_*_${pop}_${sex}_*.txt" -not -path "*/archive/*")

    # Condensed qcovar file logic
    local qcovar_file=$(find "${covar_dir}/Quantitative/${pop}/${sex}/" -maxdepth 1 -type f -name "$( [[ "$pheno_name" == "smri_vol_scs_wholeb_ROC0_2" ]] && echo "qcovar_noICV_*" || echo "qcovar_*")_${pop}_${sex}_*.txt" -not -path "*/archive/*")

    # Calculate number of samples from phenotype file
    local num_samples=$(wc -l < "$pheno_file")
    
    # Output file path
    local out_file="${results_dir}/${pop}/${sex}/${pheno_name}_${pop}_${sex}_n${num_samples}_${date}_${JOB_ID}"

    # Copy GRM and genotype files to scratch directory
    local scratch_bfile="${scratch_dir}/$(basename "${infile}")"
    local scratch_grm="${scratch_dir}/$(basename "${grm_file}")"
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
          --out "${out_file}" > "${results_dir}/${pop}/${sex}/log/${pheno_name}_${pop}_${sex}.log" 2>&1

    # Check for errors
    if [ $? -ne 0 ]; then
        echo "Error running GCTA MLMA for ${pheno_name}" >&2
        return 1
    fi

    # Clean up scratch directory
    rm -f "${scratch_dir}"/*
}

export -f run_gcta_mlma check_files_exist should_skip

# Main loop
for pop in "${pops[@]}"; do
    echo "Current population: $pop"
    for sex in "${sexes[@]}"; do
        echo "Preparing tasks for Population: $pop, Sex: $sex"
        
        sex_dir=$([ "$sex" = "F" ] && echo "females" || echo "males")
        indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}"
        grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}"
        
        out_dir="${results_dir}/${pop}/${sex}"
        mkdir -p "${out_dir}/log"

        echo "Phenotype directory: ${pheno_dir}/${pop}/${sex}"
        echo "Discrete Covariate directory: ${covar_dir}/Discrete/${pop}/${sex}"
        echo "Quantitative Covariate directory: ${covar_dir}/Quantitative/${pop}/${sex}"
        echo "bfile directory: ${indir}"
        echo "GRM directory: ${grmDir}"

        # Extract all phenotype names
        pheno_names=($(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_smri_vol_*.txt" -not -path "*/archive/*" | xargs -I{} basename {} | grep -oP 'smri_vol_scs_\w+_ROC0_2'))

        if [ ${#pheno_names[@]} -eq 0 ]; then
            echo "Error: No phenotype files found for ${pop} ${sex}. Skipping..."
            continue
        fi

        task_list=()
        for pheno_name in "${pheno_names[@]}"; do
            scratch_dir="${scratch_base}/${pop}_${sex}_${RANDOM}"
            mkdir -p "${scratch_dir}"

            task_list+=("\"$pop\" \"$sex\" \"$pheno_name\" \"$scratch_dir\"")
        done

        # Run tasks in parallel with proper quoting for arguments
        parallel --jobs 16 run_gcta_mlma ::: "${task_list[@]}"
    done
done

# Clean up scratch directory
rm -rf "${scratch_base}"
echo "Analysis complete."
