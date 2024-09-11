#!/bin/bash
#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l h_rt=15:00:00,h_data=5G,highp,arch=intel-gold*
#$ -pe shared 32
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/EUR_Male_MLMA_GWAS_$JOB_ID.$TASK_ID.out
#$ -e /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/EUR_Male_MLMA_GWAS_$JOB_ID.$TASK_ID.err
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-17:1

# Software path
gcta="/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1"

# Get current date
date=$(date +"%m%d%Y")

# Define constants
pop="EUR"
sex="M"

# Define directories
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
results_dir="${base_dir}/Results/${pop}/${sex}"
indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/males"
grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_males"

# Create results directory if it doesn't exist
mkdir -p "${results_dir}"

# Define log file
log_file="${results_dir}/GCTA_MLMA_${JOB_ID}_${SGE_TASK_ID}.log"

# Function to log messages
log_message() {
    echo "$(date): $1" >> "${log_file}"
}

# Extract all phenotype names
readarray -t phenotypes < <(find "${base_dir}/Phenotypes/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_smri_vol_*.txt" -not -path "*/archive/*" | xargs -I{} basename {} | grep -oP 'smri_vol_scs_\w+_ROC0_2' | sort)

# Get the current phenotype based on the task ID
current_phenotype=${phenotypes[$SGE_TASK_ID-1]}

log_message "Processing phenotype: ${current_phenotype}"

# Define input files
infile="${indir}/${pop}.males.genotype"
grm_file="${grmDir}/${pop}.males_GRM"
pheno_file=$(find "${base_dir}/Phenotypes/${pop}/${sex}/" -type f -name "*${current_phenotype}*.txt" -not -path "*/archive/*")
covar_file=$(find "${base_dir}/Covariates/Discrete/${pop}/${sex}/" -type f -name "covar_*${pop}_${sex}*.txt" -not -path "*/archive/*")
qcovar_prefix="qcovar_"
[ "$current_phenotype" = "smri_vol_scs_wholeb_ROC0_2" ] && qcovar_prefix="qcovar_noICV_"
qcovar_file=$(find "${base_dir}/Covariates/Quantitative/${pop}/${sex}/" -type f -name "${qcovar_prefix}*${pop}_${sex}*.txt" -not -path "*/archive/*")

# Log files found for debugging
log_message "Phenotype file found: ${pheno_file}"
log_message "Covariate file found: ${covar_file}"
log_message "Quantitative covariate file found: ${qcovar_file}"

# Create a unique scratch directory for this job
scratch_dir="/u/scratch/c/cobeaman/GCTA_GWAS_${JOB_ID}_${SGE_TASK_ID}"
mkdir -p "${scratch_dir}"
log_message "Created scratch directory: ${scratch_dir}"

# Transfer input files to scratch
log_message "Transferring input files to scratch..."
cp "${infile}.bed" "${infile}.bim" "${infile}.fam" "${scratch_dir}/"
cp "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" "${scratch_dir}/"
log_message "Transfer completed."

# Update file paths to use scratch
scratch_infile="${scratch_dir}/$(basename ${infile})"
scratch_grm_file="${scratch_dir}/$(basename ${grm_file})"

# Output file
num_samples=$(wc -l < "${pheno_file}")
out_file="${results_dir}/${current_phenotype}_${pop}_${sex}_n${num_samples}_${date}_${JOB_ID}_${SGE_TASK_ID}"

# Check if all necessary files exist in scratch
required_files=(
    "${pheno_file}" "${covar_file}" "${qcovar_file}"
    "${scratch_grm_file}.grm.bin" "${scratch_grm_file}.grm.id" "${scratch_grm_file}.grm.N.bin"
    "${scratch_infile}.bed" "${scratch_infile}.bim" "${scratch_infile}.fam"
)
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        log_message "Required file $file does not exist. Exiting."
        exit 1
    else
       log_message "Found required file: $file"
    fi
done
log_message "All required files found. Starting GCTA --mlma analysis."

# Run GCTA MLMA
$gcta --mlma \
      --bfile "${scratch_infile}" \
      --grm "${scratch_grm_file}" \
      --pheno "${pheno_file}" \
      --covar "${covar_file}" \
      --qcovar "${qcovar_file}" \
      --thread-num 32 \
      --out "${out_file}"

# Error checker
if [ $? -ne 0 ]; then
  log_message "Error running GCTA MLMA for ${pop} ${sex} ${current_phenotype}"
  exit 1
fi
log_message "Completed GCTA MLMA for ${pop} ${sex} ${current_phenotype} with ${num_samples} samples"

# Clean up scratch directory
log_message "Cleaning up scratch directory..."
rm -rf "${scratch_dir}"
log_message "Scratch directory cleaned up."

log_message "Analysis complete for phenotype: ${current_phenotype}"