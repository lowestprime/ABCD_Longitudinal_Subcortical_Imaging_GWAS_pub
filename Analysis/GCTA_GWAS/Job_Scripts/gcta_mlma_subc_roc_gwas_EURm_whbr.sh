#!/bin/sh
# Perform Genome-Wide association analysis using GCTA - Test Run - EUR M smri_vol_scs_wholeb_ROC0_2

# Set working directory
#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data

# Request resources
#$ -l highp,h_rt=70:00:00,h_data=5G
#$ -pe shared 32
#$ -l arch=intel-gold*  # Requesting nodes with Intel Gold CPUs

# Output and error handling
#$ -o GCTA_GWAS_EURM_WhBr${JOB_ID}.out
#$ -j y # join std error and std output streams, yes

# Email notifications
#$ -M cobeaman@g.ucla.edu
#$ -m bea # email when job begins, ends, and if it aborts

# Define constants
# Current date
date=$(date +"%m%d%Y")
# Population, sex, and phenotype
pop="EUR"
sex="M"
phenotype="smri_vol_scs_wholeb_ROC0_2"

# Directories
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/males"
grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_males"
results_dir="${base_dir}/Results/test_run"
mkdir -p "${results_dir}"

# Software path
gcta="/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1"

# Phenotype and covariate files
# Find phenotype file
pheno_file=$(find "${base_dir}/Phenotypes/${pop}/${sex}/" -type f -name "*${phenotype}*.txt")

# Find covar file
covar_file=$(find "${base_dir}/Covariates/Discrete/${pop}/${sex}/" -type f -name "covar_*${pop}_${sex}*.txt")

# Determine qcovar prefix
if [ "$phenotype" = "smri_vol_scs_wholeb_ROC0_2" ]; then
  qcovar_prefix="qcovar_noICV_"
else
  qcovar_prefix="qcovar_"
fi
# Find qcovar file
qcovar_file=$(find "${base_dir}/Covariates/Quantitative/${pop}/${sex}/" -type f -name "${qcovar_prefix}*${pop}_${sex}*.txt")

# GRM and bfile
infile="${indir}/${pop}.males.genotype"
grm_file="${grmDir}/${pop}.males_GRM"

# Output file
num_samples=$(wc -l < "${pheno_file}")
out_file="${results_dir}/${phenotype}_${date}_n${num_samples}"

# Check if all necessary files exist
required_files=(
    "${pheno_file}"
    "${covar_file}"
    "${qcovar_file}"
    "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin"
    "${infile}.bed" "${infile}.bim" "${infile}.fam"
)

# Loop through each required file and check if it exists
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Required file $file does not exist. Exiting."
        exit 1
    else
        echo "Found required file: $file"
    fi
done

echo "All required files found. Starting GCTA --mlma analysis."

# Run GCTA MLMA
$gcta --mlma \
      --bfile "${infile}" \
      --grm "${grm_file}" \
      --pheno "${pheno_file}" \
      --covar "${covar_file}" \
      --qcovar "${qcovar_file}" \
      --thread-num 32 \
      --out "${out_file}"

if [ $? -ne 0 ]; then
  echo "Error running GCTA MLMA for ${phenotype}" >&2
  exit 1
fi

echo "Completed GCTA MLMA for ${pop} ${sex} ${phenotype} with ${num_samples} samples"
