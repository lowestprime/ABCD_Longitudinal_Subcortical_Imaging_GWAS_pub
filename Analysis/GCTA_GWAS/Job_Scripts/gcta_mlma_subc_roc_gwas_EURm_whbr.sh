#!/bin/sh
# Perform Genome-Wide association analysis using GCTA - Test Run - EUR M smri_vol_scs_wholeb_ROC0_2

# Set working directory
#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
# Request resources
#$ -l h_rt=30:00:00,h_data=5G,highp,arch=intel-gold*
#$ -pe shared 32
# Output and error notification preferences
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run/GCTA_MLMA_$JOB_ID.log
#$ -j y # join std error and std output streams, yes
# Email notifications
#$ -M $USER@mail # do not change, can not use custom address
#$ -m bea # email when job begins, ends, and if it aborts

# Define constants
# Current date
date=$(date +"%m%d%Y")
# Population, sex, and phenotype
pop="EUR"
sex="M"
phenotype="smri_vol_scs_wholeb_ROC0_2"

# Define directories
# Base and results
base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
results_dir="${base_dir}/Results/test_run"
mkdir -p "${results_dir}"
# GRM and bfile
indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/males"
grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_males"

# Scratch directory for temporary accelerated I/O
scratch_dir="/u/scratch/c/cobeaman/GCTA_GWAS_${JOB_ID}"
mkdir -p "${scratch_dir}"

# Software path
gcta="/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1"

# Phenotype and covariate files
pheno_file=$(find "${base_dir}/Phenotypes/${pop}/${sex}/" -type f -name "*${phenotype}*.txt" -not -path "*/archive/*")
covar_file=$(find "${base_dir}/Covariates/Discrete/${pop}/${sex}/" -type f -name "covar_*${pop}_${sex}*.txt" -not -path "*/archive/*")
qcovar_prefix="qcovar_"
[ "$phenotype" = "smri_vol_scs_wholeb_ROC0_2" ] && qcovar_prefix="qcovar_noICV_"
qcovar_file=$(find "${base_dir}/Covariates/Quantitative/${pop}/${sex}/" -type f -name "${qcovar_prefix}*${pop}_${sex}*.txt" -not -path "*/archive/*")

# Log files found for debugging
echo "Phenotype file found: ${pheno_file}" >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log
echo "Covariate file found: ${covar_file}" >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log
echo "Quantitative covariate file found: ${qcovar_file}" >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log

# GRM and bfile
infile="${indir}/${pop}.males.genotype"
grm_file="${grmDir}/${pop}.males_GRM"

# Transfer GRM and bfile to the scratch directory
echo "Transferring GRM and bfile to scratch directory..." >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log
cp "${infile}.bed" "${infile}.bim" "${infile}.fam" "${scratch_dir}/"
cp "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" "${scratch_dir}/"
echo "Transfer completed." >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log

# Update file paths to use scratch directory
scratch_infile="${scratch_dir}/$(basename ${infile})"
scratch_grm_file="${scratch_dir}/$(basename ${grm_file})"

# Output file
num_samples=$(wc -l < "${pheno_file}")
out_file="${results_dir}/${phenotype}_${pop}_${sex}_n${num_samples}_${date}_$JOB_ID"

# Check if all necessary files exist in scratch directory
required_files=(
    "${pheno_file}" "${covar_file}" "${qcovar_file}"
    "${scratch_grm_file}.grm.bin" "${scratch_grm_file}.grm.id" "${scratch_grm_file}.grm.N.bin"
    "${scratch_infile}.bed" "${scratch_infile}.bim" "${scratch_infile}.fam"
)
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Required file $file does not exist. Exiting." >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log
        exit 1
    else
       echo "Found required file: $file" >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log
    fi
done
echo "All required files found in scratch directory. Starting GCTA --mlma analysis." >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log

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
  echo "Error running GCTA MLMA for ${pop} ${sex} ${phenotype}" >&2
  exit 1
fi
echo "Completed GCTA MLMA for ${pop} ${sex} ${phenotype} with ${num_samples} samples" >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log

# Clean up scratch directory (optional)
echo "Cleaning up scratch directory..." >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log
rm -rf "${scratch_dir}"
echo "Scratch directory cleaned up." >> ${results_dir}/GCTA_MLMA_${JOB_ID}.log