#!/bin/bash

#$ -wd /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data
#$ -l highp,h_rt=70:00:00,h_data=8G
#$ -pe shared 16
#$ -o /u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run/GCTA_GWAS_$JOB_ID.out
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-3:1

. /u/local/Modules/default/init/modules.sh
module load parallel

if ! command -v parallel &> /dev/null; then
  echo "Error: GNU Parallel is not installed or not functioning correctly."
  exit 1
fi

echo "GNU Parallel is working correctly. Proceeding with the analysis..."

gcta=/u/project/lhernand/sganesh/apps/gcta/gcta-1.94.1
date=$(date +"%m%d%Y")

base_dir="/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data"
pheno_dir="${base_dir}/Phenotypes"
covar_dir="${base_dir}/Covariates"
results_dir="${base_dir}/Results"
scratch_base="/u/scratch/c/cobeaman/GCTA_GWAS_${JOB_ID}"
mkdir -p "${scratch_base}"

pops=("EUR" "AMR" "AFR")
pop=${pops[$SGE_TASK_ID-1]}
sexes=("F" "M")

skip_combinations=("EUR_M_smri_vol_scs_wholeb_ROC0_2" "AMR_F_smri_vol_scs_wholeb_ROC0_2")

should_skip() {
  local combination="${1}_${2}_${3}"
  [[ " ${skip_combinations[*]} " =~ " ${combination} " ]] && return 0 || return 1
}

run_gcta_mlma() {
    local pop=$1
    local sex=$2
    local pheno_name=$3
    local date_sample=$4
    local num_samples=$5

    if should_skip "$pop" "$sex" "$pheno_name"; then
        echo "Skipping $pop $sex $pheno_name"
        return 0
    fi

    local sex_dir=$([ "$sex" = "M" ] && echo "males" || echo "females")
    local indir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/TOPMed_imputed/splitted_by_ancestry_groups/${sex_dir}"
    local grmDir="/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_5/genotype/GRM/grm_${sex_dir}"

    local pheno_file="${pheno_dir}/${pop}/${sex}/${date_sample}_pheno_${pop}_${sex}_${num_samples}_${pheno_name}.txt"
    local grm_file="${grmDir}/${pop}.${sex_dir}_GRM"
    local infile="${indir}/${pop}.${sex_dir}.genotype"
    local covar_file="${covar_dir}/Discrete/${pop}/${sex}/covar_${date_sample}.txt"
    local qcovar_file="${covar_dir}/Quantitative/${pop}/${sex}/qcovar_${date_sample}.txt"
    [[ "$pheno_name" == "smri_vol_scs_wholeb_ROC0_2" ]] && qcovar_file="${covar_dir}/Quantitative/${pop}/${sex}/qcovar_noICV_${date_sample}.txt"

    local out_file="${results_dir}/${pop}/${sex}/${pheno_name}_${pop}_${sex}_n${num_samples}_${date}_${JOB_ID}"
    local scratch_dir="${scratch_base}/${pop}_${sex}_${RANDOM}"
    mkdir -p "${scratch_dir}"

    local scratch_bfile="${scratch_dir}/$(basename ${infile})"
    local scratch_grm="${scratch_dir}/$(basename ${grm_file})"

    for file in "$pheno_file" "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" \
                "${infile}.bed" "${infile}.bim" "${infile}.fam" "$covar_file" "$qcovar_file"; do
        if [ ! -f "$file" ]; then
            echo "Missing file: $file"
            return 1
        fi
    done

    cp "${infile}.bed" "${infile}.bim" "${infile}.fam" "${scratch_dir}/"
    cp "${grm_file}.grm.bin" "${grm_file}.grm.id" "${grm_file}.grm.N.bin" "${scratch_dir}/"

    $gcta --mlma \
          --bfile "${scratch_bfile}" \
          --grm "${scratch_grm}" \
          --pheno "${pheno_file}" \
          --covar "${covar_file}" \
          --qcovar "${qcovar_file}" \
          --thread-num 16 \
          --out "${out_file}" > "${results_dir}/${pop}/${sex}/log/${pheno_name}_${pop}_${sex}.log" 2>&1

    rm -rf "${scratch_dir}"
}

export -f run_gcta_mlma should_skip

for sex in "${sexes[@]}"; do
    echo "Processing: $pop $sex"
    
    pheno_files=$(find "${pheno_dir}/${pop}/${sex}/" -maxdepth 1 -type f -name "*_pheno_*_${pop}_${sex}_*_smri_vol_*.txt" -not -path "*/archive/*")

    task_list=()
    for pheno_file in $pheno_files; do
        pheno_basename=$(basename "$pheno_file")
        date_sample=$(echo "$pheno_basename" | grep -oP '\d{8}_\w+_\w_\d+')
        num_samples=$(echo "$date_sample" | grep -oP '\d+$')
        pheno_name=$(echo "$pheno_basename" | grep -oP 'smri_vol_scs_\w+_ROC0_2')

        task_list+=("$pop $sex $pheno_name $date_sample $num_samples")
    done

    parallel --jobs 16 run_gcta_mlma ::: "${task_list[@]}"
done

rm -rf "${scratch_base}"
echo "Analysis complete."