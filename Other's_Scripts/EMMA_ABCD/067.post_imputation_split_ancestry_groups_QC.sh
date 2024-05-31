#!/bin/sh
# get allele frequency file and QC for each ancestry group

#$ -wd /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/06.post_imputation
#$ -l highp,h_rt=48:00:00,h_data=4G
#$ -pe shared 16
#$ -o joblog/120422.$JOB_ID.$TASK_ID.069.post_imputation_split_ancestry_groups_QC_MikeDB_db155.log
#$ -j y
#$ -M $USER@mail
#$ -m bea
#$ -t 1-5:1

cohort="ABCDr4"
ref="MikeDB_db155"
filebase="${cohort}_TopMed_chr1.22_rsID"
filebase=$filebase.$ref

pops=("EUR" "AMR" "AFR" "EAS" "SAS")
pop=${pops[$SGE_TASK_ID-1]}

# directories
ancestry_dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/03.ancestry_estimation/splitted_by_ancestry_groups
post_imputation_dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/06.post_imputation
ancestry_plink_dir=$post_imputation_dir/splitted_by_ancestry_groups
outdir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed
mkdir -p $ancestry_plink_dir/log $outdir/splitted_by_ancestry_groups/log

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

cd $post_imputation_dir

# get ancestry ID file and misc variables
filepath=`ls -1v $ancestry_dir/*${pop}*.rmsexmismatch*.fam | grep -v "allproblem"`
num_no_sex_mismatch=`wc -l $filepath | cut -f1 -d" "`
file=`basename $filepath | sed 's/\.fam//g'`
awk 'BEGIN {OFS = "\t"} {print $1, $2}' $filepath > $outdir/${cohort}_${pop}.${num_no_sex_mismatch}_no.sexmismatch_IDs.txt

ancestry_id_file="$outdir/${cohort}_${pop}.${num_no_sex_mismatch}_no.sexmismatch_IDs.txt"
if [ ! -e $ancestry_dir/sex.drop.ALLsamples*txt ]; then
    cat $ancestry_dir/sex.drop*txt > $ancestry_dir/sex.drop.ALLsamples.txt
    num=`wc -l $ancestry_dir/sex.drop.ALLsamples.txt | cut -f1 -d" "`
    mv $ancestry_dir/sex.drop.ALLsamples.txt $ancestry_dir/sex.drop.ALLsamples.${num}.txt
fi
num=`basename $ancestry_dir/sex.drop.ALLsamples*txt | grep -Eo "[0-9]+"`
mkdir -p $ancestry_plink_dir/log

# get MAF for each ancestry
$plink2 \
    --bfile $post_imputation_dir/$filebase \
    --keep $ancestry_id_file \
    --freq \
    --out $outdir/splitted_by_ancestry_groups/${filebase}_rm.sexmismatch.$num.$pop.$num_no_sex_mismatch

# QC
$plink2 \
    --bfile $post_imputation_dir/$filebase \
    --keep $ancestry_id_file \
    --maf 0.01 --geno 0.05 --hwe 1e-6 \
    --threads 16 \
    --memory 64000 \
    --make-bed --out $outdir/splitted_by_ancestry_groups/${filebase}_rm.sexmismatch.$num.$pop.$num_no_sex_mismatch.maf0.01.geno0.05.hwe1e6

