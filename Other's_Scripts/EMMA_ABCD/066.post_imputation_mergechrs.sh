#!/bin/sh
# merge all chr files

#$ -wd /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/06.post_imputation
#$ -l highp,h_rt=48:00:00,h_data=4G
#$ -pe shared 16
#$ -o joblog/120422.$JOB_ID.068.post_imputation_mergechrs.log
#$ -j y
#$ -M $USER@mail
#$ -m bea

# directory
post_imputation_dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/06.post_imputation

# software
. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624

cohort="ABCDr4"
ref="MikeDB_db155"
filebase="${cohort}_TopMed_chr1.22_rsID"
filebase=$filebase.$ref

cd $post_imputation_dir

# merge all chr files
touch merge_list.txt
for i in {1..22}; do
    echo chr${i}.imp.flt.r208_RSidOnly.$ref >> merge_list.txt
done

plink \
    --merge-list merge_list.txt \
    --allow-no-sex \
    --threads 16 \
    --memory 64000 \
    --keep-allele-order \
    --out $filebase

