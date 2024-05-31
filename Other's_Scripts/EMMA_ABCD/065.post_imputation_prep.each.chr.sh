#!/bin/sh
# convert individual chr_vcf to plink
# format variant IDs and remove duplicate SNPs

#$ -wd /u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed/vcf_imputed_raw
#$ -o joblog/120322_job$JOB_ID.$TASK_ID.filter.merge.index.imputed-vcfs_re.log
#$ -j y
#$ -l h_rt=80:00:00,h_data=4G,highp
#$ -pe shared 16
#$ -M $USER@mail
#$ -m bea
#$ -t 1-22:1

# directories
chr=$SGE_TASK_ID
imputeFile="chr${chr}.imp.flt.r208"
TOPMed_dir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed
indir=$TOPMed_dir/vcf_imputed_filtered
vcf2plink_dir=$TOPMed_dir/vcf_imputed_filtered_plink
post_imputation_dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/06.post_imputation
mkdir -p $vcf2plink_dir/log $post_imputation_dir/log $post_imputation_dir/log_vcf2plinklog

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

cd $post_imputation_dir

# create chr:pos:ref:alt -> chr:pos conversion file
awk '{print $2}' $vcf2plink_dir/${imputeFile}.bim > ${imputeFile}_ChrPosAllele.txt
cut -f1,2 -d':' ${imputeFile}_ChrPosAllele.txt > ${imputeFile}_ChrPos.txt
paste ${imputeFile}_ChrPosAllele.txt ${imputeFile}_ChrPos.txt > ${imputeFile}_UpdateChrPos.txt && rm ${imputeFile}_ChrPosAllele.txt ${imputeFile}_ChrPos.txt

# update variant IDs chr:pos:ref:alt -> chr:pos
$plink2 \
    --bfile $vcf2plink_dir/${imputeFile} \
    --update-name ${imputeFile}_UpdateChrPos.txt \
    --threads 16 \
    --make-bed --out ${imputeFile}_UpdateChrPos

# extract & remove duplicated SNPs
cut -f2 ${imputeFile}_UpdateChrPos.bim | sort | uniq -d > ${imputeFile}_UpdateChrPos_duplicatedsnps.snplist.txt

# exclude dup SNPs
$plink2 \
    --bfile ${imputeFile}_UpdateChrPos \
    --exclude ${imputeFile}_UpdateChrPos_duplicatedsnps.snplist.txt \
    --threads 16 \
    --make-bed --out ${imputeFile}_UpdateChrPos_NoDuplicates

# update variant IDs chr:pos -> rsID
AllChr_Sorted_Tabdelim="/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/06.post_imputation/chrpos_rsid_biallelic/chrpos_rsid_biallelic_chr1-22.txt"
file=`basename $AllChr_Sorted_Tabdelim`
refshort="MikeDB_db155"

$plink2 \
    --bfile ${imputeFile}_UpdateChrPos_NoDuplicates \
    --update-name $AllChr_Sorted_Tabdelim \
    --threads 16 \
    --make-bed --out ${imputeFile}_RSid.$refshort

# remove remaining chr:pos SNPs
grep ":" ${imputeFile}_RSid.$refshort.bim  > ${imputeFile}_RSid.${refshort}_listremainingchrpos.txt
$plink2 \
    --bfile ${imputeFile}_RSid.$refshort \
    --exclude ${imputeFile}_RSid.${refshort}_listremainingchrpos.txt \
    --threads 16 \
    --make-bed --out ${imputeFile}_RSidOnly.$refshort

