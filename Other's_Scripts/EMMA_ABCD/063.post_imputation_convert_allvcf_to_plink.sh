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
TOPMed_dir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed
indir=$TOPMed_dir/vcf_imputed_filtered
vcf2plink_dir=$TOPMed_dir/vcf_imputed_filtered_plink
mkdir -p $vcf2plink_dir

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# convert vcf to plink
chr=$SGE_TASK_ID
$plink2 \
    --vcf $indir/chr${chr}.imp.flt.r208.vcf.gz \
    --double-id \
    --threads 16 \
    --max-alleles 2 \
    --make-bed --out $vcf2plink_dir/${imputeFile}
