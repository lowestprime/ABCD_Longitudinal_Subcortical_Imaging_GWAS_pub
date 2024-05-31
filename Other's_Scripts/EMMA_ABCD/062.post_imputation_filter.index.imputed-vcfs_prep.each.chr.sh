#!/bin/sh
# 1. filter imputed data r2 > 0.8
# 2. index filtered file with tabix

#$ -wd /u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed/vcf_imputed_raw
#$ -o joblog/120322_job$JOB_ID.$TASK_ID.filter.merge.index.imputed-vcfs_re.log
#$ -j y
#$ -l h_rt=80:00:00,h_data=4G,highp
#$ -pe shared 16
#$ -M $USER@mail
#$ -m bea
#$ -t 1-22:1

# directories
indir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed/vcf_imputed_raw
outdir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed/vcf_imputed_filtered
mkdir -p $outdir

chr=$SGE_TASK_ID

# software
. /u/local/Modules/default/init/modules.sh
module load bcftools/1.11
tabixdir=/u/project/gandalm/shared/apps/tabix-0.2.6
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# filter each chr files with r2 > 0.8
bcftools view --threads 16 -i 'R2>.8' -Oz -o $outdir/chr${chr}.imp.flt.r208.vcf.gz $indir/chr${chr}.dose.vcf.gz

# create index file
$tabixdir/tabix -p vcf $outdir/chr${chr}.imp.flt.r208.vcf.gz

