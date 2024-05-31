#!/bin/sh
# prune, QC for each ancestry

#$ -wd /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/03.ancestry_estimation/splitted_by_ancestry_groups
#$ -l h_rt=24:00:00,h_data=4G,highp
#$ -pe shared 16
#$ -j y
#$ -o joblog/112922_joblog.$JOB_ID.040-050.prune_qc_GRM_PCA_each_ancestry_AFR.log
#$ -M $USER@mail
#$ -m bea
#$ -t 1-5:1

# qsub
pops=("EUR" "AMR" "AFR" "EAS" "SAS")
pop=${pops[$SGE_TASK_ID-1]}

# directories
cohort="ABCDr4"
ancestrydir=/u/project/gandalm/emmamk/ASD_GWAS/$cohort/03.ancestry_estimation/splitted_by_ancestry_groups
grmdir=/u/project/gandalm/emmamk/ASD_GWAS/$cohort/04.grm
pcdir=/u/project/gandalm/emmamk/ASD_GWAS/$cohort/05.pc_each_ancestry

cohort_shared="ABCD_Release_4"
shared_dir=/u/project/gandalm/shared/GenomicDatasets-processed/$cohort_shared/genotype
shared_genotype_qced_dir=$shared_dir/genotype_QCed_each_ancestry
shared_grmdir=$shared_dir/GRM
shared_pcdir=$shared_dir/ancestry_PCs
mkdir -p $grmdir/log $pcdir/log $shared_genotype_qced_dir $shared_grmdir $shared_pcdir $shared_pcdir_20

# software
. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2
gcta=/u/project/gandalm/shared/apps/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static

# pruning for each ancestry
file=`ls -1v $ancestrydir/*rmsexmismatch*fam | grep $pop | sed 's/\.fam//g'`
filebase=`basename $file`

# create ${filebase}_indep.50.5.2.prune.in
plink \
    --bfile $ancestrydir/$filebase \
    --indep 50 5 2 \
    --keep-allele-order \
    --out $grmdir/${filebase}_indep.50.5.2

# prune file
$plink2 \
    --bfile $ancestrydir/$filebase \
    --extract $grmdir/${filebase}_indep.50.5.2.prune.in \
    --threads 16 \
    --memory 64000 \
    --make-bed --out $grmdir/${filebase}_indep.50.5.2_pruned

# QC
$plink2 \
    --bfile $grmdir/${filebase}_indep.50.5.2_pruned \
    --geno 0.1 --maf 0.01 --hwe 1e-6 \
    --threads 16 \
    --memory 64000 \
    --make-bed --out $grmdir/${filebase}_indep.50.5.2_pruned_geno0.1.maf0.01.hwe

$plink2 \
    --bfile $grmdir/${filebase}_indep.50.5.2_pruned_geno0.1.maf0.01.hwe \
    --mind 0.1 \
    --genotyping-rate \
    --threads 16 \
    --memory 64000 \
    --make-bed --out $grmdir/${filebase}_indep.50.5.2_pruned_geno0.1.maf0.01.hwe.mind0.1

# create GRM 
file="${filebase}_indep.50.5.2_pruned_geno0.1.maf0.01.hwe.mind0.1"
$gcta \
    --thread-num 16 \
    --bfile $grmdir/$file \
    --make-grm \
    --autosome \
    --out $grmdir/${file}_GRM

# create sparse GRM for family data
$gcta \
    --thread-num 16 \
    --grm $grmdir/${file}_GRM \
    --make-bK-sparse 0.05 \
    --out $grmdir/${file}_GRMsparse   

# bK GRM for genetic correlation
$gcta \
    --thread-num 16 \
    --grm $grmdir/${file}_GRM \
    --make-bK 0.05 \
    --out $grmdir/${file}_GRMbK

# create local ancestry PCs file
$plink2 \
    --bfile $grmdir/$file \
    --pca \
    --threads 16 \
    --memory 64000 \
    --out $pcdir/${file}.PCA10

# share files
mv $grmdir/${filebase}_indep.50.5.2_pruned_geno0.1.maf0.01.hwe.mind0.1.{bed,bim,fam} $shared_genotype_qced_dir
mv $grmdir/${file}_GRM* $shared_grmdir
mv $pcdir/${file}.PCA10* $shared_pcdir

