#!/bin/sh
# QC plink

# directories
date="101422"
tempdir=$SCRATCH/ABCDr4_QC_$date
outdir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
filebase="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm.chr1.22"

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# QC
$plink2 \
	--bfile $tempdir/$filebase \
	--maf 0.01 --geno 0.05 --hwe 1e-6 \
	--genotyping-rate \
	--make-bed --out $tempdir/$filebase.maf0.01.geno0.05.hwe1e6

$plink2 \
	--bfile $tempdir/$filebase.maf0.01.geno0.05.hwe1e6 \
	--mind 0.1 \
	--genotyping-rate \
	--make-bed --out $outdir/$filebase.maf0.01.geno0.05.hwe1e6.mind0.1

