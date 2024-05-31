#!/bin/sh
# extract chr1-22 for ancestry estimation

# directories
date="101422"
tempdir=$SCRATCH/ABCDr4_QC_$date
outdir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
filebase="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm"

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

$plink2 \
	--bfile $tempdir/$filebase \
	--chr 1-22 \
	--genotyping-rate \
	--make-bed --out $tempdir/$filebase.chr1.22
