#!/bin/sh
# remove samples on plate461 (n=82)
# check if sample matching issue IDs present (n=2) -- not present in genotype data
# ref: /u/project/gamdalm/shared/GenomicDatasets/ABCD_Release_4/release_notes/9. NDA 4.0 Genetics.pdf

# directories
date="101422"
tempdir=$SCRATCH/ABCDr4_QC_$date
genotypedir=/u/project/gandalm/shared/GenomicDatasets/ABCD_Release_4/genomics_sample03
qcdir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
qclogdir=qcdir/log
filebase="ABCD_release_3.0_QCed"
mkdir -p $qclogdir $tempdir

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# Get a list of samples on plate461
# create lists of unrecommended plate sample and samples with subject matching issue
grep $'\t461' $genotypedir/ABCD_release3.0_.batch_info.txt | awk 'BEGIN {OFS = "\t"} {print $1"_"$2, $1"_"$2}' > $genotypedir/plate461.sample.NDARids.txt
bad_plate_file=$genotypedir/plate461.sample.NDARids.txt

wc -l $bad_plate_file
# 82 /u/project/gandalm/shared/GenomicDatasets/ABCD_Release_4/genomics_sample03/plate461.sample.NDARids.txt

# sanity check: if the file have inconsistent IDs
grep -nE "INVF3FYXH1G|INVPWLFYWTX" $bad_plate_file
# empty: no problem

# remove samples on plate461
$plink2 \
    --bfile $tempdir/$filebase.NDARids.2idfixed \
    --remove $bad_plate_file \
    --make-bed --out $tempdir/$filebase.NDARids.2idfixed.plate461rm

# Sanity check for 2 sample IDs with subject matching issues
# ref: /u/project/gandalm/shared/GenomicDatasets/ABCD_Release_4/genomics_sample03/imputed/SUBJ_QC_BAD.txt
awk 'BEGIN {OFS = "\t"} {print $1"_"$2, $1"_"$2}' $genotypedir/imputed/SUBJ_QC_BAD.txt > $genotypedir/imputed/SUBJ_QC_BAD.NDARids.txt
bad_sample_file=$genotypedir/imputed/SUBJ_QC_BAD.NDARids.txt

# sanity check: if the file have inconsistent IDs
grep -nE "INVF3FYXH1G|INVPWLFYWTX" $bad_sample_file
# empty: can go ahead

# remove samples with sample matching issue
$plink2 \
    --bfile $tempdir/$filebase.NDARids.2idfixed.plate461rm \
    --remove $bad_sample_file \
    --make-bed --out $qcdir/$filebase.NDARids.2idfixed.plate461rm.2badsubjects.rm
