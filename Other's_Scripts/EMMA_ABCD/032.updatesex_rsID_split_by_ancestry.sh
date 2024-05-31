#!/bin/sh
# update sex column in .bim
# chr:pos -> rsID

# directories
date="101422"
tempdir=$SCRATCH/ABCDr4_QC_$date
indir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
indir_log=$indir/log_prep_for_split_ancestry
ancestry_dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/03.ancestry_estimation
outdir=$ancestry_dir/splitted_by_ancestry_groups
mkdir -p $indir_log $outdir/log

# files
filebase="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm"
sexinfofile="$indir/sex_info.txt"
chrpos_rsid_file="/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/ABCDr4_SNPid_replacement_mikedb155_5661.txt"

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# assign self-reported sex info
$plink2 \
    --bfile $tempdir/$filebase \
    --update-sex $sexinfofile \
    --allow-extra-chr \
    --make-bed --out $indir/${filebase}_withsexinfo

# chrpos -> rsID
$plink2 --bfile $indir/${filebase}_withsexinfo \
        --update-name $chrpos_rsid_file \
        --make-bed --out $tempdir/${filebase}_withsexinfo_RSid

# extract the 2nd column, find dups, and remove
cut -f2 $tempdir/${filebase}_withsexinfo_RSid.bim | sort | uniq -d > $indir/${filebase}_withsexinfo_RSid_duplicatedsnps.snplist

# exclude dup SNPs
$plink2 \
    --bfile $tempdir/${filebase}_withsexinfo_RSid \
    --exclude $indir/${filebase}_withsexinfo_RSid_duplicatedsnps.snplist \
    --allow-extra-chr \
    --make-bed --out $tempdir/${filebase}_withsexinfo_RSid_NoDuplicates

# find remaining chr:pos SNPs and exclude
grep ":" $tempdir/${filebase}_withsexinfo_RSid_NoDuplicates.bim > $indir/${filebase}_withsexinfo_RSid_NoDuplicates_listremainingchrpos.txt
$plink2 \
    --bfile $tempdir/${filebase}_withsexinfo_RSid_NoDuplicates \
    --exclude $indir/${filebase}_withsexinfo_RSid_NoDuplicates_listremainingchrpos.txt \
    --allow-extra-chr \
    --make-bed --out $indir/${filebase}_withsexinfo_RSid_NoDuplicates_RSidOnly

# split into ancestry groups
pop_files=`ls $ancestry_dir/*knn*[0-9].txt`
for pop_file in $pop_files; do
    pop_suffix=`basename $pop_file | grep -oE "[A-Z]{3}.[0-9]+"`
    
    $plink2 \
        --bfile $indir/${filebase}_withsexinfo_RSid_NoDuplicates_RSidOnly \
        --keep $pop_file \
        --allow-extra-chr \
        --make-bed --out $outdir/${filebase}_withsexinfo_RSid_NoDuplicates_RSidOnly_$pop_suffix
done

