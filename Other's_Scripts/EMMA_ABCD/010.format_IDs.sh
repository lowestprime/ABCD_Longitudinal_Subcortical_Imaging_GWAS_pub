#!/bin/sh
# format IDs in ABCD data

# directories
date="101422"
indir=/u/project/gandalm/shared/GenomicDatasets/ABCD_Release_4/genomics_sample03/genotype_QCed
qclogdir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink/log
tempdir=$SCRATCH/ABCDr4_QC_$date
filebase="ABCD_release_3.0_QCed"
mkdir -p $qclogdir $tempdir

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# create a list for replacing family IDs with NDAR IDs
awk 'BEGIN {OFS = "\t"} {print $1, $2, $1"_"$2, $1"_"$2}' $indir/$filebase.fam > $indir/$filebase.NDARids.txt

# replace FIDs and IIDs in phenotype data format
cd $indir
$plink2 \
	--bfile $filebase \
	--update-ids $filebase.NDARids.txt \
	--make-bed --out $tempdir/$filebase.NDARids
mv $tempdir/$filebase.NDARids.log $qclogdir/01.$filebase.NDARids.log

# fix 2 sample IDs that have inconsistent format
grep -nv "_NDAR_" $tempdir/$filebase.NDARids.fam
# 2881:AB0003488_`NDAR_INVF3FYXH1G        AB0003488_`NDAR_INVF3FYXH1G     0       0       0       -9
# 6184:AB0007279_NDARINVPWLFYWTX          AB0007279_NDARINVPWLFYWTX       0       0       0       -9

grep -v "_NDAR_" $tempdir/$filebase.NDARids.fam | cut -f1,2 > id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_before.txt
grep "\`" $tempdir/$filebase.NDARids.fam | sed 's/`//g' | cut -f1,2 > id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_after.txt
echo -e "AB0007279_NDAR_INVPWLFYWTX\tAB0007279_NDAR_INVPWLFYWTX" >> id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_after.txt
paste id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_before.txt id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_after.txt > id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_replacement.txt && \
rm id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_before.txt id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_after.txt

# replace the 2 IDs using the replacement file
cd $indir
$plink2 \
	--bfile $tempdir/$filebase.NDARids \
	--update-ids id_inconsistent_INVF3FYXH1G_INVPWLFYWTX_replacement.txt \
	--make-bed --out $tempdir/$filebase.NDARids.2idfixed

# sanity check for the 2 IDs
grep -nE "INVF3FYXH1G|INVPWLFYWTX" $tempdir/$filebase.NDARids.2idfixed.fam
# 2881:AB0003488_NDAR_INVF3FYXH1G AB0003488_NDAR_INVF3FYXH1G      0       0       0       -9
# 6184:AB0007279_NDAR_INVPWLFYWTX AB0007279_NDAR_INVPWLFYWTX      0       0       0       -9
