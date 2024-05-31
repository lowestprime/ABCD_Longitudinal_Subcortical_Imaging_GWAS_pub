#!/bin/sh
# ancestry estimation
# using plink2 as much as possible if the functions available

# directories
cohort="ABCDr4"
date="112822"
tempdir=$SCRATCH/${cohort}_ancestry_estimation_${date}
lab_refdir=/u/project/gandalm/shared/refGenomes/g1000/Phase3_ALL
refdir=/u/project/gandalm/emmamk/reference/1kg/GRCh37
qced_plink_dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
ancestry_results_dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/03.ancestry_estimation
study_name="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm.chr1.22.maf0.01.geno0.05.hwe1e6.mind0.1.refallele.reset_RSidOnly"
refname='ALL.autosomes.phase3'
mkdir -p $tempdir $ancestry_results_dir/log

# software
module load plink/1.90b624
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

#### 01. Match study genotypes and reference data ####
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' $lab_refdir/$refname.bim > $refdir/$refname.ac_gt_snps
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' $qced_plink_dir/$study_name.bim > $tempdir/$study_name.ac_gt_snps

$plink2 \
    --bfile $lab_refdir/$refname \
    --exclude $refdir/$refname.ac_gt_snps \
    --make-bed --out $refdir/$refname.no_ac_gt_snps

$plink2 \
    --bfile $qced_plink_dir/$study_name \
    --exclude $tempdir/$study_name.ac_gt_snps \
    --make-bed --out $tempdir/$study_name.no_ac_gt_snps

#### 02. Prune study data ####
$plink2 \
    --bfile $tempdir/$study_name.no_ac_gt_snps \
    --indep-pairwise 50 5 0.2 \
    --out $tempdir/$study_name.no_ac_gt_snps

$plink2 \
    --bfile $tempdir/$study_name.no_ac_gt_snps \
    --extract $tempdir/$study_name.no_ac_gt_snps.prune.in \
    --make-bed --out $tempdir/$study_name.pruned

#### 03. Filter reference data for the same SNP set as in study ####
$plink2 \
    --bfile $refdir/$refname.no_ac_gt_snps \
    --extract $tempdir/$study_name.no_ac_gt_snps.prune.in \
    --make-bed --out $tempdir/$refname.pruned

#### 04. Check and correct chromosome mismatch ####
awk 'BEGIN {OFS="\t"} NR==FNR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' \
$tempdir/$study_name.pruned.bim $tempdir/$refname.pruned.bim | sed -n '/^[XY]/!p' > $tempdir/$refname.toUpdateChr
head $tempdir/$refname.toUpdateChr  # empty

## skipped because .toUpdateChr was empty
# $plink2 \
# --bfile $tempdir/$refname.pruned \
# --update-chr $tempdir/$refname.toUpdateChr 1 2 \
# --make-bed --out $tempdir/$refname.updateChr
# mv $tempdir/$refname.updateChr.log $tempdir/04.$refname.updateChr.log
# mv $tempdir/04.$refname.updateChr.log $tempdir/log/

#### 05. Upate positions and flip alleles ####
#### 05-1. Position mismatch ####
# find variants with mis-matching chromosomal positions
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} ($2 in a && a[$2] != $4) {print a[$2],$2}' \
$tempdir/$study_name.pruned.bim $tempdir/$refname.pruned.bim > $tempdir/$refname.toUpdatePos
head -5 *toUpdatePos  # empty

#### 05-2. Possible allele flips ####
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
$tempdir/$study_name.pruned.bim $tempdir/$refname.pruned.bim > $tempdir/$refname.toFlip
head -5 *toFlip
# rs60048151

#### 05-3. Upate positions and flip alleles ####
# update the mismatching positions and possible allele-flips identified above
plink \
    --bfile $tempdir/$refname.pruned \
    --flip $tempdir/$refname.toFlip \
    --keep-allele-order \
    --make-bed --out $tempdir/$refname.flipped

#### 06. Remove mismatches ####
# identify and remove any alleles that do not match after allele flipping from the reference dataset
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
$tempdir/$study_name.pruned.bim $tempdir/$refname.flipped.bim > $tempdir/$refname.mismatch
head -3 *mismatch
# rs60048151

$plink2 \
    --bfile $tempdir/$refname.flipped \
    --exclude $tempdir/$refname.mismatch \
    --make-bed --out $tempdir/$refname.clean

#### 07. Merge study genotypes and reference data ####
plink \
    --bfile $tempdir/$study_name.pruned \
    --bmerge $tempdir/$refname.clean.bed $tempdir/$refname.clean.bim $tempdir/$refname.clean.fam \
    --keep-allele-order \
    --make-bed --out $tempdir/$study_name.merge.$refname

#### 08. PCA on the merged data of Ancestry Estimation pipeline) ####
$plink2 \
    --bfile $ancestry_results_dir/$study_name.merge.$refname \
    --pca \
    --memory 64000 \
    --threads 16 \
    --out $ancestry_results_dir/$study_name.$refname

#### 09. Create population categorical data ####
# ref: https://github.com/gandallab/1kg-ancestry
awk 'BEGIN {print "FID","IID","POP"} NR==FNR {pop[$2]=$7} !($2 in pop) {pop[$2]="SET"} NR>FNR {print $1,$2,pop[$2]}' \
/u/project/gandalm/shared/refGenomes/1000genomes/chrs/20130606_g1k.ped $ancestry_results_dir/$study_name.merge.$refname.fam \
> $ancestry_results_dir/$study_name.$refname.pop

