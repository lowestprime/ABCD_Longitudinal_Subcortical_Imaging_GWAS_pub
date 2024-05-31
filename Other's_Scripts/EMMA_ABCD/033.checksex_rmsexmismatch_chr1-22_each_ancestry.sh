#!/bin/sh
# check sex mismatch and remove mismatched samples

#$ -wd /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/03.ancestry_estimation/splitted_by_ancestry_groups
#$ -l h_rt=24:00:00,h_data=4G,highp
#$ -pe shared 16
#$ -j y
#$ -o joblog/112922_joblog.$JOB_ID.$TASK_ID_036_rm.sexmismatch_chr1-22_each_ancestry.log
#$ -M $USER@mail
#$ -m bea
#$ -t 1-5:1

pops=("EUR" "AFR" "AMR" "EAS" "SAS")
pop=${pops[$SGE_TASK_ID-1]}

# directories
indir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/03.ancestry_estimation/splitted_by_ancestry_groups
cohort="ABCDr4"
date="112922"
tempdir=$SCRATCH/${cohort}_checksex_$date
mkdir -p $tempdir $indir/log

# files
file=`ls -1v $indir/*RSidOnly*$pop*fam`
filebase=`basename $file | sed 's/\.fam//g'`

# software
. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# --split-x (--split-par) to use check sex later
$plink2 \
	--bfile $indir/$filebase \
	--threads 16 \
	--memory 64000 \
	--split-par b37 \
	--genotyping-rate \
	--allow-extra-chr \
	--make-bed --out $tempdir/$filebase.splitpar

# --indep-pairphase 20000 2000 0.5 to use check sex later
plink \
	--bfile $tempdir/$filebase.splitpar \
	--memory 64000 \
	--indep-pairphase 20000 2000 0.5 \
	--allow-extra-chr \
	--keep-allele-order \
	--make-bed --out $tempdir/$filebase.splitpar.indep.pairphase.20000.2000.0.5

# check sex data in the plink file
plink \
	--bfile $tempdir/$filebase.splitpar \
	--memory 64000 \
	--extract $tempdir/$filebase.splitpar.indep.pairphase.20000.2000.0.5.prune.in \
	--check-sex \
	--keep-allele-order \
	--allow-extra-chr \
	--out $tempdir/$filebase.splitpar.indep.pairphase.20000.2000.0.5.pruned.cksex

# Remove IDs mismatched on sex
# samples reported male but F < 0.2 OR reported female but F > 0.8
mv $tempdir/$filebase.splitpar.indep.pairphase.20000.2000.0.5.pruned.cksex.sexcheck $indir
num=`awk '($3 == 1 && $4 == 2 || $3 == 2 && $4 == 1) {print}' $indir/$filebase.splitpar.indep.pairphase.20000.2000.0.5.pruned.cksex.sexcheck | wc -l | cut -f1 -d" "`
awk '($3 == 1 && $4 == 2 || $3 == 2 && $4 == 1) {print}' $indir/$filebase.splitpar.indep.pairphase.20000.2000.0.5.pruned.cksex.sexcheck > $indir/sex.drop.$pop.$num.txt

# remove sex mismatched samples $ restrict to autosomes
$plink2 \
	--bfile $indir/$filebase \
	--threads 16 \
	--memory 64000 \
	--remove $indir/sex.drop.$pop.$num.txt \
	--chr 1-22 \
	--genotyping-rate \
	--allow-extra-chr \
	--make-bed --out $indir/$filebase.rmsexmismatch.$num.chr1.22

