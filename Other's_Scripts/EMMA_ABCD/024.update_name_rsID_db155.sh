#!/bin/sh
# preparing for the subsquent ancestry estimation pipeline
# remove chr:pos and filter out low quality SNPs and samples from study data

# directories
indir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
study_name="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm.chr1.22.maf0.01.geno0.05.hwe1e6.mind0.1.refallele.reset"
reffile_dir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# update chr:pos --> rsIDs
# using reference file created with 022.create_chrpos_rsID_replacement_files_hg19_b155_Mikedb.ipynb
$plink2 \
	--bfile $indir/$study_name \
	--update-name $reffile_dir/ABCDr4_SNPid_replacement_mikedb155_5661.txt \
	--make-bed --out $indir/${study_name}_RSid

# check duplicated SNPs
cut -f2 $indir/${study_name}_RSid.bim | sort | uniq -d > $indir/${study_name}_RSid.duplicatedsnps.snplist.txt
# empty; no dup SNPs
# when not empty - remove dup SNPs

# find remaining chr:pos SNPs and exclude
cat $indir/${study_name}_RSid.bim | grep ":" > $indir/${study_name}_RSid_listremainingchrpos.txt
$plink2 \
	--bfile $indir/${study_name}_RSid \
	--exclude $indir/${study_name}_RSid_listremainingchrpos.txt \
	--make-bed --out $indir/${study_name}_RSidOnly

