#!/bin/sh
# reset REF allele

# directories
indir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
reffile_dir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype
filebase="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm.chr1.22.maf0.01.geno0.05.hwe1e6.mind0.1"

# software
plink2=/u/project/gandalm/shared/apps/plink2/v2.00a3-20220814/plink2

# using the reference file created with 022.create_chrpos_rsID_replacement_files_hg19_b155_Mikedb.ipynb
ref_file="221117_ABCDr4_originalbim_variants_rsIDs_original.flipped.unmapped.merged_mikedb.txt"
$plink2 \
	--bfile $indir/$filebase \
	--ref-allele $reffile_dir/$ref_file 4 9 \
	--threads 16 \
	--make-bed --out $indir/$filebase.refallele.reset

