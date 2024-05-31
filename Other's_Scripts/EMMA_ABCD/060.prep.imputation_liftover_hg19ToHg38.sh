#!/bin/sh
# prep for imputation
# lifting over the array data from hg19 to hg38 using CrossMap.py

# directories
ref_chain=/u/project/gandalm/shared/refGenomes/hg19/hg19ToHg38.over.chain
ref_fasta=/u/project/gandalm/shared/refGenomes/hg38/v33/GRCh38.p13.genome.fa
indir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink
outdir=$indir/liftover_hg19ToHg38
name="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm.chr1.22.maf0.01.geno0.05.hwe1e6.mind0.1.refallele.reset"
mkdir -p $outdir

# software
. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624

# plink -> vcf
plink \
    --bfile $indir/$name \
    --keep-allele-order \
    --recode vcf-iid bgz \
    --out $outdir/$name.hg19ToHg38.plink2vcf

# hg19 -> hg38
source /u/home/e/emmamk/miniconda3/etc/profile.d/conda.sh
conda activate
CrossMap.py \
    vcf $ref_chain $outdir/$name.hg19ToHg38.plink2vcf.vcf.gz \
    $ref_fasta \
    $outdir/$name.hg19ToHg38.plink2vcf.CrossMap.vcf \
    --chromid l

# vcf -> plink
plink \
    --vcf $outdir/$name.hg19ToHg38.plink2vcf.CrossMap.vcf \
    --double-id \
    --keep-allele-order \
    --make-bed --out $outdir/$name.hg19ToHg38.plink2vcf.CrossMap.vcf2plink

# compress $name.hg19ToHg38.plink2vcf.CrossMap.vcf
bgzip=/u/project/gandalm/shared/apps/tabix-0.2.6/bgzip
$bgzip $outdir/$name.hg19ToHg38.plink2vcf.CrossMap.vcf && rm $outdir/$name.hg19ToHg38.plink2vcf.CrossMap.vcf

