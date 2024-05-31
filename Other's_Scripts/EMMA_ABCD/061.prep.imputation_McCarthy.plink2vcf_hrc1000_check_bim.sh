#!/bin/sh
# prep for imputation with McCarthy Tools

# directories
indir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink/liftover_hg19ToHg38
filebase="ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm.chr1.22.maf0.01.geno0.05.hwe1e6.mind0.1.refallele.reset.hg19ToHg38.plink2vcf.CrossMap.vcf2plink"
mkdir -p $indir/log

# software
. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624
bgzip=/u/project/gandalm/shared/apps/tabix-0.2.6/bgzip

# generate .frq file
plink --bfile $indir/$filebase --freq --out $indir/$filebase

# run HRC-1000G-check-bim.pl
reffile="/u/project/gandalm/shared/refGenomes/TOPMed/ALL.TOPMed_freeze8_hg38_dbSNP.tab.gz"
perl /u/project/gandalm/shared/apps/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl \
    -b $indir/$filebase.bim \
    -f $indir/$filebase.frq \
    -r $reffile \
    -h

# run Run-plink.sh (generate QCed and updated vcf files)
sed "s/--recode vcf/--recode vcf-iid/g" $indir/Run-plink.sh > $indir/Run-plink_vcfiid.sh
chmod +x $indir/Run-plink_vcfiid.sh && $indir/Run-plink_vcfiid.sh

# move files to temp_topmed_upload
cd $indir
mkdir -p McCarthyTools/temp_topmed_upload
mv *updated*vcf McCarthyTools/temp_topmed_upload

# add "chr" to vcf files
cd $indir/McCarthyTools/temp_topmed_upload
vcf_files=`ls -1v *vcf | sed 's/.vcf$//g'`
for vcf_file in $vcf_files; do
    echo $vcf_file
    awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $vcf_file.vcf | $bgzip > $vcf_file.withChr.vcf.gz && rm $vcf_file.vcf
done

