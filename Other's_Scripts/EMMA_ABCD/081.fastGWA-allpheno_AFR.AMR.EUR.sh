#!/bin/bash
# fastGWA for All phenotypes (CBCL, SRS) for AFR, AMR, EUR

#$ -wd /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/08.fastgwa
#$ -l h_rt=4:00:00,h_data=4G,highp
#$ -pe shared 4
#$ -j y
#$ -o joblog/121422_$JOB_ID.$TASK_ID_fastGWA_age.as.is_each.pop
#$ -M $USER@mail
#$ -m bea
#$ -t 1-3:1

date="121422"
pops=("EUR" "AMR" "AFR")
pop=${phenotypes[$SGE_TASK_ID-1]}
phenotypes="cbcl srs"

individuals="AllSubj"
group="ancestry_group"
timepoint="1yr_followup"

# directories
imputedPlinkDir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed/splitted_by_ancestry_groups
grmDir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/GRM
phenodataDir=/u/project/gandalm/shared/GenomicDatasets-processed/ABCD_Release_4/phenotype

# covariate file
covarFile="$phenodataDir/covar_AllSubj_batch_gender_noNAs_baseline_11063.txt"
covar_short="batch_gender"

# software
gcta=/u/project/gandalm/shared/apps/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static

# imputed genotype file
dataFile=`ls $imputedPlinkDir/*hwe1e6.fam | sed 's/.fam//g' | grep $pop`
# grm file
grmFile=`ls $grmDir/*_GRMsparse.grm.id | sed 's/.grm.id//g' | grep $pop`

for pheno in $phenotypes; do
    if [ $pheno = "cbcl" ]; then
        score="cbcl_scr_syn_social_t_rint"
    elif [ $pheno = "srs" ]; then
        score="ssrs_p_ss_sum_rint"
    fi

    # phenotype file
    phenotypeFile=`ls $phenodataDir/Phen_${individuals}_${pheno}_${timepoint}_*${group}*`
    pheno_short="${pheno}_${score}"

    # qcovar file
    qcovarFiles=`ls $phenodataDir/qcovar_${individuals}_${pheno}_${timepoint}_age.as.is_*${group}*`
    for qcovarFile in $qcovarFiles; do
        qcovar_short=`basename $qcovarFile | sed "s/qcovar_${individuals}_${pheno}_//g" | sed -r 's/_within_ancestry_group_noNAs_[0-9]+.txt//g'`
    
        # output
        outDir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/08.fastgwa/$pheno
        outFile=$pop.$pheno_short.$covar_short.$qcovar_short.$date
        mkdir -p $outDir/log

        $gcta --fastGWA-mlm \
            --bfile $dataFile \
            --grm-sparse $grmFile \
            --pheno $phenotypeFile \
            --covar $covarFile \
            --qcovar $qcovarFile \
            --thread-num 16 \
            --out $outDir/$outFile

        # compress output
        /usr/bin/gzip $outDir/$outFile.fastGWA

        echo " "
        echo " "
    done    
done

