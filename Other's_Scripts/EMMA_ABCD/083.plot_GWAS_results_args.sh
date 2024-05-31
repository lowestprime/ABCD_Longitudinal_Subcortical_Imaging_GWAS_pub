#!/bin/bash
# qsub /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/00.scripts/08.fastGWA/qsub_qq.manhattan_plot_GMjulia_allpheno_noNAs_AFR.AMR.EUR_both.gender.together_biallelic_age.as.is_121522.sh

#$ -wd /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/08.fastgwa
#$ -l h_rt=4:00:00,h_data=4G,highp
#$ -pe shared 4
#$ -j y
#$ -o joblog/job.$JOB_ID.$TASK_ID.$JOB_NAME.log
#$ -M $USER@mail
#$ -m bea
#$ -t 1-2:1
# number of tasks: echo ${#conditions[*]}

# software $ scripts
julia=/u/project/gandalm/shared/apps/julia-1.8.1/bin/julia
scriptdir="/u/project/gandalm/emmamk/tools/GeneticsMakie_arrange"

# args
dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/08.fastgwa
conditions="$dir/conditions_plot.txt"
condition_args=`sed -n ${SGE_TASK_ID}P $conditions`
pheno=`echo $condition_args | cut -f1 -d"_"`
condition=`echo $condition_args | cut -f2 -d"_"`
keyword=`echo $condition_args | cut -f3 -d"_"`
keyword=`echo $keyword | sed -e "s/notImpaired/not_impaired/g; s/Impaired/impaired/g"`

cohort="ABCDr4"
qq_numrow=1
qq_numcol=3
projectdir="/u/project/gandalm/emmamk/ASD_GWAS"
software="fastgwa"
gwasdir="$projectdir/$cohort/08.$software/$pheno"
metaldir="NA"
plotdir=$gwasdir
samplesizedir="NA"
reverse_files="NA"

# plot script
script="manhattan_qq_locuszoom_plot.jl"
$julia $scriptdir/$script \
$cohort $pheno $condition $keyword $qq_numrow $qq_numcol $projectdir $gwasdir $metaldir $plotdir $samplesizedir $reverse_files

