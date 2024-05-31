#!/bin/bash
# create condition list for plotting

dir=/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/08.fastgwa
pheno_list="cbcl srs"
conditions_age="withAge"
conditions_sex="withSex"
conditions_pcs="PCs1-10"
conditions_batch="withBatch"
conditions_piq="noPIQ"
for pheno in $pheno_list; do
	for pcs in $conditions_pcs; do
		for batch in $conditions_batch; do
			for piq in $conditions_piq; do
				# pheno - condition - keyword
				echo "${pheno}_$conditions_age.$conditions_sex.$conditions_pcs.$conditions_batch.${conditions_piq}_${pheno}*${pcs}" >> $dir/conditions_plot.txt
			done
		done
	done
done

