#!/bin/bash

##HYPHY SCRIPT##
##all data in /n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA

WORKDIR=$(pwd -P)
SPNAMES_FILE=$1
DATE=$(date +%F)

for HOG in $(cat all_hogs)
do
    SUBDIRNUM=`expr ${HOG##HOG} % 100` 
    printf -v SUBDIR "%04d" $SUBDIRNUM
    RUNDIR=$SUBDIR/$HOG
    OUTPUT1=/n/holylfs/LABS/edwards_lab/tsackton/ratite_work/06_protein_coding/HyPhy/BS-REL/$RUNDIR/tree1/$HOG.aBSREL.OUT
    OUTPUT2=/n/holylfs/LABS/edwards_lab/tsackton/ratite_work/06_protein_coding/HyPhy/BS-REL/$RUNDIR/tree2/$HOG.aBSREL.OUT

    RESULTS_TREE1_PVAL=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --pval_type pval)
    RESULTS_TREE1_PVALHOLM=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --pval_type pvalholm)
    RESULTS_TREE2_PVAL=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE  --pval_type pval)
    RESULTS_TREE2_PVALHOLM=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --pval_type pvalholm)

    echo -e "pval\ttree1\t$RESULTS_TREE1_PVAL" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "pvalholm\ttree1\t$RESULTS_TREE1_PVALHOLM" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "pval\ttree2\t$RESULTS_TREE2_PVAL" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "pvalholm\ttree2\t$RESULTS_TREE2_PVALHOLM" >> bsrel_res_parsed_pvals_${DATE}.txt
done
