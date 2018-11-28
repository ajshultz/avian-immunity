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

    #Tree 1
    RESULTS_TREE1_PVAL=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter pval)
    RESULTS_TREE1_PVALHOLM=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter pvalholm)
    RESULTS_TREE1_OMEGA=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter omega)
    RESULTS_TREE1_WEIGHT=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter weight)
    RESULTS_TREE1_MEAN_DNDS=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter mean_dnds)
    RESULTS_TREE1_RATE_CLASSES=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter rate_classes)
    RESULTS_TREE1_BRANCH_LENGTH=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter branch_length)
    RESULTS_TREE1_LRT=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT1 --spnames_file $SPNAMES_FILE --parameter lrt)

    #Tree 2
    RESULTS_TREE2_PVAL=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE  --parameter pval)
    RESULTS_TREE2_PVALHOLM=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --parameter pvalholm)
    RESULTS_TREE2_OMEGA=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --parameter omega)
    RESULTS_TREE2_WEIGHT=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --parameter weight)
    RESULTS_TREE2_MEAN_DNDS=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --parameter mean_dnds)
    RESULTS_TREE2_RATE_CLASSES=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --parameter rate_classes)
    RESULTS_TREE2_BRANCH_LENGTH=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --parameter branch_length)
    RESULTS_TREE2_LRT=$(python parse_bsrel_hog_branch_pvals.py --hog $HOG --results_file $OUTPUT2 --spnames_file $SPNAMES_FILE --parameter lrt)

    echo -e "pval\ttree1\t$RESULTS_TREE1_PVAL" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "pvalholm\ttree1\t$RESULTS_TREE1_PVALHOLM" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "omega\ttree1\t$RESULTS_TREE1_OMEGA" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "weight\ttree1\t$RESULTS_TREE1_WEIGHT" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "mean_dnds\ttree1\t$RESULTS_TREE1_MEAN_DNDS" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "rate_classes\ttree1\t$RESULTS_TREE1_RATE_CLASSES" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "branch_length\ttree1\t$RESULTS_TREE1_BRANCH_LENGTH" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "lrt\ttree1\t$RESULTS_TREE1_LRT" >> bsrel_res_parsed_pvals_${DATE}.txt

    echo -e "pval\ttree2\t$RESULTS_TREE2_PVAL" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "pvalholm\ttree2\t$RESULTS_TREE2_PVALHOLM" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "omega\ttree2\t$RESULTS_TREE2_OMEGA" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "weight\ttree2\t$RESULTS_TREE2_WEIGHT" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "mean_dnds\ttree2\t$RESULTS_TREE2_MEAN_DNDS" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "rate_classes\ttree2\t$RESULTS_TREE2_RATE_CLASSES" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "branch_length\ttree2\t$RESULTS_TREE2_BRANCH_LENGTH" >> bsrel_res_parsed_pvals_${DATE}.txt
    echo -e "lrt\ttree2\t$RESULTS_TREE2_LRT" >> bsrel_res_parsed_pvals_${DATE}.txt
done
