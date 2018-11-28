#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
import csv

#Takes the first column of an input csv file and returns as a list, skips first line
def grab_sp(spnamesfile):
    spnames = open(spnamesfile,"r")
    sp = []
    count = 1
    for line in spnames:
        if count == 1:
            count += 1
        else:
            line = line.strip()
            line = line.split(",")
            sp.append(line[0])
            count += 1
    spnames.close()
    return(sp)


def pull_sp_pval(res_list,species,type="pval"):

    if type == "pval":
        res = ([float(res[6]) for res in res_list if species in res[0]])
    elif type == "pvalholm":
        res = ([float(res[7]) for res in res_list if species in res[0]])
    elif type == "omega":
        res = ([float(res[3]) for res in res_list if species in res[0]])
    elif type == "weight":
        res = ([float(res[4]) for res in res_list if species in res[0]])
    elif type == "mean_dnds":
        res = ([float(res[1]) for res in res_list if species in res[0]])
    elif type == "rate_classes":
        res = ([float(res[2]) for res in res_list if species in res[0]])
    elif type == "branch_length":
        res = ([float(res[8]) for res in res_list if species in res[0]])
    elif type == "lrt":
        res = ([float(res[5]) for res in res_list if species in res[0]])
    else:
        print("parameter type must be pval, pvalholm, omega, weight, mean_dnds, rate_classes, branch_length, or lrt")

    if len(res) < 1:
        return('NA')

    else:
        if type == "pval" or type == "pvalholm":
            res = str(min(res))
        else:
            res = str(max(res))
        return(res)


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--hog", help="hog to analyze")
    parser.add_argument("--results_file", help="bsrel resultsfile to analyze")
    parser.add_argument("--spnames_file", help="file of species to analyze, expects code, or abbreviation to be in the first column")
    parser.add_argument("--parameter", help="pval (raw pvalue), pvalholm (FDR adjusted), omega (OmegaOver1), weight (WtOmegaOver1), mean_dnds (Mean_dNdS), rate_classes (RateClasses), branch_length (BranchLength), lrt (LRT). Parameter names reflect aBS-REL output")

    args = parser.parse_args()

    hog = args.hog
    results_file = args.results_file
    spnames_file = args.spnames_file
    parameter = args.parameter

    if parameter != "pval" and parameter != "pvalholm" and parameter != "omega" and parameter != "weight" and parameter != "mean_dnds" and parameter != "rate_classes" and parameter != "branch_length" and parameter != "lrt":
        print("parameter must be pval (raw pvalue), pvalholm (FDR adjusted), omega (OmegaOver1), weight (WtOmegaOver1), mean_dnds (Mean_dNdS), rate_classes (RateClasses), branch_length (BranchLength), or lrt (LRT)")
        quit()

    #Get list of species names
    sp = grab_sp(spnames_file)

    #Look for the presence of a results file for a given hog, if not present, print NAs

    all_na = ["NA"]*len(sp)

    if not os.path.isfile(results_file):
        no_data = "\t".join(all_na)
        print('%s\t%s'%(hog,no_data))
        quit()

    #now read results into list
    with open(results_file, 'r') as f:
        reader=csv.reader(f)
        res_list=list(reader)

    param_list = []

    for i in range(0,len(sp)):
        species = sp[i]
        param_list.append(pull_sp_pval(res_list,species,type=parameter))

    params = "\t".join(param_list)
    print('%s\t%s'%(hog,params))


if __name__ == "__main__":
    main()
