#!/usr/bin/env python

#For use with python 2
from __future__ import print_function
import re
import sys
import os
import argparse


#This script will iterate through finished PAML directory, and extract the number of species and alignement lengths from the hog.phy files. The first argument should specify the hogfile, and the second argument should specify the desired name of the results file

def get_align_info(hog):
    #Extract relative path (needs to start from PAML working directory)
    toppath = '{:0>4}'.format(int(hog) % 100)
    # 0000/100/
    fullpath = toppath + "/" + hog + "/"
    
    align_file = '%s%s.phy'%(fullpath,hog)
    
    with open(align_file) as align:
        first_line = align.readline()
        
    hog_info = first_line.strip().split(" ")
    
    return(hog_info)
    


def main():
    hogfile = sys.argv[1]
    resfile_foroutput = sys.argv[2]
    
    #Extract hogs
    with open(hogfile) as hfile:
        hogs=[line.strip() for line in hfile]
    
    #Open res file
    resfile = open(resfile_foroutput,"w")
    resfile.write('hog\tnseq\tlength\n')
    
    res_dict = {}
    for hog in hogs:
        res_dict[int(hog)] = get_align_info(hog)
        
    
    for hog in sorted(res_dict):
        resfile.write('%d\t%s\t%s\n'%(hog,res_dict[hog][0],res_dict[hog][1]))

    resfile.close()



if __name__ == "__main__":
    main()