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
    
    try:
        with open(align_file) as align:
            first_line = align.readline()
        hog_info = first_line.strip().split(" ")
    except:
        print('no alignment file for hog %s'%hog)
        hog_info = ["",""]
    
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
    
    #Run each hog through function to grab alignment info, write to file
    for hog in hogs:
        hog_res = get_align_info(hog)
        resfile.write('%s\t%s\t%s\n'%(hog,hog_res[0],hog_res[1]))
        
        print('hog %s done'%hog)

    resfile.close()



if __name__ == "__main__":
    main()