#!/usr/bin/env python


from __future__ import print_function
import os
import sys
import io
import re
import StringIO


#Code to work with Python2
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def parse_busted_results(file):
     #takes a busted output file, splits on TREE #, returns list of output files to parse
    
    busted_results=["NA","NA","NA"]    
   
    with open(file, 'r') as pf:
        for line in pf:
        	line=line.strip()
        	if "Unrestricted class omega" in line:
        		newline = line.split("= ")
        		busted_results[1]=str(newline[2][0:15])
        		busted_results[2]=str(newline[3][0:17])
        	elif "Likelihood ratio test for episodic positive selection" in line:
        		newline = line.split("p = ")
        		busted_results[0]=str(newline[1][0:17])
        	elif "No evidence for positive selection under the unconstrained model" in line:
        		busted_results[0]=str(1)
        	else:
        		pass
    
    
    return(busted_results)


def parse_hogs(hoglist,model,input_dirs,verbose=True):
    #take list of hogs, return parsed final results dictionary
    final_results = {}
    for hog in hoglist:
        if verbose:
            print("Working on", hog)
            
        toppath = '{:0>4}'.format(int(hog) % 100)
        # 0000/100/
        fullpath = toppath + "/" + hog + "/"
        for pamldir in input_dirs:
            results_file_tree1 = pamldir + "/" + fullpath + "tree1/" + "output"
            results_file_tree2 = pamldir + "/" + fullpath + "tree2/" + "output"
                    
            try:
            	parsed_results_tree1 = parse_busted_results(results_file_tree1)
            except FileNotFoundError:
                print("Couldn't parse file for", hog, "at", pamldir + "/" + fullpath)
                continue
            try:
            	parsed_results_tree2 = parse_busted_results(results_file_tree2)
            except FileNotFoundError:
                print("Couldn't parse file for", hog, "at", pamldir + "/" + fullpath)
                continue            
                          
            final_results[hog] = {'tree1':parsed_results_tree1, 'tree2':parsed_results_tree2}
    
    return(final_results)
    
def print_results(results, handle, model):
    #results is a complicated dict, but the basic format is: hog, tree1->[tree 1 results], tree2->[tree 2 results]
    #to test let's start by printing out: hog id, tree num, pval, omega, weight
    for hog in results.keys():
    	tree1_res = "\t".join(results[hog]['tree1'])
    	tree2_res = "\t".join(results[hog]['tree2'])
    	print(hog,"tree1",tree1_res,sep="\t",end="\n",file=handle)
    	print(hog,"tree2",tree2_res,sep="\t",end="\n",file=handle)
                
                
def main():
    print("Starting to parse.")
    hogfile_toparse = sys.argv[1]
    resfile_foroutput = sys.argv[2]
    with open(hogfile_toparse) as hfile:
        hogs=[line.strip() for line in hfile]
    
    print("Done getting files.")
    with open(resfile_foroutput, 'w') as ofile:
        print("hog", "treenum", "pval", "omega", "weight", sep="\t", end="\n", file=ofile)
        results = parse_hogs(hogs,"busted",["."],True)
        print_results(results, ofile, "busted")

if __name__ == "__main__":
    main()