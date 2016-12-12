#!/usr/bin/env python


from __future__ import print_function
import os
import sys
import io
import re
import StringIO


def parse_busted_results (file):
     #takes a paml output file, splits on TREE #, returns list of output files to parse
     
    with open(file, 'r') as pf:
    	busted_results=["NA","NA","NA"]
        for line in pf:
        	if re.search("Unrestricted class omega",line):
        		newline = line.split("= ")
        		busted_results[1]=float(newline[1][0:15])
        		busted_results[2]=float(newline[2][0:17])
        	elif re. search("Likelihood ratio test for episodic positive selection",line):
        		newline = line.split("p = ")
        		busted_results[0]=float(newline[1][0:17]
        	else
        		pass
    
    
    return(busted_results)


def parse_hogs(hoglist,model,input_dirs,verbose=True,multisite=False):
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
            
            
            if hog in final_results:
                #append
                cur_len = len(final_results[hog]['trees'])
                if cur_len != len(final_results[hog]['results']):
                    print("Warning, something went wrong!!")
                    
                #update keys (tree numbers)
                new_trees = {int(x)+cur_len:parsed_trees[x] for x in parsed_trees.keys()}
                new_results = {int(x)+cur_len:parsed_results[x] for x in parsed_results.keys()}
                final_results[hog]['trees'].update(new_trees)
                final_results[hog]['results'].update(new_results)
                    
            else:       
                final_results[hog] = {'trees':parsed_trees, 'results':parsed_results}
    
    return(final_results)
    
def print_results (results, handle, model):
    #results is a complicated dict, but the basic format is: hog, tree->tree_dict, results->results_dict
    #tree_dict and results_dict are matched by key
    #to test let's start by printing out: hog id, tree num, foreground, species_tree vs gene_tree
    #original tree, then paml parameters
    for hog in results.keys():
        tree_len = len(results[hog]['trees'])
        res_len = len(results[hog]['results'])
#       print(tree_len, res_len)
#       print(results[hog]['trees'].keys())
#       print(results[hog]['results'].keys())
        for i in range(1,res_len+1):
            trees = results[hog]['trees'].get(i, {})
            res = results[hog]['results'].get(i, {})
            print(hog, model, i, trees.get('foreground'), trees.get('is_species_tree'), trees.get('original'), sep="\t", end="\t", file=handle)
            #parse results
            res_lnl = res.get('NSsites', {}).get(2, {}).get('lnL')
            res_siteclass = res.get('NSsites', {}).get(2, {}).get('parameters', {}).get('site classes')
            res_treelen = res.get('NSsites', {}).get(2, {}).get('tree length')
            print(res_lnl, res_treelen, sep="\t", end="", file=handle)
            #need to iterate over site class numbers
            try:
                for sc in sorted(res_siteclass):
                    print("\t" + res_siteclass[sc]['proportion'], res_siteclass[sc]['branch types']['foreground'], res_siteclass[sc]['branch types']['background'], sep="\t", end = "", file=handle)
                print(file=handle)
            except:
                print("\tNone", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', sep="\t", end="\n", file=handle)
                
                
def main():
    print("Starting to parse.")
    hogfile_toparse = sys.argv[1]
    resfile_foroutput = sys.argv[2]
    with open(hogfile_toparse) as hfile:
        hogs=[line.strip() for line in hfile]
    
    print("Done getting files.")
    with open(resfile_foroutput, 'w') as ofile:
        print("hog", "treenum", "pval", "omega", "weight", "model_num", "lnl", "treelen", "kappa", "omega", sep="\t", end="\n", file=ofile)
        results = parse_hogs(hogs,"sites",["."],True,multisite=False)
        print_results(results, ofile, "sites")

if __name__ == "__main__":
    main()