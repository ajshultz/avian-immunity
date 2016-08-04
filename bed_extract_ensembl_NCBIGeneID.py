#! /usr/bin/env python

import sys
import os


#This script will read in a bed file with both ensemble IDs and NCBI gene IDs merged, and output a translation table of these IDs.

bedfile = sys.argv[1]
outfile = sys.argv[2]

bed = open(bedfile,"r")
output = open(outfile,"w")

output.write("NCBIGeneID\tENSEMBLID\n")

for line in bed:
	line = line.strip()
	if "GeneID" in line:
		if "ENSGALG" in line:
			line2 = line.split("|")
			#Extract NCBI Gene IDs
			geneidpartpos = [i for i, j in enumerate(line2) if "GeneID" in j]
			geneidpart = line2[geneidpartpos[0]].split(",")
			geneidpos = [i for i, j in enumerate(geneidpart) if "GeneID" in j]
			geneidpart = geneidpart[geneidpos[0]]
			geneid = geneidpart.split("GeneID:")[-1]
			
			#Extract Ensembl IDs
			ensidpartpos = [i for i, j in enumerate(line2) if "ENSGALG" in j]
			ensidpart = line2[ensidpartpos[0]].split(";")
			ensidpos = [i for i, j in enumerate(ensidpart) if "ENSGALG" in j]
			ensidpart = ensidpart[ensidpos[0]]
			ensid = ensidpart.split("ID=gene:")[-1]
			
			
			output.write(geneid+"\t"+ensid+"\n")



bed.close()
output.close()