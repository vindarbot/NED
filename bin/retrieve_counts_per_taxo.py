#!/usr/bin/env python3
import os
import re
import sys
import gzip
import json
import inspect
import argparse
from collections import defaultdict

def read_files(abundance_file, design_file):
	abondances = open(abundance_file,'r').readlines()
	if design_file is not None:
		design_exp = open(design_file,'r').readlines()
	else:
		design_exp = None
	return abondances, design_exp

def retrieve_proportions(abundances, design, ranks, output_file):
	rank_to_lvl = {
		"Domain" : 0, 
		"Phylum" : 1, 
		"Class" : 2, 
		"Order" : 3, 
		"Family" : 4, 
		"Genus" : 5, 
		"Species" : 6
	}

	output = open(output_file,'w')

	for rank in ranks:
		try:
			rank_to_lvl[rank]
		except:
			exit(rank+' is not a valid taxonomic rank. Valid ranks: Domain, Phylum, Class, Order, Family, Genus, Species')
		output.write('#'+rank+"\n")
		total = 0

		columns = []
		if design is not None:
			samples = []
			for li in design:
				li = li.strip().split('\t')
				samples.append(li[0])
			for i in range(len(abundances[0].strip().split('\t'))):
				if abundances[0].strip().split('\t')[i] in samples:
					columns.append(i)
		else:
			for i in range(8,len(abundances[0].strip().split('\t'))):
				columns.append(i)	

		taxo_to_count = defaultdict(int)
		for li in abundances[1:]:
			li = li.strip().split('\t')
			try:
				taxo = li[0].split(';')[rank_to_lvl[rank]]
			except:
				taxo = "no data"

			for i in columns:
				total+= int(li[i])
				taxo_to_count[taxo] += int(li[i])

		for a,b in taxo_to_count.items():
			output.write(a+'\t'+str(b)+'\t'+str(b/total*100)+"\n")
		output.write('\n')

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser(description="Retrieve proportion of given taxonomic ranks.")
	# Inputs
	group_input = parser.add_argument_group('Inputs')
	group_input.add_argument('-i', '--input', required=True, help="Input tsv abundance ouput file of FROGS (affiliation_OTU step). Taxonomic affiliations must be in the first column")
	group_input.add_argument('-d', '--design', help='If this parameter is specified, taxonomic proportions are given for each group of samples. This file must indicate in the first column the name of the sample, and in the second sample the group of the sample.')
	group_input.add_argument('-r', '--ranks', default="Domain,Phylum,Class,Order,Family,Genus,Species", help='Taxonomics ranks to analyse (default: all ranks analyse.) Please separate the ranks by commas (ex Domain,Phylum,Species ).')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('-o', '--output', default='taxonomics_proportion.tsv', help='Output file.')
	args = parser.parse_args()

	abundances, design = read_files(args.input, args.design)

	ranks = args.ranks.split(',')

	retrieve_proportions(abundances, design, ranks, args.output)


	
