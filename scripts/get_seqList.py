#!/usr/bin/env python

"""
    @Arthor: Jitong Cai
    @Title: Get the sequences name
    @Description: The metatable is downloaded ncbi virus database (refseq and GenBank): https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&Completeness_s=complete.
"""

import sys
import argparse

def main():
	parser = argparse.ArgumentParser(description='Get the sequences reference name')
	parser.add_argument('--filename', default='data/viralGenomeGenbank.txt',help='the ncbi summary file for virus database')
	parser.add_argument('--family', default='Narnaviridae', help='selected virus family name')
	parser.add_argument('--ref', help='only include refseq', action='store_true')

	## read files
	args = parser.parse_args()
	fh = open(args.filename, 'r')
	ACC_COLUMN = 1
	LINEAGE_COLUMN = 5
	TYPE_COLUMN = 7
	TITLE_COLUMN = 19

	## get accession name
	ref_dict = dict()
	new_line = fh.readline()
	while True:
		new_line = fh.readline()
		if len(new_line)==0:
			break
		toks = new_line.strip().split('\t')
		lineage_family = toks[LINEAGE_COLUMN-1]
		if lineage_family.lower()!=args.family.lower():
			continue
		ref = toks[ACC_COLUMN-1]
		if args.ref:
			if toks[TYPE_COLUMN-1]!='RefSeq':
				continue
		ref_title = toks[TITLE_COLUMN-1]
		if (' partial' in ref_title) or (' cds' in ref_title):
			continue
		if ref not in ref_dict.keys():
			ref_dict[ref] = ref_title.strip("\"")
	fh.close()

	for ref, ref_title in ref_dict.items():
		print(ref+'\t'+ref_title)
	return

if __name__=='__main__':
	main()

