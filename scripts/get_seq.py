#!/usr/bin/env python

"""
    @Arthor: Jitong Cai
    @Title: Get the whole sequences
    @Description: Read the accession list and search the sequence in fasta file. Preparation for MSA.
"""

import sys
import argparse

def parse_seqlist(seqlist):
    """
        get the dictionary of sequence accessions
    """
    fh = open(seqlist,'r')
    ret = []
    while True:
        new_line = fh.readline()
        if len(new_line)==0:
            break
        toks = new_line.strip().split('\t')
        ref_name = toks[0]
        ret.append(ref_name)
    return(ret)


def main():
    parser = argparse.ArgumentParser(description='Get the start and end sequence from the database')
    parser.add_argument('--filename', default='NCBIviralGenome/viralGenomeGenbank.fasta',help='the input fasta file')
    parser.add_argument('--seqlist', default='Narnaviridae_list.txt',help='the sequence reference name list')

    ## read files
    args = parser.parse_args()
    fh = open(args.filename, 'r')
    if args.seqlist:
        ref_list = parse_seqlist(args.seqlist)

    ## parse fasta
    ## first line in fasta file
    new_line = fh.readline()
    if new_line.startswith('>'):
        name = new_line[1:].rstrip()
        seq = ''
    else:
        sys.exit('The fasta file is not in correct format')
    while True:
        new_line = fh.readline()
        ## take actions after all reading is done
        if len(new_line)==0:
            if args.seqlist:
                if any(ref_name in name for ref_name in ref_list):
                    print('>'+name)
                    print(seq)
            else:
                print('>'+name)
                print(seq)
            break
        ## take actions after detect a new sequence
        elif new_line.startswith('>'):
            if args.seqlist:
                if any(ref_name in name for ref_name in ref_list):
                    print('>'+name)
                    print(seq)
            else:
                print('>'+name)
                print(seq)
            name = new_line[1:].rstrip()
            seq = ''
        else:
            seq = seq + new_line.rstrip()

    fh.close()
    return

if __name__=='__main__':
    main()