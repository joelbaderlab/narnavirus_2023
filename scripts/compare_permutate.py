#!/usr/bin/env python

"""
    @Arthor: Jitong Cai
    @Title: Compare permutate
    @Description: Compare structure of wild-type and permutates sequences
"""

import random
import subprocess
import argparse
import sys
import os

def permuteSeq(seq):
    seqList = list(seq)
    random.shuffle(seqList)
    permuteSeq = ''.join(seqList)
    return(permuteSeq)

def extractEnergy(ctFile):
    '''
    Modified from parse_ct, only read the first line
    Only one seq should be contained in the ctFile
    '''
    fh = open(ctFile)
    header = fh.readline().strip()
    toks = header.split('dG = ')
    fh.close()
    seqLen = int(toks[0].strip())
    energy = float(toks[1].strip().split(' ')[0])
    return((seqLen, energy))
    
def collectEnergy(filePrefix):
    files = os.listdir()
    selectFiles = [x for x in files if x.startswith(filePrefix+'_') and x.endswith('.ct') ]
    energy_list = []
    if len(selectFiles)==0:
        sys.stderr.write('MFOLD makes no predictions\n')
        return((False, False))
    for ctfile in selectFiles:
        (seqLen, energy) = extractEnergy(ctfile)
        energy_list.append(energy)
    return((seqLen, energy_list))

def readFasta(file):
    '''
    only one seq should be contained in the fasta file
    '''
    fh = open(file, 'r')
    while True:
        newLine = fh.readline()
        if len(newLine)==0:
            break
        if newLine.startswith('>'):
            seq = ''
            continue
        seq = seq + newLine.strip()
    fh.close()
    return(seq)

def runMfold(filePrefix, toolPath='mfold'):
    cmd = [toolPath, 'T=30', 'SEQ='+filePrefix+'.fasta']
    subprocess.run(' '.join(cmd),shell=True)
    return

def compareRandomEnergy(filePrefix, npermute = 20):
    # wt sequence
    seq = readFasta(filePrefix+'.fasta')
    sys.stderr.write('Mfold prediction for '+filePrefix+'.fasta'+'\n')
    runMfold(filePrefix)
    (seqLen, wtEnergy_list) = collectEnergy(filePrefix)
    if not seqLen:
        return
    
    # permute sequences
    permutes = []
    for i in range(npermute):
        sys.stderr.write('Random permutation '+str(i+1)+'/'+str(npermute)+'\n')
        # write permuted sequence as a fastafile
        tempSeq = permuteSeq(seq)
        fh = open('temp.fasta', 'w')
        fh.write('>randomSeq'+str(i)+'\n')
        fh.write(tempSeq+'\n')
        fh.close()
        # predict by mfold
        runMfold('temp')
        # collect energy
        (seqLen, permuteEnergy_list) = collectEnergy('temp')
        permutes.append(permuteEnergy_list)
        # remove random file
        subprocess.run('rm temp*', shell=True)
    return([seqLen, wtEnergy_list, permutes])

def plotDat(seqLen, wtEnergy_list, permutes, outfile):
    fh = open(outfile, 'w')
    label = 'Wildtype'
    for i in range(len(wtEnergy_list)):
        fh.write(label+'\t'+str(wtEnergy_list[i])+'\t'+str(seqLen)+'\n')
    for i in range(len(permutes)):
        label = 'Random_'+str(i+1)
        for energy in permutes[i]:
            fh.write(label+'\t'+str(energy)+'\t'+str(seqLen)+'\n')
    fh.close()
    return

def main():
    parser = argparse.ArgumentParser(description = 'This program permutes given sequence and compare the folding energy')
    parser.add_argument('--inputFasta', help='input fasta file with suffix as .fasta or .fa')
    parser.add_argument('--mfoldOut', help='flag to keep mfold outputs', action='store_true')
    parser.add_argument('--outdir', help='output directory', default='./')
    args = parser.parse_args()

    sys.stderr.write('File for analysis: %s \n' % args.inputFasta)
    sys.stderr.write('Output directory: %s \n' % args.outdir)
    if (not os.path.isdir(args.outdir)):
        sys.stderr.write('creating outputs directory %s \n' % args.outdir)
        os.makedirs(args.outdir)

    # structure prediction
    if args.inputFasta.endswith('.fasta'):
        outPrefix = args.inputFasta[:-6].split('/')[-1]
    elif args.inputFasta.endswith('.fa'):
        outPrefix = args.inputFasta[:-3].split('/')[-1]
    else:
        sys.stderr.write('inputFasta has wrong format: %s \n' % args.inputFasta)
        return

    cmd = ['cp', args.inputFasta, 'mfold_'+outPrefix+'.fasta']
    subprocess.run(' '.join(cmd), shell=True)
    [seqLen, wtEnergy_list, permutes] = compareRandomEnergy('mfold_'+outPrefix)
    plotDat(seqLen, wtEnergy_list, permutes, args.outdir+outPrefix+'.dat')

    if args.mfoldOut:
        if (not os.path.isdir(args.outdir+outPrefix)):
            sys.stderr.write('creating mfold outputs directory %s \n'  % args.outdir+outPrefix)
            os.makedirs(args.outdir+outPrefix)
        cmd = ['mv', 'mfold_'+outPrefix+'*', args.outdir+outPrefix]
        subprocess.run(' '.join(cmd), shell=True)
    else:
        cmd = ['rm', 'mfold_'+outPrefix+'*']
        subprocess.run(' '.join(cmd), shell=True)
    return

if __name__=='__main__':
    main()

