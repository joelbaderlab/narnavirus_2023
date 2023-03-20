#!/usr/bin/env python

"""
    @Arthor: Jitong Cai
    @Title: Generate contact map
    @Description: Generate contact map for a given structure for visualization
"""

import sys
import os
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors

## mfold would provide no information about folding probability, here use RNAfold to capture the folding probability information

def runRNAfold(filePrefix, rnafoldPath='../tools/ViennaRNA-2.5.0/src/bin/RNAfold', utilsPath='../tools/ViennaRNA-2.5.0/src/Utils'):
    cmd = [rnafoldPath, '--noLP', '-d2', '-p', '--id-prefix='+filePrefix, '-i', filePrefix+'.fasta', '>', filePrefix+'.out']
    subprocess.run(' '.join(cmd), shell=True)
    cmd = [utilsPath+'/relplot.pl', filePrefix+'_0001_ss.ps', filePrefix+'_0001_dp.ps', '>', filePrefix+'.ps']
    subprocess.run(' '.join(cmd), shell=True)
    cmd = [utilsPath+'/mountain.pl', filePrefix+'_0001_dp.ps', '>', filePrefix+'.mountain']
    subprocess.run(' '.join(cmd), shell=True)
    return


def parseRNAfold_prob(filePrefix):
    fh = open(filePrefix+'.out', 'r')
    fh.readline()
    seq = fh.readline().strip()
    dot = fh.readline().strip().split(' ')[0]
    fh.close()

    fh = open(filePrefix+'_0001_dp.ps')
    probDict = dict()
    probRecord = []
    for i in range(len(dot)):
        probDict[i] = [0]
    while True:
        newLine = fh.readline()
        if newLine.startswith('%start of base pair probability data'):
            while True:
                newLine = fh.readline()
                if len(newLine) == 0:
                    break
                toks = newLine.strip().split(' ')
                if toks[-1]!='ubox':
                    break
                pos1 = int(toks[0])-1
                pos2 = int(toks[1])-1
                prob = float(toks[2])
                probDict[pos1].append(prob)
                probDict[pos2].append(prob)
                probRecord.append((pos1, pos2, prob))
        if len(newLine)==0:
            break
    fh.close()

    probList = []
    for i in range(len(dot)):
        if dot[i]=='.':
            probList.append(1-max(probDict[i]))
        else:
            probList.append(max(probDict[i]))
    
    #probMat = np.zeros((len(seq), len(seq)))
    #for (pos1, pos2, prob) in probRecord:
    #    probMat[pos1, pos1:(pos2+1)] = prob
    #    probMat[pos1:(pos2+1), pos2] = prob

    return([seq, dot, probList, probRecord])

def parseRNAfold_entropy(filePrefix):
    '''
    Could be a replacement for parseRNAfold_prob
    Calculate the positional entropy for each position. Well-defined regions are identified by low entropy.
    '''
    fh = open(filePrefix+'.out', 'r')
    fh.readline()
    seq = fh.readline().strip()
    dot = fh.readline().strip().split(' ')[0]
    fh.close()

    fh = open(filePrefix+'.ps', 'r')
    probList = []
    while True:
        newLine = fh.readline()
        if newLine.startswith('/pairs ['):
            while True:
                newLine = fh.readline()
                if len(newLine) == 0:
                    break
                if newLine.startswith('] def'):
                    break
                toks = newLine.strip().strip('[]').split(' ')
                left, right = int(toks[0])-1, int(toks[1])-1
                if not (dot[left]=='(' and dot[right]==')'):
                    sys.stderr.write('probability and structure don\'t match\n')
        if newLine.startswith('/S ['):
            while True:
                newLine = fh.readline()
                if len(newLine) == 0:
                    break
                if newLine.startswith('] def'):
                    break
                probList.append(newLine.strip())
        if len(newLine) == 0:
            break
    fh.close()
    return([seq, dot, probList])

def dot2pairedDict(dotSeq):
    '''modified from dot2ct'''
    bracket_dict = {'(':[],'[':[],'{':[],'<':[]}
    corresponse_dict = {')':'(',']':'[','}':'{','>':'<'}
    paired_dict = dict()
    for i in range(len(dotSeq)):
        if dotSeq[i] in bracket_dict.keys():
            bracket_dict[dotSeq[i]].append(i)
        elif dotSeq[i] in corresponse_dict.keys():
            pairedIndex = bracket_dict[corresponse_dict[dotSeq[i]]].pop()
            paired_dict[pairedIndex]=i
            paired_dict[i] = pairedIndex
        else:
            paired_dict[i]=-1
    return(paired_dict)

def dot2dist(dot):
    distList = []
    paired_dict = dot2pairedDict(dot)
    for i in range(len(dot)):
        if paired_dict[i] == -1:
            distList.append(-1)
        else:
            distList.append(abs(paired_dict[i]-i))
    return(distList)

def main():

    parser = argparse.ArgumentParser(description = 'This program generates contact map for RNAfold structure prediction')
    parser.add_argument('--inputFasta', help='input fasta file with suffix as .fasta or .fa')
    parser.add_argument('--rnafoldOut', help='flag to keep mfold outputs', action='store_true')
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
    
    fastaFilePrefix = 'rnafold_'+outPrefix
    cmd = ['cp', args.inputFasta, 'rnafold_'+outPrefix+'.fasta']
    subprocess.run(' '.join(cmd), shell=True)
    runRNAfold(fastaFilePrefix)
    [seq, dot, scores, probRecord] = parseRNAfold_prob(fastaFilePrefix)

    # setup the normalization and the colormap
    probRecord.sort(key=lambda x: x[-1])
    probs = np.array([x[-1] for x in probRecord])
    normalize = mcolors.Normalize(vmin=probs.min(), vmax=probs.max())
    colormap = cm.Blues

    # plot
    fig, ax = plt.subplots(figsize=(15,15)) 
    for (pos1, pos2, prob) in probRecord:
        ax.plot([pos1, pos1], [pos1, pos2], c=colormap(prob), alpha=0.5)
        ax.plot([pos1, pos2], [pos2, pos2], c=colormap(prob), alpha=0.5)

    # setup the colorbar
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(probs)
    #fig.colorbar(scalarmappaple, orientation="horizontal",aspect=50)

    # show the figure
    plt.savefig(args.outdir+outPrefix+'.pdf', dpi=400)

    if args.rnafoldOut:
        if (not os.path.isdir(args.outdir+outPrefix)):
            sys.stderr.write('creating rnafold outputs directory %s \n'  % (args.outdir+outPrefix))
            os.makedirs(args.outdir+outPrefix)
        cmd = ['mv', 'rnafold_'+outPrefix+'*', args.outdir+outPrefix]
        subprocess.run(' '.join(cmd), shell=True)
    else:
        cmd = ['rm', 'rnafold_'+outPrefix+'*']
        subprocess.run(' '.join(cmd), shell=True)

    return

if __name__=='__main__':
    main()