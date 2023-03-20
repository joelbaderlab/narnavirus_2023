#!/usr/bin/env python

"""
    @Arthor: Jitong Cai
    @Title: Structure segment
    @Description: Predict and average across structure ensemble.
                  Segment the sequence into regions with independent secondary structures.
"""

import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess
import sys
import argparse

def parse_ct(file):
    fh = open(file, 'r')
    header = fh.readline().strip()
    toks = header.split('dG = ')
    seqLen = int(toks[0].strip())
    energy = float(toks[1].strip().split(' ')[0])
    pairMatrix = np.zeros((seqLen, seqLen))
    while True:
        newLine = fh.readline()
        if len(newLine)==0:
            break
        toks = newLine.split(' ')
        toks = [x for x in toks if x!='']
        baseIndex = int(toks[0])
        pairIndex = int(toks[4])
        if pairIndex==0:
            continue
        if baseIndex>pairIndex:
            continue
        pairMatrix[baseIndex-1, pairIndex-1] = 1
    return(energy, pairMatrix)

def add_dispersion(pairMatrix, baseScore):
    seqLen = pairMatrix.shape[0]
    disperseMatrix = np.zeros((seqLen, seqLen))
    for i in range(seqLen):
        for j in range(seqLen):
            if pairMatrix[i, j]==1:
                disperseMatrix[i, j]=baseScore
                boundBaseScore = baseScore
                bound = 0
                while True:
                    bound = bound+1
                    boundBaseScore = boundBaseScore * 1/2
                    if boundBaseScore < 1:
                        break
                disperseMatrix[(i-bound):(i+bound+1), (j-bound):(j+bound+1)] += boundBaseScore
    return(disperseMatrix)

def detect_stem(scoreMatrix, thresh):
    indexRow, indexCol = np.where(scoreMatrix>=thresh)
    indexDot = list(zip(indexRow, indexCol))
    indexNum = indexRow-indexCol
    indexDot = [x for _,x in sorted(zip(indexNum, indexDot))]
    seqLen = scoreMatrix.shape[0]
    stemDict = dict()
    for num in range(seqLen-1, 0, -1):
        if len(indexDot)==0:
            break
        
        num2dots = []
        while True:
            if len(indexDot)==0:
                break
            newDot = indexDot.pop()
            if newDot[0]-newDot[1] < num:
                indexDot.append(newDot)
                break
            if newDot[0]-newDot[1] == num:
                num2dots.append(newDot)
        if len(num2dots)==0:
            continue
        num2dots.sort()

        num2stems = []
        for stem in stemDict.keys():
            if stemDict[stem][0] == num+1:
                num2stems.append(stem)
        
        for newDot in num2dots: 
            assign = False
            for stem in num2stems:
                farEnd = stemDict[stem][1]
                if (newDot[0]>=farEnd[0]-2) and (newDot[0]<=farEnd[1]+2):
                    assign = True
                    if stemDict[stem][0]-num==1:
                        stemDict[stem][0] = num
                        stemDict[stem][1] = (newDot[0], newDot[0])
                        stemDict[stem][2].append(newDot)
                    elif num-stemDict[stem][0]==0:
                        stemDict[stem][1] = (min(newDot[0], farEnd[0]), max(newDot[0], farEnd[1]))
                        stemDict[stem][2].append(newDot)
                    break
            if not assign:
                stemDict[newDot] = [num, (newDot[0], newDot[0]), [newDot]]

    squareDict = dict()
    for stem in stemDict.keys():
        (x, y) = stem
        (a, b) = stemDict[stem][1][0], stemDict[stem][1][0]-(stemDict[stem][0])
        dots = stemDict[stem][2]
        for dot in dots:
            if dot[0]<a:
                a = dot[0]
            if dot[0]>x:
                x = dot[0]
            if dot[1]>b:
                b = dot[1]
            if dot[1]<y:
                y = dot[1]
        squareDict[stem] = (x, y, a, b)
    return(squareDict)

def plot_structure(plotMatrix, outPrefix, stemDict=dict(), subregion=False):
    fig, ax = plt.subplots(figsize=(10,8))
    if subregion:
        im = ax.imshow(plotMatrix, vmin=0, vmax=1, cmap='Blues', interpolation='nearest', origin='lower', alpha=0.5, extent=[subregion[0], subregion[1], subregion[0], subregion[1]])
    else:  
        im = ax.imshow(plotMatrix, vmin=0, vmax=1, cmap='Blues', interpolation='nearest', origin='lower', alpha=0.5)
        subregion = [0,0]
    ax.set_xlabel('Distance from 5\' end')
    ax.set_ylabel('Distance from 5\' end')
    if len(stemDict)!=0:
        for stem, bound in stemDict.items():
            width = bound[2]-bound[0]
            height = bound[3]-bound[1]
            rect = patches.Rectangle((subregion[0]+bound[0],subregion[0]+bound[1]), width, height,linewidth=1,edgecolor='b',facecolor='none')
            ax.add_patch(rect)
    plt.savefig(outPrefix+'_dotplot.png', dpi=300)
    return
    
def detect_region(stemDict, dist=200):
    boundList = [(y,x) for (x,y) in stemDict.keys()]
    boundList.sort()
    (anchorStart, anchorEnd) = (0,0)
    assign = False
    regionList = []
    longrangeList = []
    for bound in boundList:
        (start, end) = bound
        if end-start>dist:
            longrangeList.append((start, end))
            continue
        if start > anchorEnd:
            if assign:
                regionList.append((anchorStart, anchorEnd))
            anchorStart = start
            anchorEnd = end
        elif start <= anchorEnd:
            anchorEnd = max(end, anchorEnd)
        assign = True
    if assign:
        regionList.append((anchorStart, anchorEnd))
    return(regionList, longrangeList)
            
def plot_region(plotMatrix, regionList, outPrefix, subregion=False):
    fig, ax = plt.subplots(figsize=(10,8))
    if subregion:
        im = ax.imshow(plotMatrix, vmin=0, vmax=1, cmap='Blues', interpolation='nearest', origin='lower', alpha=0.5, extent=[subregion[0], subregion[1], subregion[0], subregion[1]])
    else:
        im = ax.imshow(plotMatrix, vmin=0, vmax=1, cmap='Blues', interpolation='nearest', origin='lower', alpha=0.5)
        subregion = [0, 0]
    ax.set_xlabel('Distance from 5\' end')
    ax.set_ylabel('Distance from 5\' end')
    if len(regionList)!=0:
        for (start, end) in regionList:
            width = start - end
            height = end - start
            rect = patches.Rectangle((subregion[0]+end, subregion[0]+start), width, height,linewidth=1,edgecolor='r',facecolor='none')
            ax.add_patch(rect)
    plt.savefig(outPrefix+'_regions.png', dpi=300)
    return 

def segment_structure(ctPrefix, outPrefix, baseScore):
    files = os.listdir()
    selectFiles = [x for x in files if x.startswith(ctPrefix+'_') and x.endswith('.ct') ]
    matrixList = []
    for ctfile in selectFiles:
        energy, pairMatrix = parse_ct(ctfile)
        disperseMatrix = add_dispersion(pairMatrix, baseScore)
        matrixList.append(disperseMatrix)
    metaMatrix = np.mean(np.stack(matrixList), axis=0)
    stemDict = detect_stem(np.transpose(metaMatrix), 0.2*baseScore)

    metaMatrix = metaMatrix/baseScore
    plot_structure(metaMatrix, outPrefix, stemDict)

    regionList, longrangeList = detect_region(stemDict)
    plot_region(metaMatrix, regionList, outPrefix)
    sys.stdout.write('Region List:\n')
    sys.stdout.write('\t'.join([str(x) for x in regionList])+'\n')

    sys.stdout.write('Long-range interaction:\n')
    sys.stdout.write('\t'.join([str(x) for x in longrangeList])+'\n')
    return(metaMatrix)

def main():
    parser = argparse.ArgumentParser(description = 'This program predict secondary structure for given sequence and segment the structure into subregions')
    parser.add_argument('--inputFasta', help='input fasta file with suffix as .fasta or .fa')
    parser.add_argument('--mfold', help='mfold tool path', default='mfold', required=False)
    parser.add_argument('--mfoldOut', help='flag to keep mfold outputs', action='store_true')
    parser.add_argument('--smooth', help='base score for smoothing, smoother with larger value', type=float, default=5.0, required=False)
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
    
    sys.stderr.write('Start mfold secondary structure prediction...\n')
    cmd = ['cp', args.inputFasta, 'mfold_'+outPrefix+'.fasta']
    subprocess.run(' '.join(cmd), shell=True)
    cmd = [args.mfold, 'SEQ='+'mfold_'+outPrefix+'.fasta', 'T=30']
    subprocess.run(' '.join(cmd), shell=True)
    sys.stderr.write('mfold command:'+' '.join(cmd)+'\n')

    # structure segmentation
    print('Start structure segmentation...')
    segment_structure('mfold_'+outPrefix, outPrefix, args.smooth)
    cmd = ['mv', outPrefix+'*'+'.png', args.outdir]
    subprocess.run(' '.join(cmd), shell=True)

    if args.mfoldOut:
        if (not os.path.isdir(args.outdir+outPrefix)):
            sys.stderr.write('creating mfold outputs directory %s \n'  % (args.outdir+outPrefix))
            os.makedirs(args.outdir+outPrefix)
        cmd = ['mv', 'mfold_'+outPrefix+'*', args.outdir+outPrefix]
        subprocess.run(' '.join(cmd), shell=True)
    else:
        cmd = ['rm', 'mfold_'+outPrefix+'*']
        subprocess.run(' '.join(cmd), shell=True)
    return
    
if __name__=='__main__':
    main()

