#!/usr/bin/env python

"""
    @Arthor: Jitong Cai
    @Title: Compare gc content
    @Description: Compare structure of wild-type and permutates sequences across viridaes
"""


import subprocess
import numpy as np
import os
import random
import matplotlib.pyplot as plt

def keepCompleteGenome(fastaFile, outFile):
    fh = open(fastaFile, 'r')
    fo = open(outFile, 'w')
    while True:
        newLine = fh.readline()
        if len(newLine)==0:
            break
        if newLine.startswith('>'):
            header = newLine.strip()
            seq = fh.readline().strip()
            if ('complete genome' not in header) and ('complete sequence' not in header):
                continue
            fo.write(header+'\n')
            fo.write(seq+'\n')
    fh.close()
    fo.close()
    return

def getGCprop(seq):
    seq = seq.lower()
    g_cnt = seq.count('g')
    c_cnt = seq.count('c')
    total_cnt = len(seq)
    return(float((g_cnt+c_cnt)/total_cnt))

def statEnergy(energy_list):
    return([np.mean(energy_list), np.std(energy_list), len(energy_list), min(energy_list)])

def get_singleEnergyDist(filePrefix, npermute = 3):
    # wt sequence
    seq = readFasta(filePrefix+'.fasta')
    print('Mfold prediction for '+filePrefix+'.fasta')
    runMfold(filePrefix)
    (seqLen, wtEnergy_list) = collectEnergy(filePrefix)
    if not seqLen:
        return
    gcContent = getGCprop(seq)
    
    # permute sequences
    permutes = []
    for i in range(npermute):
        print('Random permutation '+str(i+1)+'/'+str(npermute))
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
        permutes = permutes + permuteEnergy_list
        # remove random file
        subprocess.run('rm temp*', shell=True)

    ret_list = [seqLen, gcContent] + statEnergy(wtEnergy_list) + statEnergy(permutes)
    return(ret_list)

def get_multiEnergyDist(fasta_dict, outFile):
    fo = open(outFile, 'w')
    for key, fastaFile in fasta_dict.items():
        fh = open(fastaFile, 'r')
        while True:
            newLine = fh.readline()
            if len(newLine)==0:
                break
            if newLine.startswith('>'):
                header = newLine.strip()[1:]
                seq = fh.readline().strip()
                refAcc = header.split(' ')[0]
                fh_single = open(refAcc+'.fasta', 'w')
                fh_single.write('>'+header+'\n')
                fh_single.write(seq+'\n')
                fh_single.close()
                info_list = get_singleEnergyDist(refAcc)
                print_list = [key, header] + [str(x) for x in info_list]
                fo.write('\t'.join(print_list)+'\n')
        fh.close()
    fo.close()
    return

def plot_energyDist(infoFileList):
    # load data
    dat = []
    keyList = []
    for infoFile in infoFileList:
        subprocess.run(' '.join(['cat', infoFile, '>>', 'infoFile.txt']), shell=True)
    fh = open('infoFile.txt', 'r')
    while True:
        newLine = fh.readline()
        if len(newLine)==0:
            break
        toks = newLine.strip().split('\t')
        keyList.append(toks[0])
        dat.append([float(x) for x in toks[2:]])
    # seqLen, gcContent, wtMean, wtSd, wtNumStructure, wtMinEnergy, rdMean, rdSd, rdNumStructure, rdMinEnergy 
    # columns index -> target index; e.g, 1 -> gcContent
    targetIndex = 1
    dat_array = np.array(dat)
    fh.close()

    # plot
    markers = ['o', '^', 's', 'P']
    n = 0
    keyMarker_dict = dict()
    for i in range(len(keyList)):
        key = keyList[i]
        if key not in keyMarker_dict.keys():
            keyMarker_dict[key] = [markers[n], []]
            keyMarker_dict[key][-1].append(i)
            n = n+1
        else:
            keyMarker_dict[key][-1].append(i)
    fig, ax = plt.subplots(figsize=(5,4))
    ax.set_title('Folding Energy to GC Content')
    ax.set_xlabel('Wildtype Folding Energy (kcal/mol)')
    ax.set_ylabel('Random sequence Folding Energy (kcal/mol)')
    
    # scatter point
    cmin, cmax = min(dat_array[:, targetIndex]), max(dat_array[:, targetIndex])
    for key, item in keyMarker_dict.items():
        axscatter = ax.scatter(dat_array[item[-1],2], dat_array[item[-1], 6], c=dat_array[item[-1],targetIndex], marker = keyMarker_dict[key][0], s=50, label=key, cmap='plasma')
        axscatter.set_clim(cmin, cmax)
    ax.legend(loc='upper left')
    fig.colorbar(axscatter, orientation='vertical').set_label('GC comtent')
    # x=y line
    x = np.linspace(*ax.get_xlim())
    ax.plot(x, x, '--', alpha=0.75, c='b')

    plt.savefig('FE_gcContent.png', dpi=300)
    plt.show()
    return

def main():
    fastaFile = '../data/Narnaviridae_refGenome.fasta'
    fastaFile_completeGenome = '../output/Narnaviridae_refGenome_complete.fasta'
    keepCompleteGenome(fastaFile, fastaFile_completeGenome)
    fasta_dict = {'Narna':fastaFile_completeGenome}
    outFile = '../output/Narna_energy.txt'
    get_multiEnergyDist(fasta_dict, outFile)

    fastaFile = '../data/Mitoviridae_refGenome.fasta'
    fastaFile_completeGenome = '../output/Mitoviridae_refGenome_complete.fasta'
    keepCompleteGenome(fastaFile, fastaFile_completeGenome)
    fasta_dict = {'Mito':fastaFile_completeGenome}
    outFile = '../output/Mito_energy.txt'
    get_multiEnergyDist(fasta_dict, outFile)

    #infoFileList = ['Narna_energy.txt', 'Mito_energy.txt'] 
    #plot_energyDist(infoFileList)
    return

if __name__=='__main__':
    main()