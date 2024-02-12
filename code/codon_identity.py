#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:20:09 2022

@author: marina
"""

# =============================================================================
# This script reads a fasta file with nucleotide information and calculates the
# percentage of identity in the third codon position between every pair of
# sequences in the file.
# =============================================================================

import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from Bio import AlignIO
from Bio.Seq import Seq
from statistics import mean

pathname = os.path.expanduser('~') + '/GH_project/data/codons'
prefix = ['GH70', 'GS1', 'GS2', 'BRS']
#pathname = os.path.expanduser('~') + '/codeml_analysis/subset/inputs'

transtable_RNA = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L',
              'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I',
              'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
              'GUG': 'V', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
              'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACU': 'T',
              'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A',
              'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop',
              'UAG': 'Stop', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 
              'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C',
              'UGA': 'Stop', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R',
              'CGG': 'R', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

transtable = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L',
              'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I',
              'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
              'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T',
              'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
              'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Stop',
              'TAG': 'Stop', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 
              'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C',
              'TGA': 'Stop', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
              'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

groups = {'A': 'ATU', 'T': 'ATU', 'U': 'ATU', 'C': 'CG', 'G': 'GC'}


def get_codons(sequence):                
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    return codons

for pre in prefix:
    file = f'{pre}_codon.pal2nal'
    
    pos_dict = {}
    with open(f'{pathname}/{file}') as handle:
        f = (line for line in handle)
        for line in f:
            if '_' in line:
                id_name = line.strip()
                pos_dict[id_name] = ''
            elif not any(char.isdigit() for char in line):
                pos_dict[id_name] += line.strip()
    
    pair_dict = {}
    locus_tags = list(pos_dict.keys())
    for i in range(len(locus_tags) - 1):
        for j in range(i + 1, len(locus_tags)):
            pair_dict[(locus_tags[i], locus_tags[j])] = []
            for n in range(0, len(pos_dict[locus_tags[i]]), 3):
                if pos_dict[locus_tags[i]][n:n+3] == pos_dict[locus_tags[j]][n:n+3]:
                    pair_dict[(locus_tags[i], locus_tags[j])].append(0)
                elif transtable[pos_dict[locus_tags[i]][n:n+3]] == transtable[pos_dict[locus_tags[j]][n:n+3]]:
                    pair_dict[(locus_tags[i], locus_tags[j])].append(1)
                else: pair_dict[(locus_tags[i], locus_tags[j])].append(2)
                
    df = pd.DataFrame.from_dict(pair_dict)
    if pre == 'GS1':
        y = 100
    elif pre == 'GS2':
        y = 10
    else: y = 25
    sns.set(rc={'figure.figsize':(80, y)})
    sns.heatmap(df.T, cmap=['white', '#88e5c3', '#1d955e'], cbar=False)
    plt.savefig(f'{pre}_dNdSplot.png')

            
    # for record in AlignIO.PhylipIO.SequentialPhylipIterator(handle):
    #     print('Loop!')
    #     posstr = ''
    #     for nuc in record.seq:
    #         if nuc != '-':
    #             posstr += nuc
    #         record.seq = Seq(posstr.upper()).transcribe()
    #         pos_dict[record.id] = get_codons(str(record.seq))

# percent_dict = {}                
# locus_tags = sorted(list(pos_dict.keys()))
# for n in range(len(locus_tags)-1):
#     codons_n = pos_dict[locus_tags[n]]
#     if len(codons_n[-1:]) != 3:
#         codons_n = codons_n[:-1]
#     for m in range(n+1, len(locus_tags)):
#         codons_m = pos_dict[locus_tags[m]]
#         if len(codons_m[-1:]) != 3:
#             codons_m = codons_m[:-1]
#         minlen = min(len(codons_n), len(codons_m))
#         matches = [(transtable[codons_n[i]] == transtable[codons_m[i]] and codons_n[i] != codons_m[i]) for i in range(minlen)]
#         STS = [(transtable[codons_n[i]] == transtable[codons_m[i]] and codons_n[i] != codons_m[i]) for i in range(minlen) if groups[codons_n[i][2]] == groups[codons_m[i][2]]]
#         STV = [(transtable[codons_n[i]] == transtable[codons_m[i]] and codons_n[i] != codons_m[i]) for i in range(minlen) if groups[codons_n[i][2]] != groups[codons_m[i][2]]]
#         mismatches = [transtable[codons_n[i]] != transtable[codons_m[i]] for i in range(minlen)]
#         eqsites = [(transtable[codons_n[j]] == transtable[codons_m[j]]) for j in range(minlen)]
#         samesites = [codons_n[j] == codons_m[j] for j in range(minlen)]
#         difsites = [codons_n[j] != codons_m[j] for j in range(minlen)]
#         percent = sum(matches)/minlen*100
#         percent2 = sum(mismatches)/minlen*100
#         dS = sum(matches)
#         # percent2 = sum(mismatches)/(minlen - sum(eqsites))
#         percent_dict[(locus_tags[n], locus_tags[m])] = [sum(eqsites), sum(difsites), minlen] #[percent, percent2, sum(STS)+sum(STV), sum(mismatches), sum(samesites), sum(difsites), minlen]
        


            




#             pos = 0
#             pos3str, pos2str, pos1str = '', '', ''
#             for nuc in record.seq:
#                 if nuc != '-':
#                     pos += 1
#                     if pos % 3 == 0:
#                         pos3str += nuc
#                     elif pos % 3 == 1:
#                         pos1str += nuc
#                     elif pos % 3 == 2:
#                         pos2str += nuc
#             pos_dict[record.id] = [pos1str, pos2str, pos3str]
            
# percent_dict = {}
# for locus1 in sorted(list(pos_dict.keys())):
#     for locus2 in pos_dict.keys():
#         if locus1 != locus2 and (locus2, locus1) not in list(percent_dict.keys()):
#             percent = []
#             minmax = []
#             for j in range(3):
#                 minlen = min(len(pos_dict[locus1][j]), len(pos_dict[locus2][j]))
#                 maxlen = max(len(pos_dict[locus1][j]), len(pos_dict[locus2][j]))
#                 matches = [pos_dict[locus1][j][i] == pos_dict[locus2][j][i] for i in range(minlen)]
#                 percent.append(sum(matches)/len(matches)*100)
#                 minmax.append((minlen, maxlen))
#             percent_dict[(locus1, locus2)] = percent + minmax
        