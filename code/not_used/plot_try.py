#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:07:19 2023

Earlier version of the heatmap plot that compares each codon position and
plots a different color dependent on if it's identical, synonymous, etc.

@author: Marina Mota Merlo
"""
import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from Bio import SeqIO
import matplotlib.gridspec as gridspec

pathname = os.path.expanduser('~') + '/GH_project/data/codons'
outdir = os.path.expanduser('~') + '/GH_project/plots/pairwise_substitutions'
prefix = ['GH32', 'GH70', 'GS1', 'GS2', 'BRS', 'NGB']

motif_pos = {1: (694, 701), 2: (188, 198), 3: (224, 236), 4: (311, 324), 
             5: (145, 158), 6: (626, 638), 7: (381, 390)}
GH32_motif_pos = {1: (323, 330), 2: (366, 375), 3: (427, 435), 4: (506, 514), 
             5: (565, 572), 6: (621, 630), 7: (649, 658), 8: (791, 797)}
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

cmap_dict = {'GS1': ['white', '#67e3ff', '#0079d8', '#9a9a9a'], 
             'BRS': ['white', '#67e3ff', '#0079d8', '#9a9a9a'],
             'GS2': ['white', '#67e3ff', '#0079d8'], 'NGB': ['white', '#67e3ff', '#0079d8', '#9a9a9a'], 
             'GH70': ['white', '#67e3ff', '#0079d8', '#9a9a9a', '#ff9c45', 
                      '#d64100', 'yellow'], 
             'GH32': ['white', '#67e3ff', '#0079d8', '#9a9a9a', '#ff9c45', 
                      '#d64100', 'yellow']}


if not os.path.exists(outdir):
    os.makedirs(outdir)


def get_codons(sequence):                
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    return codons

def get_catalytic_triad(df):
    positions = [191, 229, 316]
    for column in df:
        for pos in positions:
            if df.loc[pos, column] == 0:
                df.at[pos, column] = 6
                
def get_GH32_catalytic_triad(df):
    positions = [326, 508, 565]
    for column in df:
        for pos in positions:
            if df.loc[pos, column] == 0:
                df.at[pos, column] = 6


def get_motif(start, end, df):
    for column in df:
        for pos in range(start, end+1):
            if df.loc[pos, column] == 1:
                df.at[pos, column] = 4
            elif df.loc[pos, column] == 2:
                df.at[pos, column] = 5

df_aln_dict = {}
for pre in prefix:
    file = f'{pre}_codon.fna'
    
    pos_dict = {}
    with open(f'{pathname}/{file}') as handle:
        f = SeqIO.parse(handle, 'fasta')
        for record in f:
            pos_dict[record.id] = str(record.seq)
    
    pair_dict = {}
    locus_tags = list(pos_dict.keys())
    for i in range(len(locus_tags) - 1):
        for j in range(i + 1, len(locus_tags)):
            pair_dict[(locus_tags[i], locus_tags[j])] = []
            for n in range(0, len(pos_dict[locus_tags[i]]), 3):
                if '-' in pos_dict[locus_tags[i]][n:n+3] or '-' in pos_dict[locus_tags[j]][n:n+3]:
                    pair_dict[(locus_tags[i], locus_tags[j])].append(3)
                elif pos_dict[locus_tags[i]][n:n+3] == pos_dict[locus_tags[j]][n:n+3]:
                    pair_dict[(locus_tags[i], locus_tags[j])].append(0)
                elif transtable[pos_dict[locus_tags[i]][n:n+3]] == transtable[pos_dict[locus_tags[j]][n:n+3]]:
                    pair_dict[(locus_tags[i], locus_tags[j])].append(1)
                else: pair_dict[(locus_tags[i], locus_tags[j])].append(2)
                
    df_aln = pd.DataFrame.from_dict(pair_dict)
    df_aln.index += 1
    if 'GH' in pre:
        if pre == 'GH70':
            for motif in motif_pos.keys():
                get_motif(motif_pos[motif][0], motif_pos[motif][1], df_aln)
            get_catalytic_triad(df_aln)
        elif pre == 'GH32':
            for motif in GH32_motif_pos.keys():
                get_motif(GH32_motif_pos[motif][0], GH32_motif_pos[motif][1], df_aln)
            get_GH32_catalytic_triad(df_aln)
    df_aln_dict[pre] = df_aln



prefix2 = ['GH32', 'GH70', 'GS1', 'GS2', 'BRS', 'NGB', 'S1', 'S2a', 'S2b', 'S3']
codeml_file = os.path.expanduser('~') + '/codeml_analysis/core_pairwise_metrics.tsv'

# type_color = {'GS1': '#FF0000', 'GS2': '#FF009B', 'BRS': '#D930BD', 'NGB': '#B600FF',
#               'S1': '#FF9600', 'S2a': '#FFC800', 'S2b': '#FFC800', 'S3': '#FFEF00'}
type_color = {'GS1': 0, 'GS2': 1, 'BRS': 2, 'NGB': 3,
              'S1': 4, 'S2a': 5, 'S2b': 6, 'S3': 7}

def get_tuple_element(tuple_list, index_no = 0):
    # Function to get the first or second elements in tuples in a list
    try:
        element_no = [tup[index_no] for tup in tuple_list]
        return element_no
    except IndexError:
        print('Index is not a tuple index (0 or 1)!')
    except TypeError:
        print('Input is not a list of tuples!')
        
def move_GHs(lst):
  # Moves GH32 and GH70 categories to the end
  GHs = prefix2[:2]
  lst = [x for x in lst if x not in GHs] + GHs
  return lst

dict_dict = {}
for pre in prefix2:
    file = f'{pre}_codon.fna'
    
    pos_dict = {}
    with open(f'{pathname}/{file}') as handle:
        f = SeqIO.parse(handle, 'fasta')
        locus_tags = [record.id for record in f]
    
    pair_dict = {}
    with open(codeml_file) as codeml:
        codeml_df = pd.read_csv(codeml, sep = '\t')
    for i in range(len(locus_tags) - 1):
        strain1 = locus_tags[i].split('_')[0]
        if strain1 == 'K2W83':
            s1 = 'H3B104J'
        elif strain1 == 'APS55':
            strain1 = 'MP2'
            s1 = 'MP2'
        elif strain1 == 'LDX55':
            s1 = 'H4B206J'
        elif strain1 == 'FHON2':
            s1 = 'fhon2'
        else: s1 = strain1
        if s1.startswith('H') and 'B' in s1:
            s1 = s1[:4] + '-' + s1[4:] 
        for j in range(i + 1, len(locus_tags)):
            strain2 = locus_tags[j].split('_')[0]
            if strain2 == 'K2W83':
                s2 = 'H3B104J'
            elif strain2 == 'APS55':
                strain2 = 'MP2'
                s2 = 'MP2'
            elif strain2 == 'LDX55':
                s2 = 'H4B206J'
            elif strain2 == 'FHON2':
                s2 = 'fhon2'
            else: s2 = strain2
            if s2.startswith('H') and 'B' in s2:
                s2 = s2[:4] + '-' + s2[4:]
            # print(s1, s2)
            # print(strain1, strain2)
            if s1 == s2:
                pair_dict[(locus_tags[i], locus_tags[j])] = [0]
            else:
                pair_dict[(locus_tags[i], locus_tags[j])] = [float(codeml_df[((codeml_df['strain1'].replace('-', '') == s1) & (codeml_df['strain2'].replace('-', '') == s2)) | ((codeml_df['strain2'].replace('-', '') == s1) & (codeml_df['strain1'].replace('-', '') == s2))]['mean_dS'])]
    dict_dict[pre] = pair_dict
    
loctag_dict = {}
keys = list(dict_dict.keys())
keys = move_GHs(keys)
for key in keys:
    loctag_dict[key] = []
    if not key.startswith('GH'):
        for pair in list(dict_dict[key].keys()):
            dict_dict[key][pair].extend([key, key])
            [loctag_dict[key].append(element) for element in pair if element not in loctag_dict[key]]
    else:
        if key == 'GH32':
            keys2 = prefix2[6:]
        elif key == 'GH70':
            keys2 = prefix2[2:6]
        for key2 in keys2:
            for pair in list(dict_dict[key].keys()):
                if len(dict_dict[key][pair]) < 3:
                    dict_dict[key][pair].extend(['', ''])
                if pair[0] in loctag_dict[key2]:
                    dict_dict[key][pair][1] = key2
                if pair[1] in loctag_dict[key2]:
                    dict_dict[key][pair][2] = key2
                    
cmaps = {0: sns.diverging_palette(250, 30, l=65, center="light", as_cmap=True), 
         'GH70': ['#FF0000', '#FF009B', '#D930BD', '#B600FF'], 
         'GH32': ['#FF9600', '#FFC800', '#FFB600', '#FFEF00']}
                  
for GH_type in prefix: #dict_dict.keys():
    if GH_type.startswith('GH'):
        df = pd.DataFrame.from_dict(dict_dict[GH_type], orient='index', dtype=None, columns=['Identity', 'Type 1', 'Type 2'])
        df1 = df.copy()
        df1.index = [f'{tup[0]}/{tup[1]}' for tup in df1.index]
    
        for type_GH in prefix2[2:]:
            if (type_GH in set(df['Type 1'])) or (type_GH in set(df['Type 2'])):
                df1 = df1.replace(type_GH, type_color[type_GH])
                print(f'{type_GH} changed to {type_color[type_GH]}')
        if GH_type == 'GH32':
            h = 120
        else: h = 120       
        fig = plt.figure(constrained_layout=False, figsize=(120, h))
        spec = gridspec.GridSpec(ncols=3, nrows=2, figure=fig, height_ratios = [0.002, 1], width_ratios = [0.05, 0.1, 10])
        cbar_axis = fig.add_subplot(spec[0, :2])
        ax1 = fig.add_subplot(spec[1, 0])
        ax2 = fig.add_subplot(spec[1, 1])
        ax3 = fig.add_subplot(spec[1, 2])
        h1 = sns.heatmap(df1[['Type 1', 'Type 2']], cmap=cmaps[GH_type], 
                    cbar=False, ax=ax2, yticklabels=False, xticklabels=True)
        h2 = sns.heatmap(df_aln_dict[GH_type].T.reset_index(drop=True), cmap = cmap_dict[GH_type], cbar = False, ax=ax3, yticklabels=False)
        h3 = sns.heatmap(df1['Identity'].to_frame(), cmap=cmaps[0], cbar=True, ax=ax1,
                        yticklabels=df1.index, cbar_ax=cbar_axis, 
                        cbar_kws=dict(use_gridspec = True, orientation='horizontal'),
                        vmin = 0, vmax = 1)
        print(f'Made {GH_type} plot!')
    else:
        fig = plt.figure(constrained_layout=False, figsize=(120, 120))
        spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
        ax = fig.add_subplot(spec[0, 0])
        h3 = sns.heatmap(df_aln_dict[GH_type].T, cmap = cmap_dict[GH_type], 
                         cbar=False, ax=ax)
        print(f'Made {GH_type} plot!')
    # Remove None-None label, decrease space between columns and change the format in which pairs are written.

    plt.tight_layout(h_pad = 2, w_pad = -0.5)
    
    if GH_type == 'GH70':
        ax1.yaxis.set_major_locator(plt.MaxNLocator(200))
    
    h3.set_ylabel('Gene pair', fontsize = 36)
    h3.set_yticklabels(h3.get_ymajorticklabels(), fontsize = 10) 
    #plt.show()
    plt.savefig(f'{outdir}/{GH_type}_dNdSplot.png')
    print(f'{GH_type} plot saved!')

#Check that I am using the correct df_aln
            


# =============================================================================
# The code so far retrieves the pairwise dS value between strains, but now I
# should add information on gene type as a tuple and include it to a dictionary.
# The tuple types have to be in the same order as the comparison, but the tuples
# can be sorted when I assign colors.
# =============================================================================
