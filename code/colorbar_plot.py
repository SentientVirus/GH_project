#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:07:19 2023

@author: marina
"""
import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from Bio import SeqIO
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
import numpy as np
import matplotlib.gridspec as gridspec


pathname = os.path.expanduser('~') + '/GH_project/data/codons'
prefix = ['GH32', 'GH70', 'GS1', 'GS2', 'BRS', 'NGB', 'S1', 'S2a', 'S2b', 'S3']
codeml_file = os.path.expanduser('~') + '/codeml_analysis/core_pairwise_metrics.tsv'

type_color = {'GS1': '#FF0000', 'GS2': '#FF009B', 'BRS': '#D930BD', 'NGB': '#B600FF',
              'S1': '#FF9600', 'S2a': '#FFC800', 'S2b': '#FFC800', 'S3': '#FFEF00'}
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
  GHs = prefix[:2]
  lst = [x for x in lst if x not in GHs] + GHs
  return lst

dict_dict = {}
for pre in prefix:
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
            keys2 = prefix[6:]
        elif key == 'GH70':
            keys2 = prefix[2:6]
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
         'GH32': ['#FF9600', '#FFC800', '#FFC800', '#FFEF00']}
                  
for GH_type in ['GH32', 'GH70']: #dict_dict.keys():
    df = pd.DataFrame.from_dict(dict_dict[GH_type], orient='index', dtype=None, columns=['identity', 'locus1_type', 'locus2_type'])
    for type_GH in prefix[2:]:
        if (type_GH in set(df['locus1_type'])) or (type_GH in set(df['locus2_type'])):
            df = df.replace(type_GH, type_color[type_GH])
            print(f'{type_GH} changed to {type_color[type_GH]}')
    #sns.set(rc={'figure.figsize':(20, 100)})
    fig = plt.figure(constrained_layout=False, figsize=(10, 120))
    spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, height_ratios = [0.01, 1], width_ratios = [0.5, 1])
    cbar_axis = fig.add_subplot(spec[0, :])
    ax1 = fig.add_subplot(spec[1, 0])
    ax2 = fig.add_subplot(spec[1, 1])
    #cbar_ax = fig.add_axes([1.1, .4, .03, .2])
    sns.heatmap(df[['locus1_type', 'locus2_type']], cmap=cmaps[GH_type], 
                cbar=False, ax=ax2, yticklabels=False, xticklabels=False)
    s = sns.heatmap(df['identity'].to_frame(), cmap=cmaps[0], cbar=True, ax=ax1,
                    yticklabels=df.index, cbar_ax=cbar_axis, cbar_kws=dict(use_gridspec = True, orientation='horizontal'))
    
    #colorbar(ax1.get_children()[0], cax=cbar_axis, orientation="horizontal")

    plt.tight_layout(h_pad = 2)
    
    s.set_ylabel('Gene pair', fontsize = 36)
    plt.show()
    #plt.savefig(f'{outdir}/{pre}_dNdSplot.png')
    #print(df)


            
        
#Try to make one of the columns wider

# =============================================================================
# The code so far retrieves the pairwise dS value between strains, but now I
# should add information on gene type as a tuple and include it to a dictionary.
# The tuple types have to be in the same order as the comparison, but the tuples
# can be sorted when I assign colors.
# =============================================================================
