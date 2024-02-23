#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 16:54:02 2023

This is a plot to create matrices for genes GS1-2, BRS and NGB, where
the dS pairwise values are plotted on the top half and the dN pairwise
values are plotted on the bottom half.

@author: Marina Mota Merlo
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib as mpl

genes = ['GS1', 'GS2', 'BRS', 'NGB']

# GS1 = ['A1001_12310', 'A1003_12540', 'A1202_13520', 'A1401_12750', 'A1805_12820',
#        'FHON2_13540', 'G0101_12800', 'G0403_13100', 'H1B104J_13010', 
#        'H1B105A_12300', 'H1B302M_12880', 'H3B101A_13240', 'H3B104J_12990', 
#        'H3B104X_13200', 'H3B202X_12850', 'H3B203J_13370', 'H3B203M_12480', 
#        'H3B206M_12830', 'H3B206M_12840', 'H3B209X_13340', 'H4B111J_13560', 'H4B111J_13570',
#        'H4B202J_12880', 'H4B204J_13330', 'H4B205J_12990', 'H4B211M_13000',
#        'H4B402J_12600', 'H4B406M_13450', 'H4B406M_13460', 'H4B412M_13240',
#        'H4B501J_12890', 'H4B503X_12670', 'H4B504J_13460', 'H4B505J_12880',
#        'APS55_RS03850', 'LDX55_06325', 'K2W83_RS06180']
# GS2 = ['G0403_13120', 'H1B302M_12900', 'H3B101A_13260', 'H3B104J_13020_2', 
#        'H3B104X_13220', 'H3B202X_12860', 'H3B203J_13390', 'H4B202J_12890_2', 
#        'H4B204J_13350_2', 'H4B504J_13480', 'H4B505J_12900', 'APS55_RS03845', 
#        'LDX55_06335_2', 'K2W83_RS06185']
# BRS = ['A1401_12760', 'FHON2_13550', 'G0403_13110', 'H1B302M_12890', 
#        'H3B101A_13250', 'H3B104J_13020', 'H3B104X_13210', 'H3B203J_13380', 
#        'H3B209X_13350', 'H4B202J_12890', 'H4B204J_13340', 'H4B204J_13350', 
#        'H4B504J_13470', 'H4B505J_12890', 'LDX55_06330', 'LDX55_06335', 
#        'K2W83_RS06185_2']
# S2a = ['A0901_13270', 'A1001_12300', 'A1805_12810', 'H1B104J_13000', 
#        'H3B203M_12470', 'H4B206J_13400', 'H4B402J_12590', 'H4B412M_13230',
#        'H4B503X_12660']
# S2b = ['A1805_12800', 'H1B104J_12990', 'H4B412M_13220']
# S1 = ['G0101_12790', 'H4B205J_12980', 'H4B211M_12990', 'H4B501J_12880', 
#       'K2W83_RS06175']
# S3 = ['A0901_13360', 'A1001_12360', 'A1404_13450', 'H3B203M_12520', 
#       'H4B206J_13490', 'H4B402J_12660', 'H4B412M_13310', 'H4B503X_12760']

# gene_tags = {
#       'GS1':  ['A1001_12310', 'A1003_12540', 'A1202_13520', 'A1401_12750', 'A1805_12820',
#                'FHON2_13540', 'G0101_12800', 'G0403_13100', 'H1B104J_13010', 
#                'H1B105A_12300', 'H1B302M_12880', 'H3B101A_13240', 'H3B104J_12990', 
#                'H3B104X_13200', 'H3B202X_12850', 'H3B203J_13370', 'H3B203M_12480', 
#                'H3B206M_12830', 'H3B206M_12840', 'H3B209X_13340', 'H4B111J_13560', 
#                'H4B111J_13570', 'H4B202J_12880', 'H4B204J_13330', 'H4B205J_12990', 
#                'H4B211M_13000', 'H4B402J_12600', 'H4B406M_13450', 'H4B406M_13460', 
#                'H4B412M_13240', 'H4B501J_12890', 'H4B503X_12670', 'H4B504J_13460', 
#                'H4B505J_12880', 'APS55_RS03850', 'LDX55_06325', 'K2W83_RS06180'], 
#        'GS2': ['G0403_13120', 'H1B302M_12900', 'H3B101A_13260', 'H3B104J_13020_2', 
#                'H3B104X_13220', 'H3B202X_12860', 'H3B203J_13390', 'H4B202J_12890_2', 
#                'H4B204J_13350_2', 'H4B504J_13480', 'H4B505J_12900', 'APS55_RS03845', 
#                'LDX55_06335_2', 'K2W83_RS06185'],
#        'BRS': ['A1401_12760', 'FHON2_13550', 'G0403_13110', 'H1B302M_12890', 
#                'H3B101A_13250', 'H3B104J_13020', 'H3B104X_13210', 'H3B203J_13380', 
#                'H3B209X_13350', 'H4B202J_12890', 'H4B204J_13340', 'H4B204J_13350', 
#                'H4B504J_13470', 'H4B505J_12890', 'LDX55_06330', 'LDX55_06335', 
#                'K2W83_RS06185_2'],
#        'NGB': []}

gene_no = {'GS1':  39, 'GS2': 14, 'BRS': 17, 'NGB': 27}
       # 'S2a': ['A0901_13270', 'A1001_12300', 'A1805_12810', 'H1B104J_13000', 
       #         'H3B203M_12470', 'H4B206J_13400', 'H4B402J_12590', 'H4B412M_13230',
       #         'H4B503X_12660'], 
       # 'S2b': ['A1805_12800', 'H1B104J_12990', 'H4B412M_13220'],
       # 'S1': ['G0101_12790', 'H4B205J_12980', 'H4B211M_12990', 'H4B501J_12880', 
       #        'K2W83_RS06175'], 
       # 'S3': ['A0901_13360', 'A1001_12360', 'A1404_13450', 'H3B203M_12520', 
       #       'H4B206J_13490', 'H4B402J_12660', 'H4B412M_13310', 'H4B503X_12760']}

# matrix = np.random.random((50,50))

# plt.imshow(matrix)
# plt.colorbar()
# plt.show()
    
vmin = 0
vmax = 1.5
replace_str = {'K2W83_RS': 'DSM_', 'APS55_RS': 'MP2_', 'LDX55_': 'IBH001_'}

for gene in genes:
    file = f'results/{gene}/dNdS.tsv'
    with open(file) as comp_file:
        file_info = pd.read_csv(comp_file, sep = '\t')
        file_info = file_info.sort_values(by=['locus1', 'locus2'])
        file_info = file_info.reset_index()
        for key in replace_str.keys():
            file_info[['locus1', 'locus2']] = file_info[['locus1', 'locus2']].apply(lambda col: col.str.replace(key, replace_str[key]))
        dim = gene_no[gene] #len(file_info) #gene_no[gene]
        plot_matrix = np.full((dim, dim), np.nan)
        # Extract unique loci
        unique_loci = np.unique(file_info[['locus1', 'locus2']].values)
    for index, row in file_info.iterrows():
        # if any(key in row['locus1'] for key in replace_str.keys()):
        #     rep_key = row['locus1'].split('0')[0]
        #     file_info.loc[index, 'locus1'] = row['locus1'].replace(rep_key, replace_str[rep_key])
        # if any(key in row['locus2'] for key in replace_str.keys()):
        #     rep_key = row['locus2'].split('0')[0]
        #     file_info.loc[index, 'locus2'] = row['locus2'].replace(rep_key, replace_str[rep_key])
        # for i in range(len(unique_loci)):
        #     for key in replace_str.keys():
        #         unique_loci[i] = unique_loci[i].replace(key, replace_str[key])
        #         print(unique_loci[i])
    
        #plot_matrix[index, index] = 100
        n = np.where(unique_loci == row['locus1'])[0] #[0]
        m = np.where(unique_loci == row['locus2'])[0] #[0]
        plot_matrix[n, m] = row['dS']
        plot_matrix[m, n] = row['dN']
    plot_matrix[plot_matrix > 1.5] = 1.5
    plt.style.use('seaborn-pastel')
    
    #extent = [0, plot_matrix.shape[1], 0, plot_matrix.shape[0]]
    num_rows, num_cols = plot_matrix.shape
    fig_width = num_cols * 0.5
    fig_height = num_rows * 0.5
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    #fig, ax = plt.subplots(figsize=(18, 6))
    plt.imshow(plot_matrix, vmin=vmin, vmax=vmax, cmap='magma')
    row_indices = ['Row 1', 'Row 2', 'Row 3']
    col_indices = ['Col 1', 'Col 2', 'Col 3']
    plt.xticks(np.arange(len(unique_loci)), unique_loci, rotation=90)
    plt.yticks(np.arange(len(unique_loci)), unique_loci)
    colbar = plt.colorbar(ticks=np.arange(0, 1.75, 0.10),  shrink=0.8)
    #colbar.set_label(r'$d_N$ and $d_S$ value')

    #plt.show()
    mpl.rcParams.update(mpl.rcParamsDefault)

    plt.savefig(f'{gene}_dNdS_matrix.svg', format = 'svg')
    plt.show()
        
