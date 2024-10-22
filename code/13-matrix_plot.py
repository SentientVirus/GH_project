#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 16:54:02 2023

This is a plot to create matrices for genes GS1-2, BRS and NGB, where
the dS pairwise values are plotted on the top half and the dN pairwise
values are plotted on the bottom half.

@author: Marina Mota Merlo
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib as mpl

workdir = os.path.expanduser('~') + '/GH_project'
outdir = f'{workdir}/plots/heatmaps'
genes = ['GS1', 'GS2', 'BRS', 'NGB']

gene_no = {'GS1':  33, 'GS2': 14, 'BRS': 17, 'NGB': 27}

    
vmin = 0
vmax = 1.5
replace_str = {'K2W83_RS': 'DSM_', 'APS55_RS': 'MP2_', 'LDX55_': 'IBH001_'}

if not os.path.exists(outdir):
    os.makedirs(outdir)

for gene in genes:
    file = f'{workdir}/results/{gene}/dNdS.tsv'
    with open(file) as comp_file:
        file_info = pd.read_csv(comp_file, sep = '\t')
        file_info = file_info.sort_values(by=['locus1', 'locus2'])
        file_info = file_info.reset_index()
        for key in replace_str.keys():
            file_info[['locus1', 'locus2']] = file_info[['locus1', 'locus2']].apply(lambda col: col.str.replace(key, replace_str[key]))
        dim = gene_no[gene]
        plot_matrix = np.full((dim, dim), np.nan)
        # Extract unique loci
        unique_loci = np.unique(file_info[['locus1', 'locus2']].values)
        
    for index, row in file_info.iterrows():
        n = np.where(unique_loci == row['locus1'])[0] #[0]
        m = np.where(unique_loci == row['locus2'])[0] #[0]
        plot_matrix[n, m] = row['dS']
        plot_matrix[m, n] = row['dN']
    plot_matrix[plot_matrix > 1.5] = 1.5
    plt.style.use('seaborn-pastel')
    
    num_rows, num_cols = plot_matrix.shape
    fig_width = 18
    fig_height = 18
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.imshow(plot_matrix, vmin=vmin, vmax=vmax, cmap='magma')
    row_indices = ['Row 1', 'Row 2', 'Row 3']
    col_indices = ['Col 1', 'Col 2', 'Col 3']
    plt.xticks(np.arange(len(unique_loci)), unique_loci, rotation = 90, fontsize = 16)
    plt.yticks(np.arange(len(unique_loci)), unique_loci, fontsize = 16)
    plt.xlabel('Locus tag 1', size = 24)
    plt.ylabel('Locus tag 2', size = 24)
    colbar = plt.colorbar(ticks = np.arange(0, 1.75, 0.10),  shrink = 0.8)
    colbar.ax.tick_params(labelsize = 16)
    colbar.set_label(r'$d_N$ and $d_S$ value', size = 24)
    colbar.outline.set_visible(False)
    fig.frameon = True

    mpl.rcParams.update(mpl.rcParamsDefault)

    plt.savefig(f'{outdir}/{gene}_dNdS_matrix.svg', format = 'svg')
    plt.savefig(f'{outdir}/{gene}_dNdS_matrix.png', format = 'png')
    plt.savefig(f'{outdir}/{gene}_dNdS_matrix.tiff', format = 'tiff')
    plt.show()
        
