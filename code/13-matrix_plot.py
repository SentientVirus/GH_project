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
from matplotlib.transforms import ScaledTranslation

workdir = os.path.expanduser('~') + '/GH_project'
outdir = f'{workdir}/plots/heatmaps'
genes = ['GS1', 'GS2', 'BRS', 'NGB']
outplot = f'{outdir}/GH70_dNdS_matrix.svg'

gene_no = {'GS1':  33, 'GS2': 14, 'BRS': 17, 'NGB': 27}

    
vmin = 0
vmax = 1.5
replace_str = {'K2W83_RS': 'DSM_', 'APS55_RS': 'MP2_', 'LDX55_': 'IBH001_'}

if not os.path.exists(outdir):
    os.makedirs(outdir)

fig_width = 34
fig_height = 30

fig, ax = plt.subplot_mosaic([['A', 'B', ''], ['C', 'D', '']], 
                             figsize=(fig_width, fig_height),
                             gridspec_kw={'width_ratios': [14, 14, 1]})

plt.style.use('seaborn-pastel')
plt.subplots_adjust(hspace = 0.4, wspace = 0.1)

cmap = mpl.cm.magma
norm = mpl.colors.Normalize(vmin = 0, vmax = 1.5)

colbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap),
             cax = ax[''], orientation = 'vertical', 
             ticks = np.arange(0, 1.60, 0.1))
ax[''].yaxis.set_label_position('left')

colbar.ax.tick_params(labelsize = 24)
colbar.set_label(r'$d_N$ and $d_S$ value', size = 24)
colbar.outline.set_linewidth(4)

for gene in genes:
    if gene == 'GS1':
        axn = 'A'
    elif gene == 'GS2':
        axn = 'B'
    elif gene == 'BRS':
        axn = 'C'
    else: axn = 'D'
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
    
    num_rows, num_cols = plot_matrix.shape
    for axis in ['top','bottom','left','right']:
        ax[axn].spines[axis].set_linewidth(4)

    ax[axn].imshow(plot_matrix, vmin=vmin, vmax=vmax, cmap='magma')

    ax[axn].set_xticks(np.arange(len(unique_loci)))
    ax[axn].set_xticklabels(unique_loci, rotation = 90, fontsize = 16)
    ax[axn].set_yticks(np.arange(len(unique_loci)))
    ax[axn].set_yticklabels(unique_loci, fontsize = 16)
    ax[axn].set_xlabel('Locus tag 1', size = 24)
    ax[axn].set_ylabel('Locus tag 2', size = 24)
    ax[axn].set_title(gene, size = 30)
    ax[axn].text(0, 1.1, axn, transform=(
           ax[axn].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
       fontsize = 36, va = 'top', ha = 'left', bbox = dict(facecolor='none', 
                                                           edgecolor='black', 
                                                           pad = 20))


    mpl.rcParams.update(mpl.rcParamsDefault)

plt.savefig(outplot, format = 'svg')
plt.savefig(outplot.replace('svg', 'png'), format = 'png')
plt.savefig(outplot.replace('svg', 'tiff'), format = 'tiff')
plt.show()
        
