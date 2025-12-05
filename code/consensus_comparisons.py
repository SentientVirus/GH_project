#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 15:49:44 2025

Code to plot the measures of similarity from Jalview.

@author: Marina Mota-Merlo
"""

import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

colors = ['#35605A', '#74A7C1', '#D8E4FF', '#8E6F3D', '#563440']
cmap1 = LinearSegmentedColormap.from_list('christmas', list(reversed(colors)))

indir = os.path.expanduser('~') + '/GH_project/conservation_consensus'
fontsize_title = 16

conserv = sorted([f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('conservation.csv')])
consens = sorted([f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('consensus.csv')])
faa = sorted([f'{indir}/{file}' for file in os.listdir(indir) if file.endswith('.faa')])
outplot = f'{indir}/heatmap.png'

fig = plt.figure(constrained_layout=False, figsize=(10, 10))
spec = gridspec.GridSpec(nrows=5, ncols=1, figure=fig, hspace = 0.8) #, width_ratios=length_list)

for j in range(len(consens)):
    name = os.path.basename(consens[j]).replace('_consensus.csv', '').replace('_', ' ')
    with open(consens[j]) as csvfile, open(conserv[j]) as csvcomp:
        df = pd.read_csv(csvfile, sep = ',', header = None)
        df2 = pd.read_csv(csvcomp, sep = ',', header = None)
        to_plot = df2.iloc[0]
        gaps = df.iloc[1]
        to_plot_list = list(to_plot)[1:-1]
        gap_list = list(gaps)[1:-1]
        
        to_plot_float = [to_plot_list[i] if gap_list[i].strip() != '' else np.nan for i in range(len(to_plot_list))]
        to_plot_float = pd.DataFrame(to_plot_float, columns=['Conservation'])
        
        ax = fig.add_subplot(spec[j, 0]) # Set the right row and column
        ax.set_title(name, fontsize = fontsize_title, #style = 'italic',
                      fontname = 'Arial') # Add gene names at the top of the plo
        h3 = sns.heatmap(to_plot_float.T, cmap = cmap1, center = 5.5, cbar = False, ax = ax) # Plot
        ax.patch.set(lw = 1, ec = 'black') # Border of the subplots
        
plt.savefig(outplot, bbox_inches='tight')
    