#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:34:46 2024

This script will create violin plots from pairwise dS values.

@author: Marina Mota-Merlo
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import os
import numpy as np
import random

outdir = os.path.expanduser('~') + '/GH_project/plots/boxplot_dS'
outfile = 'dS_boxplot.png'

basal_strains = ['A1001', 'A1404']

inpath = os.path.expanduser('~') + '/GH_project/results'
data_file = 'dNdS.tsv'
gene_types = ['GS1', 'BRS_clade', 'GS2', 'S1', 'S2a', 'S3']
gene_types2 = ['S1', 'S2a', 'S3']
fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (9, 4))


if not os.path.exists(outdir):
    os.makedirs(outdir)


data_dict = {}
all_data = []
for j in range(len(gene_types)):
    if j == 0:
        data_dict[0] = []
    elif j == 3:
        data_dict[1] = []
    gtype = gene_types[j]
    infolder = f'{inpath}/{gtype}'
    df = pd.read_csv(f'{infolder}/{data_file}', sep = '\t')
    df['keep'] = df['locus1'].apply(lambda x: False if any(bs in x for bs in basal_strains) else True)
    df = df[df['keep']==True]
    df['keep'] = df['locus2'].apply(lambda x: False if any(bs in x for bs in basal_strains) else True)
    df = df[df['keep']==True]
    data = list(df['dS'])
    dat = [x if x < 10 else 10 for x in data]
    loc1 = list(df['locus1'])
    loc2 = list(df['locus2'])
    labels = [f'{loc1[i]}-{loc2[i]}' for i in range(len(loc1))]
    if j < 3:
        data_dict[0].append(dat)
    else:
        data_dict[1].append(dat )
    
line_color = 'black' #'#9bfffd'
median_color = 'cyan'
colors = ['#ff7575', '#e875ff', '#ff75b6', '#ffb875', '#ffd775', '#fff575']
gene_types[1] = 'BRS'
fontdict =  {'fontsize': 14} #, 'fontweight': 'bold'}

positions = [1, 2, 3]
widths = [0.3]*3
bplot = axs[0].boxplot(data_dict[0], vert=False, widths = widths, patch_artist = True)
# i = 1
# for lst in data_dict[0]:
#     for el in lst:
#         axs[0].plot(el, i+np.random.uniform(-0.05, 0.05), '.', alpha=0.1, 
#                     color = line_color, zorder = 10, markersize = 5)
#     i += 1

axs[0].yaxis.set_major_locator(ticker.FixedLocator(positions))
axs[0].yaxis.set_major_formatter(ticker.FixedFormatter(gene_types[0:3]))
for element in ['boxes','whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(bplot[element], color=line_color)
plt.setp(bplot['whiskers'], zorder = 5)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
for median in bplot['medians']:
    median.set_color(median_color)#'#00A896')#'#A7E8BD') #'#032B43')
start, end = axs[0].get_xlim()
axs[0].xaxis.set_ticks(range(0, round(end+0.5)))
axs[0].set_title('GH70 $d_{S}$', fontdict = fontdict)


positions = [1, 2, 3]
bplot2 = axs[1].boxplot(data_dict[1], vert=False, widths = widths, patch_artist = True)
# i = 1
# for lst in data_dict[1]:
#     for el in lst:
#         axs[1].plot(el, i+np.random.uniform(-0.15, 0.15), '.', alpha=0.1, 
#                     color = line_color, zorder = 10, markersize = 5) #random.uniform(-0.3, 0.3)
    # i += 1
    
axs[1].yaxis.set_major_locator(ticker.FixedLocator(positions))
axs[1].yaxis.set_major_formatter(ticker.FixedFormatter(gene_types[3:]))
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(bplot2[element], color=line_color)
plt.setp(bplot2['whiskers'], zorder = -1)
for patch, color in zip(bplot2['boxes'], colors[3:]):
    patch.set_facecolor(color)
for median in bplot2['medians']:
    median.set_color(median_color)
start, end = axs[1].get_xlim()
axs[1].xaxis.set_ticks(np.arange(0, round(end+0.25, 0), 0.25))
axs[1].set_title('GH32 $d_{S}$', fontdict = fontdict)
  
axcolor = 'lavender'
for ax in axs:
    ax.grid(color = axcolor, zorder = -2, linestyle = '--')
    ax.yaxis.grid(False)
    ax.set_facecolor('ghostwhite')
    ax.spines['bottom'].set_color(axcolor)
    ax.spines['top'].set_color(axcolor)
    ax.spines['right'].set_color(axcolor)
    ax.spines['left'].set_color(axcolor)
    ax.tick_params(axis='x', colors='black', direction='out', length=3, width=2)
    ax.tick_params(axis='y', colors='black', direction='out', length=3, width=2, labelsize = 12)
        
plt.savefig(f'{outdir}/{outfile}')
plt.savefig(f'{outdir}/{outfile.replace("png", "svg")}')
plt.savefig(f'{outdir}/{outfile.replace("png", "pdf")}')
    