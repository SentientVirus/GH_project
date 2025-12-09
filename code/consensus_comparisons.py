#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 15:49:44 2025

Code to plot the measures of similarity from Jalview.

@author: Marina Mota-Merlo
"""

#OBS! Change back to the consensus instead of conservation!

import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator

def format_x(x, pos):
    '''This function specifies how to display the tick labels of the x axis'''
    return f'{x+1:.0f}'

# colors = ['#35605A', '#74A7C1', '#D8E4FF', '#8E6F3D', '#563440']
# colors = ['#B30021', '#CF5B3E', '#CDD0D6', '#ED9F45', '#FFEA8C']
colors = ['#B30021', '#CF5B3E', '#CDD0D6', '#ED9F45', '#FFEA8C']
colors2 = ['#ED9F45', 'black']
cmap1 = LinearSegmentedColormap.from_list('red', list(reversed(colors)))
cmap2 = LinearSegmentedColormap.from_list('cyan', list(colors2))

motif_pos = {1: (693, 700), 2: (187, 197), 3: (223, 235), 4: (310, 323),
             5: (144, 157), 6: (625, 637), 7: (380, 389), 8: (637, 661)}

datasets = ['GS4', 'GS3', 'GS2', 'GS1', 'GS', 'GS_and_BrS']

motif_list = [(motif_pos[i][0], motif_pos[i][1] - motif_pos[i][0]) for i in motif_pos.keys()]

indir = os.path.expanduser('~') + '/GH_project/conservation_consensus'
fontsize_title = 64

# conserv = [f'{indir}/{dataset}_conservation.csv' for dataset in datasets]
consens = [f'{indir}/{dataset}_consensus.csv' for dataset in datasets]
outplot = f'{indir}/heatmap.png'
out_tsv = f'{indir}/conserved_positions.tsv'

class conserved:
    def __init__(self, pos, GS1, GS2, GS3, GS4):
        self.pos = pos
        self.GS1 = GS1
        self.GS2 = GS2
        self.GS3 = GS3
        self.GS4 = GS4
        self.GS_list = [GS1, GS2, GS3, GS4]
    def update_GS(self):
        self.GS_list = [self.GS1, self.GS2, self.GS3, self.GS4]
    def __str__(self):
        to_print = f'Position {self.pos}:\n'
        for i in range(len(self.GS_list)):
            if self.GS_list[i] != '':
                to_print += f'\t- GS{i+1}, {self.GS_list[i]}\n'
        return to_print

# def visualize(cell): #just so the plots show up
#     ax = plt.subplot(cell)
#     ax.plot()
#     ax.get_xaxis().set_visible(False)
#     ax.get_yaxis().set_visible(False)

fig = plt.figure(constrained_layout=False, figsize=(80, 50))

pair_spec = gridspec.GridSpec(2, 1, height_ratios = [1, 6]) 
# plt.subplots_adjust(hspace=0.0)
spec = gridspec.GridSpec(nrows=6, ncols=1, figure=fig, hspace = 0.8)
                         # height_ratios = [1, 4, 1, 4, 1, 4, 1, 4, 1, 4]) #, width_ratios=length_list)

# for i in range(5):
#     visualize(gs1[0])
#     visualize(gs2[0])

# plt.show()

dataset_dict = {}
seq_dict = {}

axs = []

for j in range(len(consens)):
    # k = j*2+1
    name = os.path.basename(consens[j]).replace('_consensus.csv', '').replace('_', ' ')
    with open(consens[j]) as csvfile: #, open(conserv[j]) as csvcomp:
        df = pd.read_csv(csvfile, sep = ',', header = None)
        # df2 = pd.read_csv(csvcomp, sep = ',', header = None)
        to_plot = df.iloc[0]
        gaps = df.iloc[1]
        to_plot_list = list(to_plot)[1:-1]
        gap_list = list(gaps)[1:-1]
        sequence = [pos.strip().split(' ')[0] for pos in gap_list]
        seq_dict[name] = sequence
        
        to_plot_float = [float(to_plot_list[i].strip()) if float(to_plot_list[i].strip()) == 100.00 else np.nan if gap_list[i].strip() != '' else 0.00 for i in range(len(to_plot_list))]

        # to_plot_float = [to_plot_list[i] if gap_list[i].strip() != '' else np.nan for i in range(len(to_plot_list))]
        dataset_dict[name] = to_plot_float
        to_plot_float = pd.DataFrame(to_plot_float, columns=['Consensus'])

        
        gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = spec[j], hspace = 0.1, height_ratios = [1, 6])
        
        bar_ax = fig.add_subplot(gs[0])
        ax = fig.add_subplot(gs[1])
        ax.set_facecolor('#FEFDED')
        
        # # bar_ax = fig.add_subplot(spec[k-1, 0])
        bar_ax.set_xlim(1, len(to_plot_float)+1)
        bar_ax.xaxis.set_major_locator(MultipleLocator(50)) #Set ticks every 50 positions
        bar_ax.xaxis.set_major_formatter(format_x) #Apply function to display the index without decimals
        bar_ax.broken_barh(motif_list, (0, 0.1), color=['black']*7 + ['slategrey'])
        # bar_ax.set_yticklabels([])
        bar_ax.get_yaxis().set_visible(False)
        bar_ax.get_xaxis().set_visible(False)
        bar_ax.set_title(name, fontsize = fontsize_title, #style = 'italic',
                      fontname = 'Arial') # Add gene names at the top of the plot
        # [bar_ax.spines[element].set_visible(False) for element in ['top', 'bottom', 'left', 'right']]
        
        # # ax = fig.add_subplot(spec[k, 0]) # Set the right row and column
        # ax.set_title(name, fontsize = fontsize_title, #style = 'italic',
        #               fontname = 'Arial') # Add gene names at the top of the plo
        htmp = sns.heatmap(to_plot_float.T, cmap = cmap1, center = 5.5, ax = ax,
                          vmin = 0, vmax = 100, linewidths = 0,
                          cbar = False, zorder = 0) # Plot
        ax.patch.set(lw = 1, ec = 'black') # Border of the subplots
        ax.xaxis.set_major_locator(MultipleLocator(50)) # Set ticks every 50 positions
        ax.xaxis.set_major_formatter(format_x) # Apply function to display the index without decimals

        plt.tick_params(axis='x', which='major', labelsize = 50,
                        length = 25, width = 5, rotation = 0) # Increase font size of ax ticks
        plt.tick_params(axis='y', which='major', labelsize = 60,
                        length = 0) # Increase font size of ax ticks
        axs.append(ax)
    
keys = list(dataset_dict.keys())[:-2]
comparison= {}
for k in range(len(keys)):
    g1 = keys[k]
    comp = [1 for n in range(len(dataset_dict[g1]))]
    for m in range(len(keys)):
        if k != m:
            g2 = keys[m]
            comp = [np.nan if ((dataset_dict[g1][n] != 100) or (dataset_dict[g2][n] != 100) or (seq_dict[g1][n] == seq_dict[g2][n])) else comp[n] for n in range(len(dataset_dict[g1]))] 
    comparison[g1] = comp
    comp_df = pd.DataFrame(comp, columns=['Unique'])
    new_ax = axs[k].twinx()
    htmp = sns.heatmap(comp_df.T, cmap = cmap2, center = 5.5, ax = new_ax,
                      vmin = 0, vmax = 1, linewidths = 0,
                      cbar = False, zorder = 0) #, cbar_kws = {'fraction': 0.05, 'aspect': 10, 'ticks': [0, 5, 10]}) # Plot
    axs[k].xaxis.set_major_locator(MultipleLocator(50)) #Set ticks every 50 positions
    axs[k].xaxis.set_major_formatter(format_x) #Apply function to display the index without decimals
    new_ax.set_yticklabels([])
    new_ax.get_yaxis().set_visible(False)
    
plt.savefig(outplot, bbox_inches='tight')

conserved_all = [n for n in range(len(comparison['GS1'])) if comparison['GS1'][n] == True == comparison['GS2'][n] == comparison['GS3'][n] == comparison['GS4'][n]]

positions = []
with open(out_tsv, 'w') as tsvfile:
    tsvfile.write('Alignment_position\tGS1_aa\tGS2_aa\tGS3_aa\tGS4_aa\n')
    for n in range(len(comparison['GS1'])):
        if any([comparison[x][n] == True for x in comparison.keys()]):
                cons = conserved(n+1, seq_dict['GS1'][n], seq_dict['GS2'][n], seq_dict['GS3'][n], seq_dict['GS4'][n])
                print(cons)
                tsvfile.write(f'{cons.pos}\t{cons.GS1}\t{cons.GS2}\t{cons.GS3}\t{cons.GS4}\n')