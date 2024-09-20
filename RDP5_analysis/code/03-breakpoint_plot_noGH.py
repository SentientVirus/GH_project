#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:46:39 2022

@author: marina
"""

import matplotlib.pyplot as plt
import pandas as pd
from statistics import mean
from math import inf
import os
import numpy as np
from scipy.signal import find_peaks

# =============================================================================
# 0. Define inputs, outputs and paths
# =============================================================================
prefix = 'all_subsets_positions_aln'
prefix2 = 'all_subsets_7methods_5+_filtered_assort'

infolder = os.path.expanduser('~') + '/GH_project/RDP5_analysis'
infile = f'{infolder}/RDP_output/{prefix2}.csv'
breakfile = f'{infolder}/RDP_output/{prefix2}.csvBreakpointPositions.csv'
outfolder = f'{infolder}/results'
gene_pos = f'{infolder}/files/tab/{prefix}.tab'
outfile = f'{infolder}/plots/{prefix2}.png'

if not os.path.exists(os.path.dirname(outfile)):
    os.makedirs(os.path.dirname(outfile))

# =============================================================================
# 1. Read RDP5 output
# =============================================================================

df = pd.read_csv(infile)

breakpoints = pd.read_csv(breakfile)
breaks = breakpoints['Breakpoint position']

gene_positions = pd.read_csv(gene_pos, sep = '\t')
genes = gene_positions[gene_positions['strain'] == 'A0901']

# =============================================================================
# 2. Set gap positions and CDS colors
# =============================================================================
gap1 = (3863, 5299)
gap2 = (16804, 18101)
gap3 = (21286, 23073)

gaps = [gap1, gap2, gap3]

colors = ['#9DAAB7', '#6586A3', '#647D88', '#626D75', '#B1BFD4']

# =============================================================================
# 4. Create variables to be plotted
# =============================================================================
x = df['Position in alignment']
y = df[' Recombination breakpoint number (200nt win)']
conf95 = df[' Upper 95% CI']
conf99 = df[' Upper 99% CI']

# =============================================================================
# Old section to use average values
# =============================================================================
# num = 10

# x200 = [x[i] for i in range(len(x)-1) if (i % num == 0) or i == x[len(x)-1]]
# y200 = [mean(y[int(num*i):int(num*(i+1))]) for i in range(int(len(y)/num))]
# y200.append(mean(y[len(y)-num:]))

# conf95_200 = [mean(conf95[int(num*i):int(num*(i+1))]) for i in range(int(len(conf95)/num))]
# conf95_200.append(mean(conf95[len(conf95)-num:]))
# conf99_200 = [mean(conf99[int(num*i):int(num*(i+1))]) for i in range(int(len(conf99)/num))]
# conf99_200.append(mean(conf99[len(conf99)-num:]))

# y = y.tolist()
# conf95 = conf95.tolist()
# conf99 = conf99.tolist()

# for i in range(1,len(y200)):
#     for j in range(len(y)):
#         if (i-1)*num <= j < i*num:
#             y[j] = y200[i]
#             conf95[j] = conf95_200[i]
#             conf99[j] = conf99_200[i]

# y = pd.Series(y)
# conf95 = pd.Series(conf95)
# conf99 = pd.Series(conf99)

# x200 = pd.Series(x200)
# y200 = pd.Series(y200)

# =============================================================================
# 5. Create a plot with 3 subplots
# =============================================================================
f, (ax2, ax, ax3) = plt.subplots(3, 1, figsize = (26, 8), 
                                 gridspec_kw={'height_ratios': [1, 8.5, 1.5]})

f.subplots_adjust(wspace=0, hspace=0) # Remove space between subplots´

# =============================================================================
# 6. Plot confidence intervals and putative breakpoints
# =============================================================================
ax.margins(x=0.01) # Decrease ax margins
ax.set_xlim(0, max(x)) # Set x axis limits

ax.plot(x, conf99, color = '#FAD4C0', alpha = 0.6, zorder = 5) # Plot 99% CI
ax.fill_between(x, conf99, 0, color = '#FEE9E1', alpha = 0.5, zorder = 0) # Fill the CI

ax.plot(x, conf95, color = '#D58870', zorder = 15) # Plot the 95% CI
ax.fill_between(x, conf95, 0, color = '#E0937A', alpha = 0.6, zorder = 10) # Fill the CI
ax.set_ylim(bottom = 0, top = max(conf99) + 0.5) # Set y axis limits

ax.plot(x, y, color = 'black', zorder = 25) # Plot breakpoint curve

# triangle_dict = {}
y_maxima = find_peaks(y, distance = 50)[0]
for yval in y_maxima:
    if y[yval] > conf99[yval]:
        ax.plot(yval*2, y[yval] + 0.4, marker=(3, 0, 180), markersize=20, 
                mew = 2, color = '#D7FA05', markeredgecolor = 'black',
                zorder = 55)
        # triangle_dict[(yval, y[yval])] = 'red'
    elif y[yval] > conf95[yval]:
        ax.plot(yval*2, y[yval] + 0.4, marker=(3, 0, 180), markersize=20, 
                mew = 2, color = '#F4FEB8', markeredgecolor = 'grey', 
                zorder = 55)
        # triangle_dict[(yval, y[yval])] = 'yellow'
        
        

ax.axvline(x=15969, color = 'r', linestyle = '--', linewidth = 3, zorder = 25) # Add horizontal line between segments

ax.tick_params(          # Remove bottom ticks from x axis
    axis = 'x',          # Axis to be modified
    which = 'both',      # Ticks to be modified (both = major + minor)
    bottom = False,      # Remove ticks along the bottom edge
    labelbottom = False) # Remove labels along the bottom edge

ax.spines['left'].set_zorder(100) # Move the ax main spines to the front
right_side = ax.spines['right']
right_side.set_visible(False) # Remove right spine
top_side = ax.spines['top']
top_side.set_visible(False) # Remove top spine

# =============================================================================
# 7. Plot breakpoint positions
# =============================================================================
ax2.set_xlim(0, max(x)) # Set x axis limits
ax.set_xticks(range(1, max(x), 2000)) # Set x axis ticks
ax2.margins(y = 0) # Set y axis margins
ax2.set_xticks(range(1, max(x), 2000)) # Ser x axis ticks
ax2.plot(breaks, [1]*len(breaks), '|', color = 'black') # Plot breakpoints
ax2.set_axis_off() # Remove axis display
ax2.set_zorder(100) # Move ax forward
ax2.set_yticks([]) # Remove ticks from y axis

# =============================================================================
# 8. Plot genes in region
# =============================================================================
ax3.set_xlim(0, max(x)) # Set ax scale
ax3.set_ylim(-200, 200)
ax3.margins(x = 0.01, y = 0) # Set ax margins
ax3.set_xticks(range(min(x), max(x), 2000)) # Set ticks
ax3.get_yaxis().set_visible(False) # Remove y axis

#Add horizontal line between genomic segments
ax3.axvline(x=15969, color='r', linestyle='--', linewidth=3, zorder = 20)

# Plot CDS
for index, gene in genes.iterrows(): # Loop through CFD
    gene_zorder = 2 # Plot them at the back
    basecolor = colors[index%5] # Get a color from a list of alternating colors
    basewidth = 300 # CDS arrow width
    headwidth = 300 # Width of the arrow head
    linewidth = 2 # Edge width
    alpha = 1 # Transparency (none)
    lc = 'black' # Edge color
    
    if gene['end'] - gene['start'] < 400:
        headwidth = 0
        linewidth = 0
        alpha = 0.8
        gene_zorder = 1
        
    if gene['strand'] == 1:
        ax3.arrow(gene['end'], 0, dx = gene['start'] - gene['end'], dy = 0, 
              facecolor = basecolor, length_includes_head = True, 
              width = basewidth, shape = 'full', head_width = headwidth, 
              edgecolor = lc, linewidth = linewidth,
              alpha = alpha, zorder = gene_zorder)
    else:
        ax3.arrow(gene['start'], 0, dx = gene['end'] - gene['start'], dy = 0, 
              facecolor = basecolor, length_includes_head = True, 
              width = basewidth, shape = 'full', head_width = headwidth, 
              edgecolor = lc, linewidth = linewidth,
              alpha = alpha, zorder = gene_zorder)

for gap in gaps:
    gene_zorder = 0
    basecolor = '#F2F2F2'
    basewidth = 600
    headwidth = 0
    linewidth = 0
    alpha = 0.75
    ax3.arrow(gap[0], 0, dx = gap[1] - gap[0], dy = 0, 
          facecolor = basecolor, length_includes_head = False, 
          width = basewidth, shape = 'full', head_width = headwidth, 
          edgecolor = 'black', linewidth = 0, alpha = alpha,
          zorder = gene_zorder)
    
    ax.arrow(gap[0], 0, dx = gap[1] - gap[0], dy = 0, 
          facecolor = basecolor, length_includes_head = False, 
          width = basewidth, shape = 'full', head_width = headwidth, 
          edgecolor = 'black', linewidth = linewidth, alpha = alpha,
          zorder = 35)
    
    ax2.arrow(gap[0], 0, dx = gap[1] - gap[0], dy = 0, 
          facecolor = basecolor, length_includes_head = False, 
          width = basewidth, shape = 'full', head_width = headwidth, 
          edgecolor = 'black', linewidth = linewidth, alpha = alpha,
          zorder = gene_zorder)


right_side3 = ax3.spines['right']
right_side3.set_visible(False)
left_side3 = ax3.spines['left']
left_side3.set_visible(False)

ax.set_ylabel('Breakpoints per 200nt window', fontsize = 16)
plt.xlabel('Position in the alignment', fontsize = 16)

plt.savefig(outfile)
plt.savefig(outfile.replace('png', 'svg'))
plt.savefig(outfile.replace('png', 'pdf'))