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

prefix = 'all_subsets_positions_aln'
prefix2 = 'all_subsets_7methods_5+_filtered_assort'

infolder = os.path.expanduser('~') + '/GH_project/RDP5_analysis'
infile = f'{infolder}/RDP_output/{prefix2}.csv'
breakfile = f'{infolder}/RDP_output/{prefix2}.csvBreakpointPositions.csv'
outfolder = f'{infolder}/results'
gene_pos = f'{infolder}/files/tab/{prefix}.tab'

df = pd.read_csv(infile)

breakpoints = pd.read_csv(breakfile)
breaks = breakpoints['Breakpoint position']

gene_positions = pd.read_csv(gene_pos, sep = '\t')
genes = gene_positions[gene_positions['strain'] == 'A0901']


#Gaps
gap1 = (3863, 5299)
gap2 = (16804, 18101)
gap3 = (21286, 23073)

gaps = [gap1, gap2, gap3]
# #subset
# GS1 = (6873, 11240)
# BRS = (11510, 16759)
# GS2 = (16943, 26071)

# name_dict = {'gtf2a': 'S2a', 'gtf2b': 'S2b', 'gtf3': 'S3', 'glu1': 'GS1',
#              'glu2': 'GS2', 'brsu': 'BRS'}
    
# GS1_list = [inf, 0]
# BRS_list = [inf, 0]
# GS2_list = [inf, 0]
# for index, row in dom.iterrows():
#     if row['start'] < 9000:
#         GS1_list[0] = min(GS1_list[0], row['start'])
#         GS1_list[1] = max(GS1_list[1], row['end'])
#     elif row['start'] < 15000:
#         BRS_list[0] = min(BRS_list[0], row['start'])
#         BRS_list[1] = max(BRS_list[1], row['end'])
#     else:
#         GS2_list[0] = min(GS2_list[0], row['start'])
#         GS2_list[1] = max(GS2_list[1], row['end'])


x = df['Position in alignment']
y = df[' Recombination breakpoint number (200nt win)']
conf95 = df[' Upper 95% CI']
conf99 = df[' Upper 99% CI']

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

f, (ax2, ax, ax3) = plt.subplots(3, 1, figsize = (26, 8), gridspec_kw={'height_ratios': [1, 8, 2]})
ax.margins(x=0.01)
ax.set_xlim(0, max(x))

ax.plot(x, conf99, color = '#FAD4C0', alpha = 0.6, zorder = 5)
ax.fill_between(x, conf99, 0, color = '#FEE9E1', alpha = 0.5, zorder = 0)

ax.plot(x, conf95, color = '#D58870', zorder = 15)
ax.fill_between(x, conf95, 0, color = '#E0937A', alpha = 0.6, zorder = 10)
ax.set_ylim(bottom=0, top = max(conf99) + 0.5)

ax.plot(x, y, color = 'black', zorder = 20)

ax.axvline(x=15969, color = 'r', linestyle = '--', linewidth = 3, zorder = 25)

# ax.plot(x[x < GS1[0]], conf95[x < GS1[0]], color = 'grey')
# ax.plot(x[(x >= GS1[0]) & (x <= GS1[1])], conf95[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF0000')
# ax.plot(x[(x > GS1[1]) & (x < BRS[0])], conf95[(x > GS1[1]) & (x < BRS[0])], color = 'grey')
# ax.plot(x[(x >= BRS[0]) & (x <= BRS[1])], conf95[(x >= BRS[0]) & (x <= BRS[1])], color = '#D930BD')
# ax.plot(x[(x > BRS[1]) & (x < GS2[0])], conf95[(x > BRS[1]) & (x < GS2[0])], color = 'grey')
# ax.plot(x[(x >= GS2[0]) & (x <= GS2[1])], conf95[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF009B')
# ax.plot(x[x > GS2[1]], conf95[x > GS2[1]], color = 'grey')

# ax.plot(x[x < GS1[0]], conf99[x < GS1[0]], color = 'grey', alpha = 0.3, zorder = -1)
# ax.plot(x[(x >= GS1[0]) & (x <= GS1[1])], conf99[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF0000', alpha = 0.3)
# ax.plot(x[(x > GS1[1]) & (x < BRS[0])], conf99[(x > GS1[1]) & (x < BRS[0])], color = 'grey', alpha = 0.3)
# ax.plot(x[(x >= BRS[0]) & (x <= BRS[1])], conf99[(x >= BRS[0]) & (x <= BRS[1])], color = '#D930BD', alpha = 0.3)
# ax.plot(x[(x > BRS[1]) & (x < GS2[0])], conf99[(x > BRS[1]) & (x < GS2[0])], color = 'grey', alpha = 0.3)
# ax.plot(x[(x >= GS2[0]) & (x <= GS2[1])], conf99[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF009B', alpha = 0.3)
# ax.plot(x[x > GS2[1]], conf99[x > GS2[1]], color = 'grey', alpha = 0.3)

# ax.fill_between(x[x < GS1[0]], conf95[x < GS1[0]], color = '#AFAFAF')
# ax.fill_between(x[(x >= GS1[0]) & (x <= GS1[1])], conf95[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF6060')
# ax.fill_between(x[(x > GS1[1]) & (x < BRS[0])], conf95[(x > GS1[1]) & (x < BRS[0])], color = 'grey')
# ax.fill_between(x[(x >= BRS[0]) & (x <= BRS[1])], conf95[(x >= BRS[0]) & (x <= BRS[1])], color = '#DC6FCA')
# ax.fill_between(x[(x > BRS[1]) & (x < GS2[0])], conf95[(x > BRS[1]) & (x < GS2[0])], color = 'grey')
# ax.fill_between(x[(x >= GS2[0]) & (x <= GS2[1])], conf95[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF60C1')
# ax.fill_between(x[x > GS2[1]], conf95[x > GS2[1]], color = '#AFAFAF')

# ax.fill_between(x[x < GS1[0]], conf99[x < GS1[0]], color = '#AFAFAF', alpha = 0.3)
# ax.fill_between(x[(x >= GS1[0]) & (x <= GS1[1])], conf99[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF6060', alpha = 0.3)
# ax.fill_between(x[(x > GS1[1]) & (x < BRS[0])], conf99[(x > GS1[1]) & (x < BRS[0])], color = 'grey', alpha = 0.3)
# ax.fill_between(x[(x >= BRS[0]) & (x <= BRS[1])], conf99[(x >= BRS[0]) & (x <= BRS[1])], color = '#DC6FCA', alpha = 0.3)
# ax.fill_between(x[(x > BRS[1]) & (x < GS2[0])], conf99[(x > BRS[1]) & (x < GS2[0])], color = 'grey', alpha = 0.3)
# ax.fill_between(x[(x >= GS2[0]) & (x <= GS2[1])], conf99[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF60C1', alpha = 0.3)
# ax.fill_between(x[x > GS2[1]], conf99[x > GS2[1]], color = '#AFAFAF', alpha = 0.3)

# ax.fill_between(x[(x >= GS1_list[0]) & (x <= GS1_list[1])], conf95[(x >= GS1_list[0]) & (x <= GS1_list[1])], color = '#FF0000')
# ax.fill_between(x[(x >= BRS_list[0]) & (x <= BRS_list[1])], conf95[(x >= BRS_list[0]) & (x <= BRS_list[1])], color = '#D930BD')
# ax.fill_between(x[(x >= GS2_list[0]) & (x <= GS2_list[1])], conf95[(x >= GS2_list[0]) & (x <= GS2_list[1])], color = '#FF009B')

# stemobj = ax.stem(x, y, markerfmt = ' ')
# plt.setp(stemobj, 'color', 'black')
# plt.setp(stemobj, 'zorder', 30)
ax2.set_xlim(0, max(x))
ax.set_xticks(range(1, max(x), 2000))
ax2.margins(y = 0)
ax2.set_xticks(range(1, max(x), 2000))
ax2.plot(breaks, [1]*len(breaks), '|', color = 'black')
ax2.set_axis_off()
f.subplots_adjust(wspace=0, hspace=0)

# ax2.axvline(x=15969, color='r', linestyle='--', linewidth=2)

right_side = ax.spines['right']
right_side.set_visible(False)
top_side = ax.spines['top']
top_side.set_visible(False)

ax3.set_xlim(0, max(x))
ax3.margins(x = 0.01)
ax3.plot(x, x/1000, alpha = 0)
ax3.set_xticks(range(min(x), max(x), 2000))
ax3.get_yaxis().set_visible(False)

ax3.axvline(x=15969, color='r', linestyle='--', linewidth=3)

colors = ['#9DAAB7', '#6586A3', '#647D88', '#626D75', '#B1BFD4']


for index, gene in genes.iterrows():
    gene_zorder = 2
    basecolor = colors[index%5] #'#AFAFAF'
    basewidth = 1000
    headwidth = 1000
    linewidth = 5
    alpha = 1
    
    basewidth = 500
    headwidth = 500
    if gene['end'] - gene['start'] < 800: # or gene['gene_name'] == 'unk.':
        headwidth = 0 #(gene['end'] - gene['start'])/1.6
        linewidth = 0
        alpha = 0.8
        gene_zorder = 1
        
    if gene['strand'] == 1:
        ax3.arrow(gene['end'], 0, dx = gene['start'] - gene['end'], dy = 0, 
              facecolor = basecolor, length_includes_head = True, 
              width = basewidth, shape = 'full', head_width = headwidth, 
              edgecolor = 'black', linewidth = linewidth, alpha = alpha,
              zorder = gene_zorder)
    else:
        ax3.arrow(gene['start'], 0, dx = gene['end'] - gene['start'], dy = 0, 
              facecolor = basecolor, length_includes_head = True, 
              width = basewidth, shape = 'full', head_width = headwidth, 
              edgecolor = 'black', linewidth = linewidth, alpha = alpha,
              zorder = gene_zorder)
        
# brs_start =  min(genes[genes['name'] == 'brsu']['start'])
# brs_end = max(genes[genes['name'] == 'brsu']['end'])
# ax3.arrow(brs_end, 0, dx = brs_start - brs_end, dy = 0, 
#       facecolor = '#DC6FCA', length_includes_head = True, 
#       width = 1000, shape = 'full', head_width = 1000, 
#       edgecolor = 'black', linewidth = 5, alpha = 1)  
# leg = ax3.text(brs_end-(brs_end-brs_start)/(1+2/3), -60, 
#                'BRS', fontsize = 18)

for gap in gaps:
    gene_zorder = 0
    basecolor = '#F2F2F2'
    basewidth = 500
    headwidth = 0
    linewidth = 0
    alpha = 1
    ax3.arrow(gap[0], 0, dx = gap[1] - gap[0], dy = 0, 
          facecolor = basecolor, length_includes_head = False, 
          width = basewidth, shape = 'full', head_width = headwidth, 
          edgecolor = 'black', linewidth = 0, alpha = alpha,
          zorder = gene_zorder)
    
    ax.arrow(gap[0], 0, dx = gap[1] - gap[0], dy = 0, 
          facecolor = basecolor, length_includes_head = False, 
          width = basewidth, shape = 'full', head_width = headwidth, 
          edgecolor = 'black', linewidth = linewidth, alpha = 0.75,
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