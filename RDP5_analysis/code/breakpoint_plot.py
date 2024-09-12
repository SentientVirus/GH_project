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

GH32 = 'other_GS1_S2-3'
GH70 = 'GH70_subset' #'all_BRS'

infolder = 'RDP5_results'
infile = f'{infolder}/{GH70}_breakplot.csv'
outfolder = 'RDP5_results'
GH_pos = f'GH_positions/{GH70}_GH_positions_all.tab'
GH_domains = f'GH_positions/{GH70}_GH_domains.tab'

df = pd.read_csv(infile)

#all_BRS
# GS1 = (6873, 11225)
# BRS = (11495, 16594)
# GS2 = (16776, 23642)

#subset
GS1 = (6873, 11240)
BRS = (11510, 16759)
GS2 = (16943, 26071)

name_dict = {'gtf2a': 'S2a', 'gtf2b': 'S2b', 'gtf3': 'S3', 'glu1': 'GS1',
             'glu2': 'GS2', 'brsu': 'BRS'}

with open(GH_pos) as positions:
    genes = pd.read_csv(positions, delimiter = '\t')
    
with open(GH_domains) as domains:
    dom = pd.read_csv(domains, delimiter = '\t')
    
GS1_list = [inf, 0]
BRS_list = [inf, 0]
GS2_list = [inf, 0]
for index, row in dom.iterrows():
    if row['start'] < 9000:
        GS1_list[0] = min(GS1_list[0], row['start'])
        GS1_list[1] = max(GS1_list[1], row['end'])
    elif row['start'] < 15000:
        BRS_list[0] = min(BRS_list[0], row['start'])
        BRS_list[1] = max(BRS_list[1], row['end'])
    else:
        GS2_list[0] = min(GS2_list[0], row['start'])
        GS2_list[1] = max(GS2_list[1], row['end'])

breakpoints = pd.read_csv(f'{infile}BreakpointPositions.csv')
breaks = breakpoints['Breakpoint position']

x = df['Position in alignment']
y = df[' Recombination breakpoint number (200nt win)']
conf95 = df[' Upper 95% CI']
conf99 = df[' Upper 99% CI']

num = 200

x200 = [x[i] for i in range(len(x)-1) if (i % num == 0) or i == x[len(x)-1]]
y200 = [mean(y[int(num*i):int(num*(i+1))]) for i in range(int(len(y)/num))]
y200.append(mean(y[len(y)-num:]))

conf95_200 = [mean(conf95[int(num*i):int(num*(i+1))]) for i in range(int(len(conf95)/num))]
conf95_200.append(mean(conf95[len(conf95)-num:]))
conf99_200 = [mean(conf99[int(num*i):int(num*(i+1))]) for i in range(int(len(conf99)/num))]
conf99_200.append(mean(conf99[len(conf99)-num:]))

y = y.tolist()
conf95 = conf95.tolist()
conf99 = conf99.tolist()

for i in range(1,len(y200)):
    for j in range(len(y)):
        if (i-1)*num <= j < i*num:
            y[j] = y200[i]
            conf95[j] = conf95_200[i]
            conf99[j] = conf99_200[i]

y = pd.Series(y)
conf95 = pd.Series(conf95)
conf99 = pd.Series(conf99)

x200 = pd.Series(x200)
y200 = pd.Series(y200)

f, (ax2, ax, ax3) = plt.subplots(3, 1, figsize = (18, 8), gridspec_kw={'height_ratios': [1, 8, 2]})
ax.margins(x=0.01)

ax.plot(x[x < GS1[0]], conf95[x < GS1[0]], color = 'grey')
ax.plot(x[(x >= GS1[0]) & (x <= GS1[1])], conf95[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF0000')
ax.plot(x[(x > GS1[1]) & (x < BRS[0])], conf95[(x > GS1[1]) & (x < BRS[0])], color = 'grey')
ax.plot(x[(x >= BRS[0]) & (x <= BRS[1])], conf95[(x >= BRS[0]) & (x <= BRS[1])], color = '#D930BD')
ax.plot(x[(x > BRS[1]) & (x < GS2[0])], conf95[(x > BRS[1]) & (x < GS2[0])], color = 'grey')
ax.plot(x[(x >= GS2[0]) & (x <= GS2[1])], conf95[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF009B')
ax.plot(x[x > GS2[1]], conf95[x > GS2[1]], color = 'grey')

ax.plot(x[x < GS1[0]], conf99[x < GS1[0]], color = 'grey', alpha = 0.3, zorder = -1)
ax.plot(x[(x >= GS1[0]) & (x <= GS1[1])], conf99[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF0000', alpha = 0.3)
ax.plot(x[(x > GS1[1]) & (x < BRS[0])], conf99[(x > GS1[1]) & (x < BRS[0])], color = 'grey', alpha = 0.3)
ax.plot(x[(x >= BRS[0]) & (x <= BRS[1])], conf99[(x >= BRS[0]) & (x <= BRS[1])], color = '#D930BD', alpha = 0.3)
ax.plot(x[(x > BRS[1]) & (x < GS2[0])], conf99[(x > BRS[1]) & (x < GS2[0])], color = 'grey', alpha = 0.3)
ax.plot(x[(x >= GS2[0]) & (x <= GS2[1])], conf99[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF009B', alpha = 0.3)
ax.plot(x[x > GS2[1]], conf99[x > GS2[1]], color = 'grey', alpha = 0.3)

ax.fill_between(x[x < GS1[0]], conf95[x < GS1[0]], color = '#AFAFAF')
ax.fill_between(x[(x >= GS1[0]) & (x <= GS1[1])], conf95[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF6060')
ax.fill_between(x[(x > GS1[1]) & (x < BRS[0])], conf95[(x > GS1[1]) & (x < BRS[0])], color = 'grey')
ax.fill_between(x[(x >= BRS[0]) & (x <= BRS[1])], conf95[(x >= BRS[0]) & (x <= BRS[1])], color = '#DC6FCA')
ax.fill_between(x[(x > BRS[1]) & (x < GS2[0])], conf95[(x > BRS[1]) & (x < GS2[0])], color = 'grey')
ax.fill_between(x[(x >= GS2[0]) & (x <= GS2[1])], conf95[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF60C1')
ax.fill_between(x[x > GS2[1]], conf95[x > GS2[1]], color = '#AFAFAF')

ax.fill_between(x[x < GS1[0]], conf99[x < GS1[0]], color = '#AFAFAF', alpha = 0.3)
ax.fill_between(x[(x >= GS1[0]) & (x <= GS1[1])], conf99[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF6060', alpha = 0.3)
ax.fill_between(x[(x > GS1[1]) & (x < BRS[0])], conf99[(x > GS1[1]) & (x < BRS[0])], color = 'grey', alpha = 0.3)
ax.fill_between(x[(x >= BRS[0]) & (x <= BRS[1])], conf99[(x >= BRS[0]) & (x <= BRS[1])], color = '#DC6FCA', alpha = 0.3)
ax.fill_between(x[(x > BRS[1]) & (x < GS2[0])], conf99[(x > BRS[1]) & (x < GS2[0])], color = 'grey', alpha = 0.3)
ax.fill_between(x[(x >= GS2[0]) & (x <= GS2[1])], conf99[(x >= GS2[0]) & (x <= GS2[1])], color = '#FF60C1', alpha = 0.3)
ax.fill_between(x[x > GS2[1]], conf99[x > GS2[1]], color = '#AFAFAF', alpha = 0.3)

ax.fill_between(x[(x >= GS1_list[0]) & (x <= GS1_list[1])], conf95[(x >= GS1_list[0]) & (x <= GS1_list[1])], color = '#FF0000')
ax.fill_between(x[(x >= BRS_list[0]) & (x <= BRS_list[1])], conf95[(x >= BRS_list[0]) & (x <= BRS_list[1])], color = '#D930BD')
ax.fill_between(x[(x >= GS2_list[0]) & (x <= GS2_list[1])], conf95[(x >= GS2_list[0]) & (x <= GS2_list[1])], color = '#FF009B')

stemobj = ax.stem(x200, y200, markerfmt = ' ')
plt.setp(stemobj, 'color', 'black')
plt.setp(stemobj, 'zorder', 20)
ax2.set_xlim(0, max(x))
ax.set_xticks(range(0, max(x), 2000))
ax.margins(y = 0)
ax2.set_xticks(range(0, max(x), 2000))
ax2.plot(breaks, [7]*len(breaks), '|', color = 'black')
ax2.set_axis_off()
f.subplots_adjust(wspace=0, hspace=0)

right_side = ax.spines['right']
right_side.set_visible(False)
top_side = ax.spines['top']
top_side.set_visible(False)

ax3.margins(x = 0.01)
ax3.plot(x, x/1000, alpha = 0)
ax3.set_xticks(range(0, max(x), 2000))
ax3.get_yaxis().set_visible(False)


for index, gene in genes.iterrows():
    gene_zorder = 2
    basecolor = '#AFAFAF'
    basewidth = 1000
    headwidth = 1000
    linewidth = 5
    alpha = 1
    if 'brs' in gene['name']:
        basecolor = '#DC6FCA'
    elif gene['name'] == 'glu1':
        basecolor = '#FF6060'
    elif 'glu2' in gene['name']:
        basecolor = '#FF60C1'
    else:
        basewidth = 500
        headwidth = 500
        if gene['end'] - gene['start'] < 1000 or gene['name'] == 'hpr':
            headwidth = 0
            linewidth = 0
            alpha = 0.1
            gene_zorder = 1
    if gene['strand'] == '-' and 'brs' not in gene['name']:
        ax3.arrow(gene['end'], 0, dx = gene['start'] - gene['end'], dy = 0, 
              facecolor = basecolor, length_includes_head = True, 
              width = basewidth, shape = 'full', head_width = headwidth, 
              edgecolor = 'black', linewidth = linewidth, alpha = alpha,
              zorder = gene_zorder)
    elif 'brs' not in gene['name']:
        ax3.arrow(gene['start'], 0, dx = gene['end'] - gene['start'], dy = 0, 
              facecolor = basecolor, length_includes_head = True, 
              width = basewidth, shape = 'full', head_width = headwidth, 
              edgecolor = 'black', linewidth = linewidth, alpha = alpha,
              zorder = gene_zorder)
    if 'glu' in gene['name'] and '_' not in gene['name']:
        leg = ax3.text(gene['end']-(gene['end']-gene['start'])/(1+2/3), -60, 
                       name_dict[gene['name']], fontsize = 18)
        
brs_start =  min(genes[genes['name'] == 'brsu']['start'])
brs_end = max(genes[genes['name'] == 'brsu']['end'])
ax3.arrow(brs_end, 0, dx = brs_start - brs_end, dy = 0, 
      facecolor = '#DC6FCA', length_includes_head = True, 
      width = 1000, shape = 'full', head_width = 1000, 
      edgecolor = 'black', linewidth = 5, alpha = 1)  
leg = ax3.text(brs_end-(brs_end-brs_start)/(1+2/3), -60, 
               'BRS', fontsize = 18)


right_side3 = ax3.spines['right']
right_side3.set_visible(False)
left_side3 = ax3.spines['left']
left_side3.set_visible(False)

ax.set_ylabel('Breakpoints per 200nt window', fontsize = 16)
plt.xlabel('Position in the alignment', fontsize = 16)