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

GH32 = 'all_S2-3'
GH70 = 'all_BRS'

infolder = 'RDP5_results'
infile = f'{infolder}/{GH32}_breakplot.csv'
outfolder = 'RDP5_results'
GH_pos = f'GH_positions/{GH32}_GH_positions_all.tab'
GH_domains = f'GH_positions/{GH32}_GH_domains.tab'


df = pd.read_csv(infile)

#all_GS1_S2-3
S2b = (7006, 10845)
S2a = (10927, 14988)
GS1 = (15232, 19548)
S3 = (26482, 29335)

#other_GS1_S2-3
S2b = (7001, 10840)
S2a = (10922, 14983)
GS1 = (15145, 19501)
S3 = (26482, 29321)

name_dict = {'gtf2a': 'S2a', 'gtf2b': 'S2b', 'gtf3': 'S3', 'glu1': 'GS1',
             'glu2': 'GS2', 'brsu': 'BRS'}

with open(GH_pos) as positions:
    genes = pd.read_csv(positions, delimiter = '\t')
    
with open(GH_domains) as domains:
    dom = pd.read_csv(domains, delimiter = '\t')
    
S2b_list = [inf, 0]
S2a_list = [inf, 0]
GS1_list = [inf, 0]
S3_list = [inf, 0]
for index, row in dom.iterrows():
    if row['start'] < 9000:
        S2b_list[0] = min(S2b_list[0], row['start'])
        S2b_list[1] = max(S2b_list[1], row['end'])
    elif row['start'] < 15000:
        S2a_list[0] = min(S2a_list[0], row['start'])
        S2a_list[1] = max(S2a_list[1], row['end'])
    elif row['start'] < 20000:
        GS1_list[0] = min(GS1_list[0], row['start'])
        GS1_list[1] = max(GS1_list[1], row['end'])
    else:
        S3_list[0] = min(S3_list[0], row['start'])
        S3_list[1] = max(S3_list[1], row['end'])


breakpoints = pd.read_csv(f'{infile}BreakpointPositions.csv')
breaks = breakpoints['Breakpoint position']

x = df['Position in alignment']
y = df[' Recombination breakpoint number (200nt win)']
conf95 = df[' Upper 95% CI']
conf99 = df[' Upper 99% CI']

num = 100

x200 = [x[i] for i in range(len(x)-1) if (i % num == 0) or i == x[len(x)-1]]
y200 = [mean(y[int(num*i):int(num*(i+1))]) for i in range(int(len(y)/num))]
y200.append(mean(y[len(y)-num:]))
# x = [x[i] for i in range(len(x)) if x[i] % num == 0]

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

# x = pd.Series(x)
y = pd.Series(y)
conf95 = pd.Series(conf95)
conf99 = pd.Series(conf99)

x200 = pd.Series(x200)
y200 = pd.Series(y200)

f, (ax2, ax, ax3) = plt.subplots(3, 1, figsize = (18, 8), gridspec_kw={'height_ratios': [1, 8, 2]})
ax.margins(x=0.01)

ax.plot(x[x < S2b[0]], conf95[x < S2b[0]], color = 'grey')
ax.plot(x[(x >= S2b[0]) & (x <= S2b[1])], conf95[(x >= S2b[0]) & (x <= S2b[1])], color = '#FFAA00') #'#F19A00')
ax.plot(x[(x > S2b[1]) & (x < S2a[0])], conf95[(x > S2b[1]) & (x < S2a[0])], color = 'grey')
ax.plot(x[(x >= S2a[0]) & (x <= S2a[1])], conf95[(x >= S2a[0]) & (x <= S2a[1])], color = '#FFAA00')
ax.plot(x[(x > S2a[1]) & (x < GS1[0])], conf95[(x > S2a[1]) & (x < GS1[0])], color = 'grey')
ax.plot(x[(x >= GS1[0]) & (x <= GS1[1])], conf95[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF0000')
ax.plot(x[(x > GS1[1]) & (x < S3[0])], conf95[(x > GS1[1]) & (x < S3[0])], color = 'grey')
ax.plot(x[(x >= S3[0]) & (x <= S3[1])], conf95[(x >= S3[0]) & (x <= S3[1])], color = '#FFD500') #'#EDB700')
ax.plot(x[x > S3[1]], conf95[x > S3[1]], color = 'grey')

ax.plot(x[x < S2b[0]], conf99[x < S2b[0]], color = 'grey', alpha = 0.3)
ax.plot(x[(x >= S2b[0]) & (x <= S2b[1])], conf99[(x >= S2b[0]) & (x <= S2b[1])], color = '#FFAA00', alpha = 0.3)
ax.plot(x[(x > S2b[1]) & (x < S2a[0])], conf99[(x > S2b[1]) & (x < S2a[0])], color = 'grey', alpha = 0.3)
ax.plot(x[(x >= S2a[0]) & (x <= S2a[1])], conf99[(x >= S2a[0]) & (x <= S2a[1])], color = '#FFAA00', alpha = 0.3)
ax.plot(x[(x > S2a[1]) & (x < GS1[0])], conf99[(x > S2a[1]) & (x < GS1[0])], color = 'grey', alpha = 0.3)
ax.plot(x[(x >= GS1[0]) & (x <= GS1[1])], conf99[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF0000', alpha = 0.3)
ax.plot(x[(x > GS1[1]) & (x < S3[0])], conf99[(x > GS1[1]) & (x < S3[0])], color = 'grey', alpha = 0.3)
ax.plot(x[(x >= S3[0]) & (x <= S3[1])], conf99[(x >= S3[0]) & (x <= S3[1])], color = '#FFD500', alpha = 0.3)
ax.plot(x[x > S3[1]], conf99[x > S3[1]], color = 'grey', alpha = 0.3)

ax.fill_between(x[x < S2b[0]], conf95[x < S2b[0]], color = '#AFAFAF')
ax.fill_between(x[(x >= S2b[0]) & (x <= S2b[1])], conf95[(x >= S2b[0]) & (x <= S2b[1])], color = '#FFC800')
ax.fill_between(x[(x > S2b[1]) & (x < S2a[0])], conf95[(x > S2b[1]) & (x < S2a[0])], color = '#AFAFAF')
ax.fill_between(x[(x >= S2a[0]) & (x <= S2a[1])], conf95[(x >= S2a[0]) & (x <= S2a[1])], color = '#FFC800')
ax.fill_between(x[(x > S2a[1]) & (x < GS1[0])], conf95[(x > S2a[1]) & (x < GS1[0])], color = '#AFAFAF')
ax.fill_between(x[(x >= GS1[0]) & (x <= GS1[1])], conf95[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF6060')
ax.fill_between(x[(x > GS1[1]) & (x < S3[0])], conf95[(x > GS1[1]) & (x < S3[0])], color = '#AFAFAF')
ax.fill_between(x[(x >= S3[0]) & (x <= S3[1])], conf95[(x >= S3[0]) & (x <= S3[1])], color = '#FFEF00')
ax.fill_between(x[x > S3[1]], conf95[x > S3[1]], color = '#AFAFAF')

ax.fill_between(x[x < S2b[0]], conf99[x < S2b[0]], color = '#AFAFAF', alpha = 0.3)
ax.fill_between(x[(x >= S2b[0]) & (x <= S2b[1])], conf99[(x >= S2b[0]) & (x <= S2b[1])], color = '#FFC800', alpha = 0.3)
ax.fill_between(x[(x > S2b[1]) & (x < S2a[0])], conf99[(x > S2b[1]) & (x < S2a[0])], color = '#AFAFAF', alpha = 0.3)
ax.fill_between(x[(x >= S2a[0]) & (x <= S2a[1])], conf99[(x >= S2a[0]) & (x <= S2a[1])], color = '#FFC800', alpha = 0.3)
ax.fill_between(x[(x > S2a[1]) & (x < GS1[0])], conf99[(x > S2a[1]) & (x < GS1[0])], color = '#AFAFAF', alpha = 0.3)
ax.fill_between(x[(x >= GS1[0]) & (x <= GS1[1])], conf99[(x >= GS1[0]) & (x <= GS1[1])], color = '#FF6060', alpha = 0.3)
ax.fill_between(x[(x > GS1[1]) & (x < S3[0])], conf99[(x > GS1[1]) & (x < S3[0])], color = '#AFAFAF', alpha = 0.3)
ax.fill_between(x[(x >= S3[0]) & (x <= S3[1])], conf99[(x >= S3[0]) & (x <= S3[1])], color = '#FFEF00', alpha = 0.3)
ax.fill_between(x[x > S3[1]], conf99[x > S3[1]], color = '#AFAFAF', alpha = 0.3)

ax.fill_between(x[(x >= GS1_list[0]) & (x <= GS1_list[1])], conf95[(x >= GS1_list[0]) & (x <= GS1_list[1])], color = '#FF0000')
ax.fill_between(x[(x >= S2a_list[0]) & (x <= S2a_list[1])], conf95[(x >= S2a_list[0]) & (x <= S2a_list[1])], color = '#FFAA00')
ax.fill_between(x[(x >= S3_list[0]) & (x <= S3_list[1])], conf95[(x >= S3_list[0]) & (x <= S3_list[1])], color = '#FFD500')

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

stemobj = ax.stem(x200, y200, markerfmt = ' ') #, color = 'black')
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
    gene_zorder = 1
    basecolor = '#AFAFAF'
    basewidth = 1000
    headwidth = 1000
    linewidth = 5
    alpha = 1
    if 'gtf2' in gene['name']:
        basecolor = '#FFC800'
    elif gene['name'] == 'glu1':
        basecolor = '#FF6060'
    elif 'gtf3' in gene['name']:
        basecolor = '#FFEF00'
    else:
        basewidth = 500
        headwidth = 500
        if gene['end'] - gene['start'] < 1000 or gene['name'] == 'hpr':
            headwidth = 0
            linewidth = 0
            alpha = 0.1
            gene_zorder = 0
    # if gene['name'] != 'hpr':
    ax3.arrow(gene['end'], 0, dx = gene['start'] - gene['end'], dy = 0, 
          facecolor = basecolor, length_includes_head = True, 
          width = basewidth, shape = 'full', head_width = headwidth, 
          edgecolor = 'black', linewidth = linewidth, alpha = alpha,
          zorder = gene_zorder)
    if 'gtf' in gene['name'] or 'glu' in gene['name']:
        leg = ax3.text(gene['end']-(gene['end']-gene['start'])/(1+2/3), -60, 
                       name_dict[gene['name']], fontsize = 18)
    
# ax3.arrow(GS1[0]+GS1[1]-GS1[0], 0.5, dx = -(GS1[1]-GS1[0]), dy = 0, facecolor = '#FF6060',
#           length_includes_head = True, width = 1000, shape = 'full', 
#           head_width = 1000, edgecolor = 'black', linewidth = 5, label = 'GS1')
# leg = ax3.text(GS1[0]+GS1[1]-GS1[0]-2500, -60, 'GS1', fontsize = 18)

# ax3.arrow(S2a[0]+S2a[1]-S2a[0], 0.5, dx = -(S2a[1]-S2a[0]), dy = 0, facecolor = '#FFC800',
#           length_includes_head = True, width = 1000, shape = 'full', 
#           head_width = 1000, edgecolor = 'black', linewidth = 5, label = 'S2a')
# leg = ax3.text(S2a[0]+S2a[1]-S2a[0]-2500, -60, 'S2a', fontsize = 18)

# ax3.arrow(S2b[0]+S2b[1]-S2b[0], 0.5, dx = -(S2b[1]-S2b[0]), dy = 0, facecolor = '#FFC800',
#           length_includes_head = True, width = 1000, shape = 'full', 
#           head_width = 1000, edgecolor = 'black', linewidth = 5, label = 'S2b')
# leg = ax3.text(S2b[0]+S2b[1]-S2b[0]-2500, -60, 'S2b', fontsize = 18)

# ax3.arrow(S3[0]+S3[1]-S3[0], 0.5, dx = -(S3[1]-S3[0]), dy = 0, facecolor = '#FFEF00',
#           length_includes_head = True, width = 1000, shape = 'full', 
#           head_width = 1000, edgecolor = 'black', linewidth = 5, label = 'S3')
# leg = ax3.text(S3[0]+S3[1]-S3[0]-1500, -60, 'S3', fontsize = 18)

right_side3 = ax3.spines['right']
right_side3.set_visible(False)
left_side3 = ax3.spines['left']
left_side3.set_visible(False)

ax.set_ylabel('Breakpoints per 200nt window', fontsize = 16)
plt.xlabel('Position in the alignment', fontsize = 16)