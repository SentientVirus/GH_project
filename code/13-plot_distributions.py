#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 16:02:26 2025

Test script to create a distribution of dN and dS values for the pairs of
strains that feature GS2.

@author: Marina Mota-Merlo
"""
# =============================================================================
# 0. Import required packages
# =============================================================================

import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import ks_2samp, pearsonr
import numpy as np
from matplotlib.ticker import MultipleLocator

# =============================================================================
# 1. Define output and formatting variables
# =============================================================================

gtypes = ['GS1', 'GS2', 'BRS', 'NGB'] #Gene subtypes to plot
outplot = os.path.expanduser('~') + '/GH_project/plots/distributions/GH70.png' #Path to output image

#Create output directory if it doesn't exist
if not os.path.exists(os.path.dirname(outplot)):
    os.makedirs(os.path.dirname(outplot))

tick_fontsize = 14 #Fontsize (ax tick labels)
fontsize = 16 #Fontsize (legend)
fontsize_title = 24 #Fontsize (titles)
fontsize_subtitle = 20 #Fontsize (ax labels)
plt.rcParams['font.family'] = 'Arial' #Font type

i = 0 #Variable to loop through the axes

dNdS = {}
avg_dNdS = {}

ticks = [0.0, 0.5, 1.0, 1.5]
tick_labels = ['0.0', '0.5', '1.0', r'$\geq1.5$']

# =============================================================================
# 2. Create the plot layout
# =============================================================================

with plt.style.context('fivethirtyeight'): #Ass plot style 
    fig, axs = plt.subplots(2, 4, constrained_layout = False, figsize=(18, 6), #Set no. of frames, layout and figure size
                            gridspec_kw={'width_ratios': [1, 1, 1, 1]}) #Adjust width ratios of the plots

fig.patch.set_facecolor('white') #Make the background white

plt.subplots_adjust(hspace = 0.2, top = 1, bottom = 0.03, wspace = 0.4) #Adjust spacing around the plots

#Create lists to color the plot under x == y
x_line = [x for x in range(3)]
y_line = [y for y in range(3)]
where_param = [y_line[i] == x_line[i] for i in range(len(y_line))] #Comparison between x and y

fig.text(0.51, -0.08, 'Pairwise gene ${d_S}$', ha = 'center', fontsize = fontsize_subtitle) #Add title of the x axis of dN vs dS scatter plots

# =============================================================================
# 3. Go through the input files and generate the plots
# =============================================================================

for gtype in gtypes: #Loop through gene types
    print(gtype) #Print gene types
    dNdS[gtype] = []

    indir = os.path.expanduser('~') + f'/GH_project/all_core/{gtype}/results' #Path to input directory
    
    #List of input files with core gene pairwise metrics
    infiles = [f'{indir}/{subdir}/dNdS.tsv' for subdir in os.listdir(indir) if os.path.isdir(f'{indir}/{subdir}') and subdir.startswith('A1401')]
    #Input file with the pairwise metrics for each gene
    type_comparison = os.path.expanduser('~') + f'/GH_project/results/{gtype}/dNdS.tsv'
    
    
    #Create lists to store core dN, dS and dN/dS (w)
    dN_list = []
    dS_list = []
    w_list = []
    
    #Create lists to store the gene subtype dN, dS and dN/dS (w)
    type_dN = []
    type_dS = []
    type_w = []
    filtered_dN = []
    filtered_dS = []
    
    
    for file in infiles: #Loop through input files
        with open(file) as tsv: #Open the file
            reader = tsv.readlines()[1:] #Read file contents skipping the first line
            for line in reader: #Loop through file contents
                values = line.split('\t') #Retrieve the tab-separated values as a list
                dN = float(values[2]) #Retrieve dN
                dS = float(values[3]) #Retrieve dS
                w = float(values[4]) #Retrieve dN/dS

                #Add the values to the corresponding list if the value is > 0.01 and < 1.5
                if dN < 1.5: #If the dN falls in the desired range
                    dN_list.append(dN) #Add it to the list
                else:
                    dN_list.append(1.5)
                if dS < 1.5: #Same, but for the dS
                    dS_list.append(dS)
                else:
                    dS_list.append(1.5)
                w_list.append(w) #Retrieve the dN/dS
             
    #Do the same, but for the pairwise comparisons within the gene subtype
    with open(type_comparison) as tsv:
        reader = tsv.readlines()[1:] 
        for line in reader:
            values = line.split('\t')
            dN = float(values[2])
            dS = float(values[3])
            w = float(values[4])
            
            if dN < 1.5:
                type_dN.append(dN)
            else:
                type_dN.append(1.5)
            if dS < 1.5:
                type_dS.append(dS)
                filtered_dN.append(dN)
                filtered_dS.append(dS)
            else:
                type_dS.append(1.5)
            type_w.append(w)
    
    print(f'Length of core gene lists: dS - {len(dS_list)}, dN - {len(dN_list)} - w - {len(w_list)}')
    print(f'Length of {gtype} lists: dS - {len(type_dS)}, dN - {len(type_dN)} - w - {len(type_w)}')

    core = dS_list
    subtype = type_dS
    max_x = 1.5
    
    #Perform a Kolmogorov-Smirnov test to check if the distributions of core genes and GH70 subtypes are significantly different
    ks_stat = ks_2samp(core, subtype)
    print(ks_stat)
              
    for k in range(2):
        axs[k, i].set_facecolor('#FFF9EF') #Color the axis background
        axs[k, i].grid(False) #Remove grid
        axs[k, i].set_xlim(0, max_x+0.02) #Set the range of the x axis


    #Plot the core gene dS values as a histogram with 50 bins
    axs[1, i].hist(core, bins = 50, range = (0, max_x), linewidth = 0.5, color = '#6B818C', edgecolor = 'black', density = False, alpha = 1, zorder = 15)
    axs[1, i].set_xticks(ticks, labels = tick_labels)

    ax2 = axs[1, i].twinx() #Create a twin axis
    #Plot a histogram again, this type with the pairwise dS from the gene subtype
    ax2.hist(subtype, bins = 50, range = (0, max_x), linewidth = 0.5, color = '#D8E4FF', edgecolor='black', density = False, alpha = 0.8, zorder = 15)
    ax2.tick_params(axis='both', which='major', labelsize = tick_fontsize, length = 0) #Remove the ticks, but keep the labels
    
    max_y = max(ax2.get_yticks()) #Get the highest tick value
    ax2.set_ylim(0, max_y) #Set the scale to the highest tick value
    x_pos = max_x*1.3/1.5

    #Add text with the KS test results
    if ks_stat[1] < 1E-10:
        ax2.text(x_pos, 0.9*max_y, f'$KS = {ks_stat[0]*ks_stat.statistic_sign:.3f}$\n' + r'$p < 10^{-10}$', fontsize = 14, 
                horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)
    elif ks_stat[1] < 1E-5:
        ax2.text(x_pos, 0.9*max_y, f'$KS = {ks_stat[0]*ks_stat.statistic_sign:.3f}$\n' + r'$p < 10^{-5}$', fontsize = 14, 
                horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)
    elif ks_stat[1] < 1E-3:
        ax2.text(x_pos, 0.9*max_y, f'$KS = {ks_stat[0]*ks_stat.statistic_sign:.3f}$\n' + r'$p < 10^{-3}$', fontsize = 14, 
                horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)
    else:
        ax2.text(x_pos, 0.9*max_y, f'$KS = {ks_stat[0]*ks_stat.statistic_sign:.3f}$\n$p = {ks_stat[1]:.3f}$', fontsize = 14, 
                horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)

    if i == 0:
        axs[1, i].set_ylabel('No. of core pairwise ${d_S}$', fontsize = fontsize_subtitle) #Add the label for the first histogram
        axs[0, i].set_ylabel('Pairwise gene ${d_N}$', fontsize = fontsize_subtitle) #Add the label for the first histogram, second row
    elif i == 3: #If it is the last plot in the row
        ax2.set_ylabel('No. of gene pairwise ${d_S}$', fontsize = fontsize_subtitle, rotation = 270, labelpad = 25) #Add the label for the second histogram
    if i != 3: #If it is not the last plot
        ax2.set_ylabel('') #Remove any label for the second histogram
        
    axs[0, i].scatter(type_dS, type_dN, alpha = 0.8, s = 20, c = ['#D8E4FF']*len(type_dS), edgecolors = 'black', zorder = 10) #Add scatter plot (dS vs dN)
    r1 = pearsonr(filtered_dS, filtered_dN) #Calculate the Pearson correlation coefficient
    slope1, intercept1 = np.polyfit(filtered_dS, filtered_dN, 1) #Perform linear regression
    print(r1[1])
    if r1[1] < 1E-10:
        pval = ' $p < 10^{-10}$'
    elif r1[1] < 1E-5:
        pval = ' $p < 10^{-5}$'
    elif r1[1] < 1E-3:
        pval = ' $p < 10^{-3}$'
    else: pval = ' $p = {r1[1]:.3f}$'
        
    
    #Add linear regression to the plot
    axs[0, i].axline(xy1=(0, intercept1), color = 'darkorange', slope=slope1, label=f'$y = {slope1:.3f}x {intercept1:+.3f}$\n$R² = {r1[0]**2:.3f},$' + pval, zorder = 5)
    axs[0, i].legend(loc = 'upper right', frameon = False) #Add legend to the linear regression at the upper right corner and without a frame
    
    #Format the scatter plot axes
    axs[0, i].plot([0, 1], [0, 1], color = 'lightgrey') #Plot a line where x = y
    axs[0, i].set_ylim(0, 0.25) #Set y ax limits (based on the highest dN in the dataset)
    axs[0, i].fill_between(y_line, x_line, where = where_param, 
                           interpolate = False, color = 'white', alpha = 0.5) #Fill everything under the x = y line
    axs[0, i].set_xticks(ticks, labels = tick_labels)
    
    minor_ticks = MultipleLocator(0.25) #Create a locator for minor ticks
    axs[0, i].xaxis.set_minor_locator(minor_ticks) #Apply it to the scatterplot ax
    axs[0, i].grid(color = '#F7F0D5', linestyle = '-', linewidth = 2, which = 'both') #Add a grid to the plot, on both major and minor ticks
    axs[0, i].set_title(f'${gtype.replace("R", "r")}$', fontsize = fontsize_title) #Add the gene name as plot title

    i += 1 #Increase the count variable by 1
    
[spine.set_edgecolor('black') for n in range(0, 4) for m in range(0, 2) for spine in axs[m, n].spines.values()]
[axs[m, n].tick_params(axis = 'both', which = 'major', labelsize = tick_fontsize) for n in range(0, 4) for m in range(0, 2)] #Adjust fontsize of the tick labels

#Add scatter plot legend
subtitle1 = Line2D([], [], marker = '', color = 'black', label = 'Scatter plots', 
                    ms = 0, ls = 'none') #Title of the scatter plot legends

sub_patch = Line2D([0], [0], marker = 'o', color = 'black', label = '(x, y)',
                        markerfacecolor = '#D8E4FF', markersize = 10, 
                        alpha = 0.8, linewidth = 0) #Circular marker for the dots in the scatter plot

line_patch = Line2D([0, 1], [1, 0], marker = '$/$', color = 'lightgrey', 
                    label = 'x = y', markersize = 20, alpha = 1, linewidth = 0) #Diagonal line marker

line_patch2 = Line2D([0, 1], [1, 0], marker = '$/$', color = 'darkorange', 
                    label = 'Linear\nregressions', markersize = 20, alpha = 1, linewidth = 0) #Diagonal line marker

patches = [subtitle1, sub_patch, line_patch, line_patch2] #List with all elements of the legend

legend = fig.legend(handles = patches,  handlelength = 0.6, handleheight = 1.5, 
                    loc = 'upper right', framealpha = 0, 
                    frameon = False, fontsize = fontsize, bbox_to_anchor = (1.13, 1)) #Add legend to plot
           
#Set legend headers to bold
legend.get_texts()[0].set_fontweight('bold')
legend.get_texts()[0].set_fontsize(fontsize_subtitle)
legend.get_texts()[0].set_position((-25, 0))

#Adjust legend font size
for i in range(1, 3):
    legend.get_texts()[i].set_fontsize(fontsize_subtitle) #Increase the legend fontsize
legend._legend_box.align = 'left' #Align legend text to the left

#Add histogram plot legend
subtitle2 = Line2D([], [], marker = '', color = 'black', label = 'Histograms   ', 
                    ms = 0, ls = 'none') #Title of the histogram plot legends

sub_patch2 = Line2D([0], [0], marker = 's', color = 'black', label = 'Core ${d_S}$',
                        markerfacecolor = '#6B818C', markersize = 10, 
                        alpha = 0.8, linewidth = 0) #Square marker for core gene histograms

sub_patch3 = Line2D([0], [0], marker = 's', color = 'black', label = 'Gene ${d_S}$',
                        markerfacecolor = '#D8E4FF', markersize = 10, 
                        alpha = 0.8, linewidth = 0) #Square marker for gene subtype histograms


patches2 = [subtitle2, sub_patch2, sub_patch3] #List with all elements of the legend

legend2 = fig.legend(handles = patches2,  handlelength = 0.6, handleheight = 1.5, 
                    loc = 'upper right', framealpha = 0, 
                    frameon = False, fontsize = fontsize, bbox_to_anchor = (1.13, 0.5)) #Add legend to plot
           
#Set legend headers to bold
legend2.get_texts()[0].set_fontweight('bold')
legend2.get_texts()[0].set_fontsize(fontsize_subtitle)
legend2.get_texts()[0].set_position((-25, 0))

#Adjust legend font size
for i in range(1, 3):
    legend2.get_texts()[i].set_fontsize(fontsize_subtitle) #Increase the legend fontsize
legend2._legend_box.align = 'left' #Align legend text to the left

plt.savefig(outplot, bbox_inches = 'tight') #Save the image as PNG
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches = 'tight') #Save the image as TIFF
plt.savefig(outplot.replace('.png', '.pdf'), bbox_inches = 'tight') #Save the image as PDF
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches = 'tight') #Save the image as SVG
