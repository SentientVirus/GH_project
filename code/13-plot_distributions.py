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
from scipy.stats import ks_2samp

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

# =============================================================================
# 2. Create the plot layout
# =============================================================================

with plt.style.context('fivethirtyeight'): #Ass plot style 
    fig, axs = plt.subplots(1, 4, constrained_layout = False, figsize=(18, 3), #Set no. of frames, layout and figure size
                            gridspec_kw={'width_ratios': [1, 1, 1, 1]}) #, #Adjust width ratios of the plots
                            # sharey = 'row')
fig.patch.set_facecolor('white') #Make the background white

plt.subplots_adjust(hspace = 0, top = 1, bottom = 0.05, wspace = 0.4) #Adjust spacing around the plots

# =============================================================================
# 3. Go through the input files and generate the plots
# =============================================================================

for gtype in gtypes: #Loop through gene types
    print(gtype) #Print gene types

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
    
    for file in infiles: #Loop through input files
        with open(file) as tsv: #Open the file
            reader = tsv.readlines()[1:] #Read file contents skipping the first line
            for line in reader: #Loop through file contents
                values = line.split('\t') #Retrieve the tab-separated values as a list
                dN = float(values[2]) #Retrieve dN
                dS = float(values[3]) #Retrieve dS
                w = float(values[4]) #Retrieve dN/dS
                
                #Add the values to the corresponding list if the value is > 0.01 and < 1.5
                if 0.01 < dN < 1.5:
                    dN_list.append(dN)
                if 0.01 < dS < 1.5:
                    dS_list.append(dS)
                if 0.01 < w < 1.5:
                    w_list.append(w)
             
    #Do the same, but for the pairwise comparisons within the gene subtype
    with open(type_comparison) as tsv:
        reader = tsv.readlines()[1:] 
        for line in reader:
            values = line.split('\t')
            dN = float(values[2])
            dS = float(values[3])
            w = float(values[4])
            
            if 0.01 < dN < 1.5:
                type_dN.append(dN)
            if 0.01 < dS < 1.5:
                type_dS.append(dS)
            if 0.01 < w < 1.5:
                type_w.append(w)
                
    #Perform a Kolmogorov-Smirnov test to check if the distributions of core genes and GH70 subtypes are significantly different
    ks_stat = ks_2samp(type_dS, dS_list)
    print(ks_stat)
                
    axs[i].set_facecolor('#FFF9EF') #Color the axis background
    axs[i].grid(False) #Remove grid
    axs[i].set_xlim(0, 1.5) #Set the range of the x axis
    axs[i].set_xticks([0, 0.5, 1, 1.5]) #Set axis ticks
    #Plot the core gene dS values as a histogram with 50 bins
    axs[i].hist(dS_list, bins = 50, range = (0, 1.5), linewidth = 0.5, color = '#6B818C', edgecolor='black', density = False, alpha = 1, zorder = 15)
    
    ax2 = axs[i].twinx() #Create a twin axis
    #Plot a histogram again, this type with the pairwise dS from the gene subtype
    ax2.hist(type_dS, bins=50, range = (0, 1.5), linewidth=0.5, color = '#D8E4FF', edgecolor='black', density = False, alpha = 0.8, zorder = 15)
    ax2.tick_params(axis='both', which='major', labelsize = tick_fontsize, length = 0) #Remove the ticks, but keep the labels
    
    max_y = max(ax2.get_yticks()) #Get the highest tick value
    ax2.set_ylim(0, max_y) #Set the scale to the highest tick value
    #Add text with the KS test results
    ax2.text(1.45, max_y, f'$KS = {ks_stat[0]:.3f}$, $p$-$value < 10⁻⁵$', fontsize = 14, 
             horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)
    
    # axs[i].set_ylim(0, 12.5)
    # ax2.set_ylim(0, 12.5)
    
    if i == 0: #If it is the first plot in the row
        axs[i].set_ylabel('No. of core pairwise ${d_S}$', fontsize = fontsize_subtitle) #Add the label for the first histogram
    elif i == 3: #If it is the last plot in the row
        ax2.set_ylabel('No. of subtype pairwise ${d_S}$', fontsize = fontsize_subtitle, rotation = 270, labelpad = 25) #Add the label for the second histogram
    if i != 3: #If it is not the last plot
        ax2.set_ylabel('') #Remove any label for the second histogram
    
    i += 1 #Increase the count variable by 1 
    
    plt.title(f'${gtype}$', fontsize = fontsize_title) #Add the gene name as plot title
    
[axs[n].tick_params(axis='both', which='major', labelsize = tick_fontsize) for n in range(0, 4)] #Adjust fontsize of the tick labels

plt.savefig(outplot, bbox_inches = 'tight') #Save the image as PNG
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches = 'tight') #Save the image as TIFF
plt.savefig(outplot.replace('.png', '.pdf'), bbox_inches = 'tight') #Save the image as PDF
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches = 'tight') #Save the image as SVG