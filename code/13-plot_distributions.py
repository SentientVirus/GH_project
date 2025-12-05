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
from math import sqrt

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
    fig, axs = plt.subplots(3, 4, constrained_layout = False, figsize=(18, 9), #Set no. of frames, layout and figure size
                            gridspec_kw={'width_ratios': [1, 1, 1, 1]}) #, #Adjust width ratios of the plots
                            # sharey = 'row')
fig.patch.set_facecolor('white') #Make the background white

plt.subplots_adjust(hspace = 0.2, top = 1, bottom = 0.03, wspace = 0.4) #Adjust spacing around the plots

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
                dS_check = False #Boolean to check if the dS falls in the desired range
                dN_check = False #Boolean to check if the dN falls in the desired range
                
                #Add the values to the corresponding list if the value is > 0.01 and < 1.5
                if 0.01 < dN < 1.5: #If the dN falls in the desired range
                    dN_list.append(dN) #Add it to the list
                    dN_check = True #Set the dN boolean to true
                if 0.01 < dS < 1.5: #Same, but for the dS
                    dS_list.append(dS)
                    dS_check = True
                if dN_check and dS_check: #If both the dN and dS fall in the desired range
                    w_list.append(w) #Retrieve the dN/dS
             
    #Do the same, but for the pairwise comparisons within the gene subtype
    with open(type_comparison) as tsv:
        reader = tsv.readlines()[1:] 
        for line in reader:
            values = line.split('\t')
            dN = float(values[2])
            dS = float(values[3])
            w = float(values[4])
            check1 = False #Boolean to check if the dN falls within the desired range
            check2 = False #Boolean to check if the dS falls within the desired range
            
            if 0.01 < dN < 1.5:
                type_dN.append(dN)
                check1 = True
            if 0.01 < dS < 1.5:
                type_dS.append(dS)
                check2 = True
            if check1 and check2:
                type_w.append(w)
    
    #Test to standardize the values for core genes    
    mean_core = sum(dS_list)/len(dS_list)
    sd_core = [(n-mean_core)**2 for n in dS_list]
    sd_core = sqrt(sum(sd_core)/(len(sd_core)-1))
    std_core_dS = [(n-mean_core)/sd_core for n in dS_list]

    #Same, but for gene comparisons
    mean_type = sum(type_dS)/len(type_dS)
    sd_type = [(n-mean_type)**2 for n in type_dS]
    sd_type = sqrt(sum(sd_type)/(len(sd_type)-1))
    std_type_dS = [(n-mean_type)/sd_type for n in type_dS]

    for k in range(3):
        print(k, i)
        if k == 0:
            core = dS_list
            subtype = type_dS
        elif k == 1:
            core = dN_list
            subtype = type_dN
        elif k == 2:
            core = w_list
            subtype = type_w


        #Perform a Kolmogorov-Smirnov test to check if the distributions of core genes and GH70 subtypes are significantly different
        ks_stat = ks_2samp(core, subtype)
        print(ks_stat)
                    
        axs[k, i].set_facecolor('#FFF9EF') #Color the axis background
        axs[k, i].grid(False) #Remove grid
        axs[k, i].set_xlim(0, 1.5) #Set the range of the x axis
        axs[k, i].set_xticks([0, 0.5, 1, 1.5]) #Set axis ticks

        #Plot the core gene dS values as a histogram with 50 bins
        axs[k, i].hist(core, bins = 50, range = (0, 1.5), linewidth = 0.5, color = '#6B818C', edgecolor='black', density = False, alpha = 1, zorder = 15)

        ax2 = axs[k, i].twinx() #Create a twin axis
        #Plot a histogram again, this type with the pairwise dS from the gene subtype
        ax2.hist(subtype, bins=50, range = (0, 1.5), linewidth=0.5, color = '#D8E4FF', edgecolor='black', density = False, alpha = 0.8, zorder = 15)
        ax2.tick_params(axis='both', which='major', labelsize = tick_fontsize, length = 0) #Remove the ticks, but keep the labels
        
        max_y = max(ax2.get_yticks()) #Get the highest tick value
        ax2.set_ylim(0, max_y) #Set the scale to the highest tick value

        #Add text with the KS test results
        if ks_stat[1] < 1E-10:
            ax2.text(1.45, max_y, f'$KS = {ks_stat[0]:.3f}$, ' + r'$p$-$value < 10^{-10}$', fontsize = 14, 
                    horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)
        elif ks_stat[1] < 1E-5:
            ax2.text(1.45, max_y, f'$KS = {ks_stat[0]:.3f}$, ' + r'$p$-$value < 10^{-5}$', fontsize = 14, 
                    horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)
        elif ks_stat[1] < 1E-3:
            ax2.text(1.45, max_y, f'$KS = {ks_stat[0]:.3f}$, ' + r'$p$-$value < 10^{-3}$', fontsize = 14, 
                    horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)
        else:
            ax2.text(1.45, max_y, f'$KS = {ks_stat[0]:.3f}$, $p$-$value = {ks_stat[1]:.3f}$', fontsize = 14, 
                    horizontalalignment = 'right', verticalalignment = 'top', zorder = 30)

        if k == 0: #If it is the first row
            plt.title(f'${gtype.replace("R", "r")}$', fontsize = fontsize_title) #Add the gene name as plot title

        if i == 0 and k == 0: #If it is the first plot in the row
            axs[k, i].set_ylabel('No. of core pairwise ${d_S}$', fontsize = fontsize_subtitle) #Add the label for the first histogram
        elif i == 3 and k == 0: #If it is the last plot in the row
            ax2.set_ylabel('No. of gene pairwise ${d_S}$', fontsize = fontsize_subtitle, rotation = 270, labelpad = 25) #Add the label for the second histogram
        elif i == 0 and k == 1:
            axs[k, i].set_ylabel('No. of core pairwise ${d_N}$', fontsize = fontsize_subtitle) #Add the label for the first histogram
        elif i == 3 and k == 1:
            ax2.set_ylabel('No. of gene pairwise ${d_N}$', fontsize = fontsize_subtitle, rotation = 270, labelpad = 25) #Add the label for the first histogram
        elif i == 0 and k == 2:
            axs[k, i].set_ylabel(r'No. of core pairwise $\omega$', fontsize = fontsize_subtitle) #Add the label for the first histogram
        elif i == 3 and k == 2:
            ax2.set_ylabel(r'No. of gene pairwise $\omega$', fontsize = fontsize_subtitle, rotation = 270, labelpad = 25) #Add the label for the first histogram
        if i != 3: #If it is not the last plot
            ax2.set_ylabel('') #Remove any label for the second histogram

    i += 1 #Increase the count variable by 1
    
[axs[m, n].tick_params(axis='both', which='major', labelsize = tick_fontsize) for n in range(0, 4) for m in range(0, 3)] #Adjust fontsize of the tick labels

plt.savefig(outplot, bbox_inches = 'tight') #Save the image as PNG
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches = 'tight') #Save the image as TIFF
plt.savefig(outplot.replace('.png', '.pdf'), bbox_inches = 'tight') #Save the image as PDF
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches = 'tight') #Save the image as SVG
