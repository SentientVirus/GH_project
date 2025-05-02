#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 14:24:33 2024

@author: marina
"""
import os
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import pearsonr

# =============================================================================
# 1. Define the paths to input and output files
# =============================================================================
workdir = os.path.expanduser('~') + '/GH_project' #Working directory
pairwise_ds = f'{workdir}/results' #Directory with input pairwise dS data
core_ds = f'{workdir}/all_core/results/global/core_pairwise_metrics.tsv' #'core_genes/core_pairwise_metrics.tsv' #File with core dS data
outplot = f'{workdir}/plots/scatter/GH_types_suppl_scatterplot.png' #Output file
coredir = f'{workdir}/all_core'

GH_types = ['GS1', 'GS2', 'BRS', 'NGB'] #, 'S2a', 'S3'] #Gene types to be plotted
to_exclude = ['A1001', 'A1404']

tick_fontsize = 14 #Fontsize (ax tick labels)
fontsize = 16 #Fontsize (legend)
fontsize_title = 24 #Fontsize (titles)
fontsize_subtitle = 20 #Fontsize (ax labels)
plt.rcParams['font.family'] = 'Arial'

dN_ylims = {0: (0, 0.015), 1: (0, 0.015), 2: (0, 0.015), 3: (0, 0.02)}
dN_xlims = {0: (0, 0.15), 1: (0, 0.2), 2: (0, 0.25), 3: (0, 0.025)}

core_ylims = {}

def replace_strain_name(locus_tag):
    """Function to change locus tags that do not correspond with strain names
    to strain names. The function also removes RS from RefSeq locus tags to
    make the final tags that are plotted more readable. It also removes AKU
    from the non-RefSeq locus tags, as it indicates the species, which is
    the same for all the strains, hence redundant.
    
    Parameters
    ----------
    locus_tag : str
        The locus tag to be overwritten."""
     
    #Modify locus tags to fit strain names    
    locus_tag = locus_tag.replace('LDX55', 'IBH001').replace('APS55', 'MP2').replace('K2W83', 'DSM').replace('RS', '').replace('FHON', 'fhon').replace('AKU', '')
    
    #Get strain from locus tag
    strain = locus_tag.split('_')[0] #Get the text before _ in the locus tag
    if strain.startswith('H') and '-' not in strain: #If the strain name starts by H (i.e. H4B402J)
        strain = strain[:4] + '-' + strain[4:] #Add - to the name (i.e. H4B4-02J)
        
    return strain

# =============================================================================
# Create a dictionary with each pair of strains as keys and the average
# pairwise dS as a value
# =============================================================================
if not os.path.exists(os.path.dirname(outplot)): #If the output directory doesn't exits
    os.makedirs(os.path.dirname(outplot)) #Create the output directory
    
core_ds_dict = {} #Create a dictionary to store the pairwise core dS values
core_dn_dict = {}
with open(core_ds) as core: #Open the file with core dS values
    df = pd.read_csv(core, sep = '\t') #Load the file as a dataframe
    for index, row in df.iterrows(): #Loop through the rows in the dataframe
        strain1 = replace_strain_name(row['strain1']) #Retrieve the first strain in the pair
        strain2 = replace_strain_name(row['strain2']) #Retrieve the second strain in the pair
        if not strain1[0] == 'L' and not strain2[0] == 'L': #Skip strains with names starting with L
            core_ds_dict[tuple(sorted((strain1, strain2)))] = row['mean_dS'] #Use sorted strain names as key and core pairwise dS as value
            core_dn_dict[tuple(sorted((strain1, strain2)))] = row['mean_dN'] #Use sorted strain names as key and core pairwise dN as value
        
#Create a plot with two rows and six columns
with plt.style.context('fivethirtyeight'):
    fig, axs = plt.subplots(2, 8, constrained_layout = False, figsize=(36, 9),
                            gridspec_kw={'width_ratios': [1, 0.75, 1, 0.75, 1, 0.75, 1, 0.75]})
    
fig.patch.set_facecolor('white')

plt.subplots_adjust(hspace = 0.3, top = 1, bottom = 0.07, wspace = 0.25)

axs[0, 0].set_ylabel('${d_N}$ of each core gene', fontsize = fontsize_subtitle) #Add a title to the y axis of the core dN vs dS scatter plot
axs[1, 0].set_ylabel('Core ${d_N}$', fontsize = fontsize_subtitle) #Add a title to the y axis of the core dN vs pairwise dN scatter plot
fig.text(0.51, 0.54, '${d_S}$ of each core gene', ha = 'center', fontsize = fontsize_subtitle) #Add title of the x axis of core dN vs dS scatter plots
fig.text(0.51, 0.01, '${d_N}$', ha = 'center', fontsize = fontsize_subtitle) #Add title of the x axis of core dN vs dN scatter plots

for i in range(len(GH_types)): #Loop through gene types
    GH = GH_types[i] #Name of the gene type
    core_dir = f'{coredir}/{GH}/results'
    axs[0, i].set_title(GH, fontsize = fontsize_title) #Use the name of the gene as the title of the scatter plots
    if i == 0:
        color = '#FF707C'
    elif i == 1:
        color = '#FFE570'
    elif i == 2:
        color = '#84E8A7'
    elif i == 3:
        color = '#336E74'
    
    label_list = [] #Create an empty list to store labels (not used)
    # x_list = [] #Create a list to store pairwise dS values
    y_list = [] #Create a list to store core pairwise dS values
    dN_list = [] #List to store pairwise dN
    dN_y_list = [] #List to store core pairwise dN values
    
    with open(f'{pairwise_ds}/{GH}/dNdS.tsv') as gene: #Open file with dS values
        df2 = pd.read_csv(gene, sep = '\t') #Load file as a dataframe
        for index2, row2 in df2.iterrows(): #Loop through rows in the dataframe
            label = f'{row2["locus1"]}-{row2["locus2"]}' #Create a label with the locus tags in the comparison
            x = row2['dS'] #Get pairwise dS
            dN = row2['dN']
            strain1 = replace_strain_name(row2['locus1']) #Get the name of the first strain
            strain2 = replace_strain_name(row2['locus2']) #Get the name of the second strain
            if strain1 == strain2: #If both genes are in the same strain
                y = 0 #Core dS is 0
                dN_y = 0
            else: 
                y = core_ds_dict[tuple(sorted((strain1, strain2)))] #Otherwise, store core dS in a dictionary
                dN_y = core_dn_dict[tuple(sorted((strain1, strain2)))]
            
            # #Set all dS values above 1.5 to 1.5
            # if x > 1.5:
            #     x = 1.5
            # if y > 1.5:
            #     y = 1.5
                
            # if x <= 1.5 and y <= 1.5:
            label_list.append(label) #Add label to the label list
            # x_list.append(x) #Add x to the pairwise dS list
            # y_list.append(y) #Add y to the core pairwise dS list
            dN_list.append(dN)
            dN_y_list.append(dN_y)
            
    
    x_pos = [x for x in range(120)]
    y_pos = [y for y in range(120)]
    
    dS_corelist = [] #Create an empty list to store core dS values
    dN_corelist = [] #Create an empty list to store core dN values
    high_dS_corelist = []
    high_dN_corelist = []
    for core_gene in os.listdir(core_dir): #Loop through files in the input directory
        if os.path.isdir(f'{core_dir}/{core_gene}'):
            indir = f'{core_dir}/{core_gene}'
            for file in os.listdir(indir):
                if file == 'dNdS.tsv': #If the file includes dN and dS calues
                    print(f'{indir}/{file}')
                    df = pd.read_csv(f'{indir}/{file}', sep = '\t') #Read the file as a dataframe
                    for index, row in df.iterrows(): #Loop through the dataframe
                        check = any(strain in row['locus1'] for strain in to_exclude) or any(strain in row['locus2'] for strain in to_exclude) #Check that basal strains are not present
                        if not check and not (GH == 'GS1' and ('A1003' in row['locus1'] or 'A1003' in row['locus2'])): #Check that gene GS4 is not present
                            dS = row['dS'] #min(1.5, row['dS']) #Get dS (set maximum of 1.5)
                            dN = row['dN'] #min(1.5, row['dN']) #Get dN (set maximum of 1.5)
                            if dS <= 1.5:
                                dS_corelist.append(dS)
                                dN_corelist.append(dN)
                            else:    
                                high_dS_corelist.append(dS) #Append dS to list
                                high_dN_corelist.append(dN) #Append dN to list
    
    axs[0, i*2].scatter(dS_corelist, dN_corelist, alpha = 0.8, s = 20, c = [color]*len(dS_corelist), edgecolors = 'black', zorder = 10)
    r1 = pearsonr(dS_corelist, dN_corelist)
    slope1, intercept1 = np.polyfit(dS_corelist, dN_corelist, 1)
    axs[0, i*2].axline(xy1=(0, intercept1), color = 'k', slope=slope1, label=f'$y = {slope1:.3f}x {intercept1:+.3f}$\n$R² = {r1[0]**2:.3f}, p$-$value = {r1[1]:.3f}$', zorder = 5)
    axs[0, i*2].legend(loc = 'upper right', frameon = False)
    axs[0, i*2].plot([0, 1], [0, 1], color = 'lightgrey')
    axs[0, i*2].set_xlim(0, 1.5)
    axs[0, i*2].set_ylim(0, 0.5)
    axs[0, i*2].fill_between(y_pos, x_pos, where=y_pos == x_pos, interpolate=False, color='white', alpha = 0.5)
    
    axs[0, i*2+1].scatter(high_dS_corelist, high_dN_corelist, alpha = 0.8, s = 20, c = [color]*len(high_dS_corelist), edgecolors = 'black', zorder = 10)
    r1_2 = pearsonr(high_dS_corelist, high_dN_corelist)
    slope1_2, intercept1_2 = np.polyfit(high_dS_corelist, high_dN_corelist, 1)
    axs[0, i*2+1].axline(xy1=(0, intercept1_2), color = 'k', slope=slope1_2, label=f'$y = {slope1_2:.3f}x {intercept1_2:+.3f}$\n$R² = {r1_2[0]**2:.3f}, p$-$value = {r1_2[1]:.3f}$', zorder = 5)
    axs[0, i*2+1].legend(loc = 'upper right', frameon = False)
    axs[0, i*2+1].plot([0, 1], [0, 1], color = 'lightgrey')
    axs[0, i*2+1].set_xlim(1.5, 80)
    axs[0, i*2+1].set_ylim(0, 0.5)
    axs[0, i*2+1].fill_between(y_pos, x_pos, where=y_pos == x_pos, interpolate=False, color='white', alpha = 0.5)
    
    axs[1, i*2].scatter(dN_list, dN_y_list, alpha = 0.8, s = 20, c = [color]*len(dN_list), edgecolors = 'black', zorder = 10)
    r2 = pearsonr(dN_list, dN_y_list)
    slope2, intercept2 = np.polyfit(dN_list, dN_y_list, 1)
    axs[1, i*2].axline(xy1=(0, intercept2), color = 'k', slope=slope2, label=f'$y = {slope2:.3f}x {intercept2:+.3f}$\n$R² = {r2[0]**2:.3f}, p$-$value = {r2[1]:.3f}$', zorder = 5)
    axs[1, i*2].legend(loc = 'upper right', frameon = False)
    axs[1, i*2].plot([0, 1], [0, 1], color = 'lightgrey')
    axs[1, i*2].set_xlim(dN_xlims[i])
    axs[1, i*2].set_ylim(dN_ylims[i])
    axs[1, i*2].fill_between(y_pos, x_pos, where=y_pos == x_pos, interpolate=False, color='white', alpha = 0.5)
    axs[1, i*2+1].axis('off')
    
    
[axs[m, k].tick_params(axis='both', which='major', labelsize = tick_fontsize) for m in range(0, 2) for k in range(0, 4)] #Increase tick label font size

# Add global legend
subtitle1 = Line2D([], [], marker = '', color = 'black', label = 'Scatter plots', 
                    ms = 0, ls = 'none') #Title of the scatter plot legends

sub_patch = Line2D([0], [0], marker = 'o', color = 'black', label = 'Parwise ${d_S}$',
                        markerfacecolor = 'slategrey', markersize = 10, 
                        alpha = 0.8, linewidth = 0) #Circular marker in blue

line_patch = Line2D([0, 1], [1, 0], marker = '$/$', color = 'lightgrey', 
                    label = 'x = y', markersize = 20, alpha = 1, linewidth = 0) #Diagonal line marker



patches = [subtitle1, sub_patch, line_patch] #, subtitle2, rect_sub, rect_over, line_hist] #List with all elements of the legend

legend = fig.legend(handles = patches,  handlelength = 0.6, handleheight = 1.5, 
                    loc = 'upper right', framealpha = 0, 
                    frameon = False, fontsize = fontsize, bbox_to_anchor = (1.03, 1)) #Add legend to plot
           
#Set legend headers to bold
legend.get_texts()[0].set_fontweight('bold')
legend.get_texts()[0].set_fontsize(fontsize_title)
legend.get_texts()[0].set_position((-20, 0))

#Adjust legend font size
for i in range(1, 3):
    legend.get_texts()[i].set_fontsize(fontsize_subtitle)
# # legend.get_texts()[4].set_fontweight('bold')
# # legend.get_texts()[4].set_fontsize(fontsize_title)
# # legend.get_texts()[4].set_position((-20, 0))
legend._legend_box.align = 'left' #Align legend text to the left
    
plt.savefig(outplot, bbox_inches = 'tight') #Save the image as PNG
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches = 'tight') #Save the image as TIFF
plt.savefig(outplot.replace('.png', '.pdf'), bbox_inches = 'tight') #Save the image as PDF
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches = 'tight') #Save the image as SVG
    