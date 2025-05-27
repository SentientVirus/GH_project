#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 16:54:02 2023

This is a plot to create matrices for genes GS1-2, BRS and NGB, where
the dS pairwise values are plotted on the top half and the dN pairwise
values are plotted on the bottom half.

@author: Marina Mota-Merlo
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from matplotlib.transforms import ScaledTranslation

# =============================================================================
# 1. Set the name of input and output files and define input variables
# =============================================================================
workdir = os.path.expanduser('~') + '/GH_project' #Project directory
outdir = f'{workdir}/plots/heatmaps' #Directory to save the outputs
outplot = f'{outdir}/GH70_dNdS_matrix.svg' #Path to save output file

gene_no = {'GS1':  33, 'GS2': 14, 'BRS': 17, 'NGB': 27} #Number of pairs per gene

BRSb = ['H3B101A_13250', 'H3B203J_13380', 'H4B204J_13340', 'IBH001_06330'] #List of BRSb genes

vmin = 0 #Minimum dN or dS
vmax = 1.5 #Maximum (saturated) dN or dS
fig_width = 34 #Width of the output figure
fig_height = 30 #Height of the output figure

replace_str = {'K2W83_RS': 'DSMZ_', 'APS55_RS': 'MP2_', 'LDX55_': 'IBH001_'} #Dictionary to edit the locus tags

# =============================================================================
# 2. Create output directory, output figure and colorbar
# =============================================================================
if not os.path.exists(outdir): #Create the output directory if it doesn't exist
    os.makedirs(outdir)

# Create output figure (five plots, where the last is 14 times narrower)
fig, ax = plt.subplot_mosaic([['A', 'B', ''], ['C', 'D', '']], 
                             figsize = (fig_width, fig_height),
                             gridspec_kw = {'width_ratios': [14, 14, 1]})

plt.style.use('seaborn-pastel') #Add style to plot
plt.subplots_adjust(hspace = 0.35, wspace = 0.1) #Adjust space between subplots

cmap = mpl.cm.magma #Get the colormap for the plot
norm = mpl.colors.Normalize(vmin = vmin, vmax = vmax) #Set the maximum and minimum of the colorbar 

colbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap), #Create a color bar
             cax = ax[''], orientation = 'vertical', #Assign it to the narrow ax and make it vertical
             ticks = np.arange(0, 1.60, 0.1)) #Add ticks at 0.1 intervals until 1.5

ax[''].yaxis.set_label_position('left') #Write the label of the colorbar to the left

colbar.ax.tick_params(labelsize = 24) #Set the fontsize of the colorbar ticks
colbar.set_label(r'$d_N$ and $d_S$ value', size = 48) #Set the fontsize of the colorbar label
colbar.outline.set_linewidth(4) #Set the width of the colorbar borders

# =============================================================================
# 3. Loop through gene types and generate a heatmap for each type
# =============================================================================
for gene in gene_no.keys():
    #Assign one of the suplots to each gene type
    if gene == 'GS1':
        axn = 'A'
    elif gene == 'GS2':
        axn = 'B'
    elif gene == 'BRS':
        axn = 'C'
    else: axn = 'D'
    
    file = f'{workdir}/results/{gene}/dNdS.tsv' #Get the file with input dN and dS values
    with open(file) as comp_file: #Open input file
        file_info = pd.read_csv(comp_file, sep = '\t') #Read the file as a dataframe
        file_info = file_info.sort_values(by=['locus1', 'locus2']) #Sort pairs by locus tags
        file_info = file_info.reset_index() #Reset the index of the dataframe to match the order
        # file_info = file_info.sort_values('w')
        for key in replace_str.keys(): #Loop through text to replace in the locus tags
            file_info[['locus1', 'locus2']] = file_info[['locus1', 'locus2']].apply(lambda col: col.str.replace(key, replace_str[key])) #Replace said text
        dim = gene_no[gene] #Get the dimennsions of the substitution matrix
        plot_matrix = np.full((dim, dim), np.nan) #Create an empty matrix of those dimensions
        unique_loci = np.unique(file_info[['locus1', 'locus2']].values) #Extract unique loci from the input data
        # mask = np.isin(unique_loci, BRSb, invert = True)
        # end = np.isin(unique_loci, BRSb)
        # unique_loci = np.concatenate([unique_loci[mask], unique_loci[end]]) 
        unique_loci2 = np.array([locus.split('_')[0] if locus not in BRSb else locus.split('_')[0] + 'b' for locus in unique_loci])
        
    for index, row in file_info.iterrows(): #Loop through the dictionary
        l1 = np.where(unique_loci == row['locus1'])[0] #Get position of the matrix of the first locus tag
        l2 = np.where(unique_loci == row['locus2'])[0] #Get position in the matrix of the second locus tag
        plot_matrix[l1, l2] = row['dS'] #Assign the dS value to the right pair
        plot_matrix[l2, l1] = row['dN'] #Assign the dN value to the right pair
    plot_matrix[plot_matrix > vmax] = vmax #If the value > 1.5, convert it to 1.5
    
    for axis in ['top','bottom','left','right']: #Loop through the four sides of the figure frame
        ax[axn].spines[axis].set_linewidth(4) #Set the width of the border

    ax[axn].imshow(plot_matrix, vmin=vmin, vmax = vmax, cmap ='magma') #Plot the data as a heatmap

    ax[axn].set_xticks(np.arange(len(unique_loci))) #Use the unique loci as x axis ticks
    ax[axn].set_xticklabels(unique_loci2, rotation = 90, fontsize = 24) #Format the ticks (vertical position, increase font size)
    ax[axn].set_yticks(np.arange(len(unique_loci))) #Use the unique loci as y axis ticks
    ax[axn].set_yticklabels(unique_loci2, fontsize = 24) #Format the ticks (increase font size)
    ax[axn].set_xlabel('Locus tag 1', size = 36) #Set title of x axis
    ax[axn].set_ylabel('Locus tag 2', size = 36) #Set title of y axis
    ax[axn].set_title(gene, size = 30) #Set title of plot (gene type)
    ax[axn].text(0, 1.1, axn, transform=( #Add plot letter (A-D)
           ax[axn].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)), #Convert position to figure scale
       fontsize = 48, va = 'top', ha = 'left', bbox = dict(facecolor='none',  #Plot to the top left of the plot, add box around
                                                           edgecolor='black', #Add a black edfe to the box
                                                           pad = 20)) #Padding between the text and the box

# =============================================================================
# 4. Save the figure in different formats
# =============================================================================
plt.savefig(outplot, format = 'svg') #Save the figure to SVG format
plt.savefig(outplot.replace('svg', 'png'), format = 'png') #Save the figure to PNG format
plt.savefig(outplot.replace('svg', 'pdf'), format = 'pdf') #Save the figure to PDF format
plt.savefig(outplot.replace('svg', 'tiff'), format = 'tiff') #Save the figure to TIFF format
plt.show() #Show the figure
        
