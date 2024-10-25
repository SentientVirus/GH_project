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

# =============================================================================
# 1. Define the paths to input and output files
# =============================================================================
workdir = os.path.expanduser('~') + '/GH_project' #Working directory
pairwise_ds = f'{workdir}/results' #Directory with input pairwise dS data
core_ds = f'{pairwise_ds}/core_genes/core_pairwise_metrics.tsv' #File with core dS data
outplot = f'{workdir}/plots/scatter/GH_types_scatterplot.png' #Output file
coredir = f'{workdir}/results/core_genes/GH_repr'

GH_types = ['GS1', 'GS2', 'BRS', 'NGB', 'S2a', 'S3'] #Gene types to be plotted
to_exclude = ['A1001', 'A1404']

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
        
#Create a plot with two rows and six columns, where all columns in a row share y axis
fig, axs = plt.subplots(6, 6, constrained_layout = False, sharey = 'row', 
                        figsize=(18, 18))

axs[0, 0].set_ylabel('Core ${d_S}$', fontsize = 12) #Add a title to the y axis of the dS scatter plot
axs[1, 0].set_ylabel('Count', fontsize = 12) #Add a title to the y axis of the dS histogram
axs[2, 0].set_ylabel('Core ${d_N}$', fontsize = 12) #Add a title to the y axis of the dN scatter plot
axs[3, 0].set_ylabel('Count', fontsize = 12) #Add a title to the y axis of the dN histogram
axs[4, 0].set_ylabel('${d_N}$', fontsize = 12) #Add a title to the y axis of the dN vs dS scatter plot
axs[5, 0].set_ylabel('Core ${d_N}$') #Add a title to the y axis of the dN vs dS scatter plot
fig.text(0.5, 0.754, 'Pairwise ${d_S}$', ha = 'center', fontsize = 12) #Add title of the x axis of dS scatter plots
fig.text(0.5, 0.623, 'Difference to core ${d_S}$', ha = 'center', fontsize = 12) #Add title of the x axis of the dS histograms
fig.text(0.5, 0.4925, 'Pairwise ${d_N}$', ha = 'center', fontsize = 12) #Add title of the x axis of dN scatter plots
fig.text(0.5, 0.362, 'Difference to core ${d_N}$', ha = 'center', fontsize = 12) #Add title of the x axis of the dN histograms
fig.text(0.5, 0.231, '${d_S}$', ha = 'center', fontsize = 12) #Add title of the dN vs dS scatterplot
fig.text(0.5, 0.1, 'Core ${d_S}$', ha = 'center', fontsize = 12) #Add title of the core dN vs dS scatterplot
plt.subplots_adjust(hspace = 0.3, wspace = 0.2) #Adjust horizontal and vertical padding between plots

for i in range(len(GH_types)): #Loop through gene types
    GH = GH_types[i] #Name of the gene type
    axs[0, i].set_title(GH, fontsize = 14) #Use the name of the gene as the title of the scatter plots
    
    label_list = [] #Create an empty list to store labels (not used)
    x_list = [] #Create a list to store pairwise dS values
    y_list = [] #Create a list to store core pairwise dS values
    dif_list = [] #Create a list to store the difference between x and y
    dN_list = []
    dN_y_list = []
    dN_dif_list = []
    
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
            
            #Set all dS values above 1.5 to 1.5
            if x > 1.5:
                x = 1.5
            if y > 1.5:
                y = 1.5
                
            label_list.append(label) #Add label to the label list
            x_list.append(x) #Add x to the pairwise dS list
            y_list.append(y) #Add y to the core pairwise dS list
            dif_list.append(x - y) #Store the difference between x and y in a list
            dN_list.append(dN)
            dN_y_list.append(dN_y)
            dN_dif_list.append(dN - dN_y) #Store the difference between dN (pairwise) and dN (core) in a list
            
    #Retrieve values < core dS and values > core dS separately
    subx = [x_list[j] for j in range(len(y_list)) if x_list[j] <= y_list[j]]
    suby = [y_list[j] for j in range(len(y_list)) if x_list[j] <= y_list[j]]
    overx = [x_list[j] for j in range(len(y_list)) if x_list[j] > y_list[j]]
    overy = [y_list[j] for j in range(len(y_list)) if x_list[j] > y_list[j]]
    
    sub_dN = [dN_list[j] for j in range(len(dN_y_list)) if dN_list[j] <= dN_y_list[j]]
    sub_dN_y = [dN_y_list[j] for j in range(len(dN_y_list)) if dN_list[j] <= dN_y_list[j]]
    over_dN = [dN_list[j] for j in range(len(dN_y_list)) if dN_list[j] > dN_y_list[j]]
    over_dN_y = [dN_y_list[j] for j in range(len(dN_y_list)) if dN_list[j] > dN_y_list[j]]
    
    sub_x = [x_list[j] for j in range(len(dN_list)) if x_list[j] <= dN_list[j]]
    sub_y = [dN_list[j] for j in range(len(dN_list)) if x_list[j] <= dN_list[j]]
    over_x = [x_list[j] for j in range(len(dN_list)) if x_list[j] > dN_list[j]]
    over_y = [dN_list[j] for j in range(len(dN_list)) if x_list[j] > dN_list[j]]
    
    #Plot values < core dS and values > core dS separately to color them differently
    axs[0, i].scatter(subx, suby, alpha = 0.2)
    axs[0, i].scatter(overx, overy, alpha = 0.2)
    axs[0, i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    axs[1, i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    axs[2, i].scatter(sub_dN, sub_dN_y, alpha = 0.2)
    axs[2, i].scatter(over_dN, over_dN_y, alpha = 0.2)
    axs[2, i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    axs[3, i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    axs[4, i].scatter(sub_x, sub_y, alpha = 0.2)
    axs[4, i].scatter(over_x, over_y, alpha = 0.2)
    axs[4, i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    #Set the x axes limits of the scatter plots to 0-1.5 and the y axis lower limit to 0
    axs[0, i].set_xlim(0, 1.5)
    axs[0, i].set_ylim(0)
    
    axs[4, i].set_xlim(0, 1.5)
    axs[4, i].set_ylim(0, 0.45)
    
    if GH not in ['GS1', 'GS2', 'BRS']:
        axs[2, i].set_xlim(0, 0.075)
        # axs[2, i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    else: axs[2, i].set_xlim(0, 0.25)
    axs[2, i].set_ylim(0)
    
    #Plot a diagonal line where y = x
    axs[0, i].plot([0, 1], [0, 1], color = 'black')
    axs[2, i].plot([0, 1], [0, 1], color = 'black')
    axs[4, i].plot([0, 1], [0, 1], color = 'black')
    # axs[4, i].plot([0, 1], [0, 1], color = 'black')
    #Plot the differences between x and y separately depending of if they are > 0
    sub_dif = [dif for dif in dif_list if dif <= 0]
    over_dif = [dif for dif in dif_list if dif > 0]
    
    sub_dN_dif = [dif for dif in dN_dif_list if dif <= 0]
    over_dN_dif = [dif for dif in dN_dif_list if dif > 0]
    
    #Set the limit of the x axis of the histogram plots
    axs[1, i].set_xlim(-0.3, 1.5)
    axs[3, i].set_xlim(-0.01, 0.24)
    #Plot values >0 and <0 separately, set number of bars in the histogram
    axs[1, i].hist(sorted(sub_dif), bins = int((max(sub_dif)-min(sub_dif))//0.0625))
    axs[1, i].hist(sorted(over_dif), bins = int((max(over_dif)-min(over_dif))//0.0625))
    axs[3, i].hist(sorted(sub_dN_dif), bins = 1)
    axs[3, i].hist(sorted(over_dN_dif), bins = int(((max(over_dN_dif)-min(over_dN_dif))*8)//0.0625))
    #Add a vertical line over 0
    axs[1, i].axvline(0, color = 'black')
    axs[3, i].axvline(0, color = 'black')
    
    dS_corelist = []
    dN_corelist = []
    for file in os.listdir(coredir):
        if file.startswith(GH) and file.endswith('dNdS.tsv'):
            df = pd.read_csv(f'{coredir}/{file}', sep = '\t')
            for index, row in df.iterrows():
                check = any(strain in row['locus1'] for strain in to_exclude) or any(strain in row['locus2'] for strain in to_exclude)
                if not check and not (GH == 'GS1' and ('A1003' in row['locus1'] or 'A1003' in row['locus2'])):
                    dS = min(1.5, row['dS'])
                    dN = min(1.5, row['dN'])
                    dS_corelist.append(dS)
                    dN_corelist.append(dN)
                    
    axs[5, i].set_xlim(0, 1.5)
    axs[5, i].set_ylim(0, 0.45)
    
    sub_xcore = [dS_corelist[j] for j in range(len(dN_corelist)) if dS_corelist[j] <= dN_corelist[j]]
    sub_ycore = [dN_corelist[j] for j in range(len(dN_corelist)) if dS_corelist[j] <= dN_corelist[j]]
    over_xcore = [dS_corelist[j] for j in range(len(dN_corelist)) if dS_corelist[j] > dN_corelist[j]]
    over_ycore = [dN_corelist[j] for j in range(len(dN_corelist)) if dS_corelist[j] > dN_corelist[j]]
    axs[5, i].scatter(sub_xcore, sub_ycore, alpha = 0.2, s = 15)
    axs[5, i].scatter(over_xcore, over_ycore, alpha = 0.2, s = 15)
    axs[5, i].plot([0, 1], [0, 1], color = 'black')
    axs[5, i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    
# Add global legend
sub_patch = Line2D([0], [0], marker = 'o', color = '#1F77B4', label = 'x ≤ y',
                        markerfacecolor = '#1F77B4', markersize = 10, 
                        alpha = 0.2, linewidth = 0)

rect_sub = mpatches.Patch(color = '#1F77B4', label = 'x - y ≤ 0', ec = 'white', 
                          lw = 0)

over_patch = Line2D([0], [0], marker = 'o', color = '#FF7F0E', label = 'x > y',
                        markerfacecolor = '#FF7F0E', markersize = 10, 
                        alpha = 0.2, linewidth = 0, markeredgewidth = 2)

line_patch = Line2D([0, 1], [1, 0], marker = '|', color = 'black', label = 'x = y',
                    markersize = 10, alpha = 1, linewidth = 0)

rect_over = mpatches.Patch(color = '#FF7F0E', label = 'x - y > 0', ec = 'white', 
                           lw = 0)

patches = [sub_patch, over_patch, line_patch, rect_sub, rect_over]

legend = plt.legend(handles = patches,  handlelength = 0.6, handleheight = 1.5, 
                    title = 'Legend', loc = 'upper right', framealpha = 0, 
                    frameon = False, fontsize = 12, bbox_to_anchor = (1.7, 7.71))
           
plt.setp(legend.get_title(),fontsize = 14)
legend._legend_box.align = 'left'
    
plt.savefig(outplot, bbox_inches = 'tight') #Save the image as PNG
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches = 'tight') #Save the image as TIFF
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches = 'tight') #Save the image as SVG
    