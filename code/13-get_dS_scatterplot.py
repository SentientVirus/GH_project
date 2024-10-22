#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 14:24:33 2024

@author: marina
"""
import os
import pandas as pd
from matplotlib import pyplot as plt

# =============================================================================
# 1. Define the paths to input and output files
# =============================================================================
workdir = os.path.expanduser('~') + '/GH_project' #Working directory
pairwise_ds = f'{workdir}/results' #Directory with input pairwise dS data
core_ds = f'{pairwise_ds}/core_genes/core_pairwise_metrics.tsv' #File with core dS data
outplot = f'{workdir}/plots/scatter/GH_types_scatterplot.png' #Output file

GH_types = ['GS1', 'GS2', 'BRS', 'NGB', 'S2a', 'S3'] #Gene types to be plotted

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
with open(core_ds) as core: #Open the file with core dS values
    df = pd.read_csv(core, sep = '\t') #Load the file as a dataframe
    for index, row in df.iterrows(): #Loop through the rows in the dataframe
        strain1 = replace_strain_name(row['strain1']) #Retrieve the first strain in the pair
        strain2 = replace_strain_name(row['strain2']) #Retrieve the second strain in the pair
        if not strain1[0] == 'L' and not strain2[0] == 'L': #Skip strains with names starting with L
            core_ds_dict[tuple(sorted((strain1, strain2)))] = row['mean_dS'] #Use sorted strain names as key and core pairwise dS as value
        
#Create a plot with two rows and six columns, where all columns in a row share y axis
fig, axs = plt.subplots(2, 6, constrained_layout = False, sharey = 'row', 
                        figsize=(18, 6))

axs[0, 0].set_ylabel('Core ${d_S}$') #Add a title to the y axis of the first scatter plot
axs[1, 0].set_ylabel('Count') #Add a title to the y axis of the first histogram 
fig.text(0.5, 0.47, 'Pairwise ${d_S}$', ha = 'center') #Add title of the x axis of all scatter plots
fig.text(0.5, 0.03, 'Difference to core ${d_S}$', ha = 'center') #Add title of the x axis of all histograms
plt.subplots_adjust(hspace = 0.3, wspace = 0.2) #Adjust horizontal and vertical padding between plots

for i in range(len(GH_types)): #Loop through gene types
    GH = GH_types[i] #Name of the gene type
    axs[0, i].set_title(GH) #Use the name of the gene as the title of the scatter plots
    
    label_list = [] #Create an empty list to store labels (not used)
    x_list = [] #Create a list to store pairwise dS values
    y_list = [] #Create a list to store core pairwise dS values
    dif_list = [] #Create a list to store the difference between x and y
    
    with open(f'{pairwise_ds}/{GH}/dNdS.tsv') as gene: #Open file with dS values
        df2 = pd.read_csv(gene, sep = '\t') #Load file as a dataframe
        for index2, row2 in df2.iterrows(): #Loop through rows in the dataframe
            label = f'{row2["locus1"]}-{row2["locus2"]}' #Create a label with the locus tags in the comparison
            x = row2['dS'] #Get pairwise dS
            strain1 = replace_strain_name(row2['locus1']) #Get the name of the first strain
            strain2 = replace_strain_name(row2['locus2']) #Get the name of the second strain
            if strain1 == strain2: #If both genes are in the same strain
                y = 0 #Core dS is 0
            else: y = core_ds_dict[tuple(sorted((strain1, strain2)))] #Otherwise, store core dS in a dictionary
            
            #Set all dS values above 1.5 to 1.5
            if x > 1.5:
                x = 1.5
            if y > 1.5:
                y = 1.5
                
            label_list.append(label) #Add label to the label list
            x_list.append(x) #Add x to the pairwise dS list
            y_list.append(y) #Add y to the core pairwise dS list
            dif_list.append(x - y) #Store the difference between x and y in a list
            
    #Retrieve values < core dS and values > core dS separately
    subx = [x_list[j] for j in range(len(y_list)) if x_list[j] <= y_list[j]]
    suby = [y_list[j] for j in range(len(y_list)) if x_list[j] <= y_list[j]]
    overx = [x_list[j] for j in range(len(y_list)) if x_list[j] > y_list[j]]
    overy = [y_list[j] for j in range(len(y_list)) if x_list[j] > y_list[j]]
    #Plot values < core dS and values > core dS separately to color them differently
    axs[0, i].scatter(subx, suby)
    axs[0, i].scatter(overx, overy)
    #Set the x axes limits of the scatter plots to 0-1.5 and the y axis lower limit to 0
    axs[0, i].set_xlim(0, 1.5)
    axs[0, i].set_ylim(0)
    #Plot a diagonal line where y = x
    axs[0, i].plot([0, 1], [0, 1], color = 'black')
    #Plot the differences between x and y separately depending of if they are > 0
    sub_dif = [dif for dif in dif_list if dif <= 0]
    over_dif = [dif for dif in dif_list if dif > 0]
    #Set the limit of the x axis of the histogram plots
    axs[1, i].set_xlim(-0.5, 1.5)
    #Plot values >0 and <0 separately, set number of bars in the histogram
    axs[1, i].hist(sorted(sub_dif), bins = int((max(sub_dif)-min(sub_dif))//0.0625))
    axs[1, i].hist(sorted(over_dif), bins = int((max(over_dif)-min(over_dif))//0.0625))
    #Add a vertival line over 0´
    axs[1, i].axvline(0, color = 'black')
    
plt.savefig(outplot, bbox_inches = 'tight') #Save the image as PNG
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches = 'tight') #Save the image as TIFF
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches = 'tight') #Save the image as SVG
    