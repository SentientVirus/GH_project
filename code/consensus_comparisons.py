#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 15:49:44 2025

Code to plot the measures of similarity from Jalview.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================

import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator

# =============================================================================
# 1. Define classes and functions
# =============================================================================

class conserved:
    '''Class to store the consensus information about positions in the alignment'''
    def __init__(self, pos, GS1, GS2, GS3, GS4):
        self.pos = pos # Position
        self.GS1 = GS1 # Amino acid in GS1
        self.GS2 = GS2 # Amino acid in GS2
        self.GS3 = GS3 # Amino acid in GS3
        self.GS4 = GS4 # Amino acid in GS4
        self.GS_list = [GS1, GS2, GS3, GS4] #List with all the amino acids
    def update_GS(self): #Function to update the list (not needed here)
        self.GS_list = [self.GS1, self.GS2, self.GS3, self.GS4]
    def __str__(self): #Function to define how to print the object
        to_print = f'Position {self.pos}:\n'
        for i in range(len(self.GS_list)):
            if self.GS_list[i] != '':
                to_print += f'\t- GS{i+1}, {self.GS_list[i]}\n'
        return to_print

def format_x(x, pos):
    '''This function specifies how to display the tick labels of the x axis'''
    return f'{x+1:.0f}'

# =============================================================================
# 2. Define formatting parameters and inputs
# =============================================================================

# Set up input parameters
datasets = ['GS4', 'GS3', 'GS2', 'GS1', 'GS', 'GS_and_BrS'] #Datasets to analyse

fontsize_title = 64

motif_pos = {1: (693, 700), 2: (187, 197), 3: (223, 235), 4: (310, 323), #Dictionary with motif positions (8 is loop A2)
             5: (144, 157), 6: (625, 637), 7: (380, 389), 8: (637, 661)}

motif_list = [(motif_pos[i][0], motif_pos[i][1] - motif_pos[i][0]) for i in motif_pos.keys()] #Convert motif information to a list

# Paths to input files
indir = os.path.expanduser('~') + '/GH_project/conservation_consensus' # Input directory
consens = [f'{indir}/{dataset}_consensus.csv' for dataset in datasets] # Files with consensus sequences
outplot = f'{indir}/heatmap.png' # Output heatmap
out_tsv = f'{indir}/conserved_positions.tsv' # Output tab file

# Set up the formatting of the plot
colors = ['#B30021', '#CDD0D6', '#FFEA8C'] # Colors to plot the regions that are consistent within types
colors2 = ['#ED9F45', 'black'] # Colors for the regions that are similar within, but not between, subtypes
cmap1 = LinearSegmentedColormap.from_list('red', list(reversed(colors))) #Make a colormap based on the first list
cmap2 = LinearSegmentedColormap.from_list('cyan', list(colors2)) #Make a colormap based on the second list

fig = plt.figure(constrained_layout=False, figsize=(80, 50)) # Create the figure layout
spec = gridspec.GridSpec(nrows = 6, ncols = 1, figure = fig, hspace = 0.8) # Create a grid with six plots (one per dataset)

# Create empty dictionaries and lists to store information
dataset_dict = {}
seq_dict = {}
axs = []
comparison = {}
positions = []

# =============================================================================
# 3. Create the plot of consensus regions
# =============================================================================

for j in range(len(consens)): # Loop through files with consensus sequences
    name = os.path.basename(consens[j]).replace('_consensus.csv', '').replace('_', ' ') #Retrieve the name of the dataset (plot titles)
    with open(consens[j]) as csvfile: # Open input file
        df = pd.read_csv(csvfile, sep = ',', header = None) # Read the file as a dataframe
        to_plot = df.iloc[0] # Retrieve the percentage of identity
        to_plot_list = list(to_plot)[1:-1] # Convert to list and remove empty/header positions
        
        gaps = df.iloc[1] # Retrieve the amino acid sequences
        gap_list = list(gaps)[1:-1] # Convert to list and remove empty/header positions
        sequence = [pos.strip().split(' ')[0] for pos in gap_list] # Remove % identity and clean the formatting
        seq_dict[name] = sequence # Save the complete consensus sequence to a dictionary
        
        # Convert the percentages of identity to floats and keep only positions where identity is 100% and set gap positions to 0% identity
        to_plot_float = [float(to_plot_list[i].strip()) if float(to_plot_list[i].strip()) == 100.00 else np.nan if gap_list[i].strip() != '' else 0.00 for i in range(len(to_plot_list))]

        dataset_dict[name] = to_plot_float # Save the list with this information to a dictionary
        
        to_plot_float = pd.DataFrame(to_plot_float, columns=['Consensus']) # Convert the list to a dataframe

        # Create grids of two plots (one for the motifs and the other for sequence consensus)
        gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = spec[j], hspace = 0.1, height_ratios = [1, 6])
        
        bar_ax = fig.add_subplot(gs[0]) # Use one of the subplots to make a barplot (motifs)
        ax = fig.add_subplot(gs[1]) # Use the other subplot as the main frame (consensus information)
        # ax.set_facecolor('#FEFDED') # Set the background color of the main ax
        
        bar_ax.set_xlim(0, len(to_plot_float)) # Set the x axis to cover the length of the alignment
        bar_ax.xaxis.set_major_locator(MultipleLocator(50)) # Set ticks every 50 positions
        bar_ax.xaxis.set_major_formatter(format_x) # Apply function to display the index without decimals
        bar_ax.broken_barh(motif_list, (0, 0.1), color=['black']*7 + ['slategrey']) # Plot the motifs as bars, with a different color for loop A2
        bar_ax.get_yaxis().set_visible(False) # Hide y axis ticks
        bar_ax.get_xaxis().set_visible(False) # Hide x axis ticks
        bar_ax.set_title(name, fontsize = fontsize_title,
                      fontname = 'Arial') # Add gene names at the top of the plot
        
        htmp = sns.heatmap(to_plot_float.T, cmap = cmap1, center = 5.5, ax = ax,
                          vmin = 0, vmax = 100, linewidths = 0,
                          cbar = False, zorder = 0) # Plot the heatmap in the main ax
        ax.set_facecolor('#FEFDED') # Set the background color of the main ax
        ax.patch.set(lw = 1, ec = 'black') # Border of the subplots
        ax.xaxis.set_major_locator(MultipleLocator(50)) # Set ticks every 50 positions
        ax.xaxis.set_major_formatter(format_x) # Apply function to display the index without decimals

        plt.tick_params(axis='x', which='major', labelsize = 50, # Increase font size of tick labels
                        length = 25, width = 5, rotation = 0) # Make ticks larger and horizontal
        plt.tick_params(axis='y', which='major', labelsize = 60, # Increase font size of tick labels
                        length = 0) # Hide ticks
        axs.append(ax) # Add the ax to a list
    
# =============================================================================
# 4. Overlay regions that are only consistent within subtypes    
# =============================================================================

keys = list(dataset_dict.keys())[:-2] # Get the names of GS subtypes
for k in range(len(keys)): # Loop through the GS subtypes
    g1 = keys[k] # Get first subtype name
    comp = [1 for n in range(len(dataset_dict[g1]))] # Create a list of ones
    for m in range(len(keys)): # Loop through the subbtypes again
        if k != m: # If the two subtypes are different
            g2 = keys[m] # Get the second subtype
            # Convert all the positions where the conservation is not 100% or where there are similarities between subtypes to NaN
            comp = [np.nan if ((dataset_dict[g1][n] != 100) or (dataset_dict[g2][n] != 100) or (seq_dict[g1][n] == seq_dict[g2][n])) else comp[n] for n in range(len(dataset_dict[g1]))] 
    comparison[g1] = comp # Save the comparison to a dictionary
    comp_df = pd.DataFrame(comp, columns = ['Unique']) # Convert the comparison to a dataframe to plot it
    
    new_ax = axs[k].twinx() # Create an ax on top of the main ax
    htmp = sns.heatmap(comp_df.T, cmap = cmap2, center = 5.5, ax = new_ax, # Overlay the positions that are only consistent within subtypes
                      vmin = 0, vmax = 1, linewidths = 0, cbar = False, zorder = 0)
    axs[k].xaxis.set_major_locator(MultipleLocator(50)) # Set ticks every 50 positions
    axs[k].xaxis.set_major_formatter(format_x) # Apply function to display the index without decimals
    new_ax.set_yticklabels([]) # Remove the labels of the y axis
    new_ax.get_yaxis().set_visible(False) # Remove the ticks of the y axis
    
    
# =============================================================================
# 5. Create a tab file with positions that are conserved within subtypes
# =============================================================================

# Print positions that are unique in each subtype
conserved_all = [print(f'Position {n} is conserved in all GS subtypes!') for n in range(len(comparison['GS1'])) if comparison['GS1'][n] == True == comparison['GS2'][n] == comparison['GS3'][n] == comparison['GS4'][n]]

with open(out_tsv, 'w') as tsvfile: # Open output file
    tsvfile.write('Alignment_position\tGS1_aa\tGS2_aa\tGS3_aa\tGS4_aa\n') # Write headers
    for n in range(len(comparison['GS1'])): # Loop through the length of the alignment
        if any([comparison[x][n] == True for x in comparison.keys()]): # If any conserved position differs between subtypes
                cons = conserved(n+1, seq_dict['GS1'][n], seq_dict['GS2'][n], seq_dict['GS3'][n], seq_dict['GS4'][n]) # Store the information in an object
                print(cons) # Print the object
                tsvfile.write(f'{cons.pos}\t{cons.GS1}\t{cons.GS2}\t{cons.GS3}\t{cons.GS4}\n') # Write the information to file
    
# =============================================================================
# 6. Save the plot
# =============================================================================
    
# Save the plot in multiple formats
plt.savefig(outplot, bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.pdf'), bbox_inches='tight')
