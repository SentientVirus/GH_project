#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates a plot of the phylogeny of A. kunkeei strains, with
a representation of the GH70 and GH32 genes that are present in a particular
region of the genome. The phylogenies used include representative strains only, 
and only the GS1, GS2 and BRS phylogenies were included.

@author: Marina Mota Merlo
"""

from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace
#Use PhyloTree
import csv, os, logging, traceback
import yaml
import numpy as np
from collections import Counter

# =============================================================================
# Logging
# =============================================================================

# logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
#                     format = '%(asctime)s %(message)s',
#                     datefmt = '%Y-%m-%d %H:%M:%S')

# def handle_exception(exc_type, exc_value, exc_traceback):
#     if issubclass(exc_type, KeyboardInterrupt):
#         sys.__excepthook__(exc_type, exc_value, exc_traceback)
#         return

#     logger.error(''.join(["Uncaught exception: ",
#                           *traceback.format_exception(exc_type, exc_value, exc_traceback)
#                           ]))

# sys.excepthook = handle_exception

# sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# Set colors for strains and CDS
# =============================================================================

#Dictionary assigning colors to each strain
leaf_color = {'A0901': '#D55E00', 'A1001': '#771853', 'A1002': '#D55E00',
              'A1003': '#0072B2', 'A1201': '#D55E00', 'A1202': '#33B18F',
              'A1401': '#33B18F', 'A1404': '#FF74D6', 'A1802': '#D55E00',
              'A1803': '#D55E00', 'A1805': '#33B18F', 'A2001': '#D55E00',
              'A2002': '#771853', 'A2003': '#771853', 'A2101': '#33B18F',
              'A2102': '#D55E00', 'A2103': '#0072B2', 'G0101': '#0072B2', 
              'G0102': '#0072B2', 'G0103': '#0072B2', 'G0401': '#0072B2',
              'G0402': '#0072B2', 'G0403': '#33B18F', 'G0404': '#0072B2',
              'G0405': '#0072B2', 'G0406': '#33B18F', 'G0407': '#0072B2',
              'G0408': '#0072B2', 'G0410': '#0072B2', 'G0412': '#0072B2',
              'G0414': '#0072B2', 'G0415': '#0072B2', 'G0417': '#0072B2',
              'G0420': '#33B18F', 'G0601': '#0072B2', 'G0602': '#0072B2',
              'G0702': '#0072B2', 'G0801': '#0072B2', 'G0802': '#0072B2',
              'G0803': '#0072B2', 'G0804': '#0072B2', 'Fhon2': '#0072B2',
              'H1B104J': '#0072B2', 'H1B105A': '#0072B2', 'H1B302M': '#0072B2',
              'H2B105J': '#D55E00', 'H3B101A': '#0072B2', 'H3B101J': '#0072B2',
              'H3B101X': '#D55E00', 'H3B102A': '#0072B2', 'H3B102X': '#D55E00',
              'H3B103J': '#0072B2', 'H3B103X': '#D55E00', 'H3B103M': '#0072B2', 
              'H3B104J': '#0072B2', 'H3B104X': '#0072B2', 'H3B107A': '#0072B2', 
              'H3B109M': '#0072B2', 'H3B110M': '#0072B2', 'H3B111A': '#0072B2', 
              'H3B111M': '#0072B2', 'H3B202M': '#D55E00', 'H3B202X': '#0072B2', 
              'H3B203J': '#0072B2', 'H3B203M': '#D55E00', 'H3B204J': '#0072B2', 
              'H3B204M': '#D55E00', 'H3B205J': '#0072B2', 'H3B206M': '#D55E00', 
              'H3B207X': '#0072B2', 'H3B208X': '#0072B2', 'H3B209X': '#33B18F', 
              'H4B101A': '#D55E00', 'H4B102A': '#D55E00', 'H4B103J': '#D55E00', 
              'H4B104A': '#D55E00', 'H4B111J': '#D55E00', 'H4B114J': '#D55E00', 
              'H4B116J': '#D55E00', 'H4B202J': '#0072B2', 'H4B203M': '#D55E00', 
              'H4B204J': '#0072B2', 'H4B205J': '#0072B2', 'H4B206J': '#D55E00', 
              'H4B210M': '#D55E00', 'H4B211M': '#0072B2', 'H4B303J': '#D55E00', 
              'H4B402J': '#0072B2', 'H4B403J': '#D55E00', 'H4B404J': '#D55E00', 
              'H4B405J': '#D55E00', 'H4B406M': '#D55E00', 'H4B410M': '#D55E00', 
              'H4B411M': '#D55E00', 'H4B412M': '#0072B2', 'H4B501J': '#0072B2', 
              'H4B502X': '#0072B2', 'H4B503X': '#0072B2', 'H4B504J': '#33B18F',
              'H4B505J': '#33B18F', 'H4B507J': '#0072B2', 'H4B507X': '#0072B2', 
              'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': '#D55E00', 
              'DSM': '#0072B2'}

# =============================================================================
# Set target region and load input files
# =============================================================================
# TO DO: Integrate into Snakemake
outdir =  os.path.expanduser('~') + '/GH_project/plots/trees' #Directory to store output plots
GH_types = ['BRS', 'GS1', 'GS2', 'NGB', 'S2a'] #Types of GH to plot

#Create output directory if it doesn't exist
if not os.path.exists(outdir): 
    os.makedirs(outdir)

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
        
    locus_tag = locus_tag.replace('LDX55', 'IBH001').replace('APS55', 'MP2').replace('K2W83', 'DSM').replace('RS', '').replace('FHON', 'Fhon').replace('AKU', '')
    return locus_tag

def remove_minus(strains):
    """"Simple function to remove the - symbol from any string in a list of
    strings."""
    new_list = [strain.replace('-', '') for strain in strains]
    return new_list

def fix_strand(my_info_list):
    """My genes of interest are in the reverse strand, but the locus tags are
    ordered based on the forward strand. This function reverts the order of
    the genes before plotting.
    
    Parameters
    ----------
    my_info_list : 
        Dictionary of tuples, where the keys are the gene names and the values
        are the start and end positions."""
        
    for el in reversed(my_info_list): #Loop through elements in list
        my_index = my_info_list.index(el) #Get the index of the element
        if my_index == len(my_info_list)-1: #If the element is at the back of the list
            start = padding + gapscale #Set the start at the beginning
        else: #Otherwise
            start = pass_next - el[1] + end  #Set the position as the end of the previous element plus the distance from this element to the previous one
        end = start + (el[1] - el[0]) #Set the end position of the current element as the start plus the length of the element
        pass_next = el[0] #Start of the previous element of the list
        my_info_list[my_index][0] = start #Save the start position to a dictionary
        my_info_list[my_index][1] = end #Save the end position to a dictionary
    return my_info_list


for GH in GH_types:
    treefile =  os.path.expanduser('~') + f'/GH_project/data/fasta/GH70/trees/{GH}_repset.mafft.faa.treefile' #Load tree file
    if GH == 'S2a':
        treefile = treefile.replace('70', '32')
    outplot = f'{GH}_phylogeny.png' #Specify name of the output plot
    
# =============================================================================
# Create rooted tree from treefile
# =============================================================================
    t = Tree(treefile, format = 0) #Read tree file
    print(GH)
    #Root trees
    if GH == 'GS1':
        outnode = 'A1001_12310' #Root of GS1
    elif GH == 'GS2':
        outnode = t.get_common_ancestor('H3B104X_13220', 'H4B505J_12900') #Root of GS2
    elif GH == 'BRS':
        outnode = t.get_common_ancestor('LDX55_06330', 'H4B204J_13340') #Root of BRS
    elif GH == 'NGB':
        outnode = 'A0901_05380' #Root of NGB
    else: outnode = 'A1001_12300' #Root of S2a
    t.set_outgroup(outnode) #Set the root
    
# =============================================================================
# Set style of tree
# =============================================================================
    ts = TreeStyle() #Create default tree style
    ts.show_branch_length = False # Hide support values
    ts.scale = 2000 #General tree scale
    ts.scale_length = 0.05
    ts.branch_vertical_margin = -5 # Space between branches
    ts.show_branch_support = False
    ts.show_leaf_name = False # Hide unformatted leaf names
    ts.show_scale = True # Show tree scale
    
    #Add legend
    ts.legend_position = 3 #Bottom left
    ts.legend.add_face(TextFace(' '*58, fsize = 10), column = 0) #Legend formatting
    
    #Formatting of nodes (and branch lines)
    ns = NodeStyle()
    ns['size'] = 0 #Nodes are not visible
    ns['vt_line_width'] = 2 #Width of vertical lines
    ns['hz_line_width'] = 2 #Width of horizontal lines
    ns['hz_line_type'] = 0 #Solid horizontal lines
    for n in t.traverse(): #Loop through the nodes in the tree
        n.set_style(ns) #Set the node style
        if n not in t.get_leaves() and n.support > 1: #If there's any node support
            if n.support >= 95: #If the support is over 95%
                color = 'black' #Color the support value in black
            else: color = 'dimgrey' #Otherwise, color the support value in grey
            if n.support >= 80: #If the support value is above 80%
                support_face = TextFace(int(n.support), fgcolor = color, fsize = 12) #Format support value text
                n.add_face(support_face, column=0, position='branch-top') #Add the text to the tree
    
# =============================================================================
# Save CDS information for the main plot
# =============================================================================
    #Get strain names
    leaves = t.get_leaves() #Get tree leaves sorted by phylogeny
    lnames = [replace_strain_name(leaf.name) for leaf in leaves] #Get the strain names to be plotted 
    lnames.reverse() #Reverse list with leaf names

# =============================================================================
# Create the plot
# =============================================================================
    for leaf in leaves: #Loop through leaves in the tree
        lname = replace_strain_name(leaf.name) #Get the right name for each leaf
        strain = lname.split('_')[0] #Get the name of the strain
        color = leaf_color.get(strain, None) #Color leaves according to strain phylogroup
        name_face = TextFace(lname, fgcolor = color, fsize = 20) #Create text with corrected locus tag
        leaf.add_face(name_face, column=0, position='branch-right') #Add formatted leaf names
        
# =============================================================================
# Print tree
# =============================================================================
    t.ladderize(1)
    t.render(f'{outdir}/{outplot}', tree_style = ts) #Render tree in PNG format
    t.render(f'{outdir}/{outplot.replace(".png", ".tiff")}', tree_style = ts) #Render tree in TIFF format
    t.render(f'{outdir}/{outplot.replace(".png", ".svg")}', tree_style = ts) #Render tree in SVG format
