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
              'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': 'black', 
              'DSMZ12361': 'black'}

#Highlighted gene colors for GS1-2 and BRS, faint colors for the rest
# gene_colors = {'GS1': '#FF7575', 'GS2': '#FF75B6', 'BRS': '#E875FF',
#                'S1': '#FFA175', 'S2': '#FFB875', 'S2a': '#FFB875', 'S2b': '#FFD775', 
#                'S3': '#FFF575', 'mhpD': '#CCCCCC', 'hpcG': '#CCCCCC', 
#                'oppA': '#CCCCCC', 'GS3': '#FFC2D0', 'GS4': '#FF8A75'}

# #Faint colors for GS1-2 and BRS
# alt_colors = {'GS1': '#FFC0C0', 'GS2': '#FFC0DE', 'BRS': '#F4C0FF', 
#               'S1': '#FFD3BD', 'S2': '#FFE0C3', 'S2a': '#FFE0C3', 'S2b': '#FFF0CB', 
#               'S3': '#FFFCD0'}

gene_colors = {'GS1': '#FF7575', 'GS2': '#ca7afe', 'BRS': '#7a84fe',
                'S1': '#FFF575', 'S2': '#FFB875', 'S2a': '#b5ff70', 'S2b': '#74ff85', 
                'S3': '#71f6ff', 'mhpD': '#CCCCCC', 'hpcG': '#CCCCCC', 
                'oppA': '#CCCCCC', 'nox': '#CCCCCC', 'GS3': '#ffc0f9', 'GS4': '#FF8A75'}

#Faint colors for GS1-2 and BRS
alt_colors = {'GS1': '#FFC0C0', 'GS2': '#e6c0ff', 'BRS': '#b5bbff', 
              'S1': '#fcffc2', 'S2': '#d7feb2', 'S2a': '#d7feb2', 'S2b': '#b8ffc0', 
              'S3': '#bafbff'}

# =============================================================================
# Set target region and load input files
# =============================================================================
# TO DO: Annotate script and integrate into Snakemake

segment_length = 29000 #25000 #Length of the graphical representation of the CDS in the plot
gapscale = 1000 #Gap added for padding at the beginning
padding = 500 #Padding on the horizontal line with the CDS
config_file = os.path.expanduser('~') + '/GH_project/config.yml' #File with gene groups
indir =  os.path.expanduser('~') + '/GH_project/plots/tabfiles' #Directory with input tabs
outdir =  os.path.expanduser('~') + '/GH_project/plots/trees' #Directory to store output plots
GH70_types = ['BRS', 'GS1', 'GS2'] #Types of GH to plot
GH32_types = ['S1', 'S2', 'S3'] #Types of GH to plot

domain_path =  os.path.expanduser('~') + '/interproscan' #Path to domain annotations

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

# =============================================================================
# Read types from Snakemake
# =============================================================================
# Read locus tags of different gene subsets
with open(config_file) as conf: #Open config file with dataset informatin
    py_config = yaml.safe_load(conf) #Load the information in the file
    GS1_repr = py_config['GS1_repr'] #Retrieve the representative strains with GS1
    GS2_repr = py_config['GS2_repr'] #Retrieve the representative strains with GS2
    BRS_repr = py_config['BRS_repr'] #Retrieve the representative strains with BRS
    GH70_repr = list(np.unique(GS1_repr + GS2_repr + BRS_repr)) #Create a list with them
    
    S1_repr = remove_minus(py_config['S1_repr']) #Retrieve the representative strains with S1
    S2a_repr = remove_minus(py_config['S2a_repr']) #Retrieve the representative strains with S2a
    S2b_repr = remove_minus(py_config['S2b_repr']) #Retrieve the representative strains with S2b
    S3_repr = remove_minus(py_config['S3_repr']) #Retrieve the representative strains with S3
    
    
    GS1 = py_config['GS1'] #Retrieve the GS1 locus tags
    GS2 = py_config['GS2'] #Retrieve the GS2 locus tags
    GS3 = py_config['GS3'] #Retrieve the GS3 locus tags
    GS3_repr = [tag.split('_')[0] for tag in GS3] #Retrieve the GS3 representative strains
    S1 = py_config['S1'] #Retrieve the S1 locus tags
    S2a = py_config['S2a'] #Retrieve the S2a locus tags
    S2b = py_config['S2b'] #Retrieve the S2b locus tags
    S3 = py_config['S3'] #Retrieve the S3 locus tags
    BRS = ['A1401_12770'] + py_config['BRS'] #Retrieve the BRS locus tags
    GH70_doms = GS1 + GS2 + BRS #Put the GS1-2 and BRS locus tags together
    GH32_doms = S1 + S2a + S2b + S3
    GH_doms = GH70_doms + GH32_doms
    strains = py_config['representatives'] #Get all representative strains
    s_nominus = remove_minus(strains) #Remove minus symbol from strain names
    GS2_BRS =  [GH_gene.replace('_2', '') for GH_gene in GS2 if GH_gene.replace('_2', '_1') in BRS] # Retrieve BRS domains of the GS2_BRS proteins

# =============================================================================
# Load tree file and save gene types to dictionary
# =============================================================================
domain_dict = {} #Create empty dictionary
for GH in GH70_types + GH32_types: #Loop through gene types
    if GH in GH70_types:
        print(f'\n\nGH70 {GH} genes...', end = '\n') #Print gene type
        treefile =  os.path.expanduser('~') + f'/GH_project/data/fasta/GH70/trees/{GH}_repset.mafft.faa.treefile' #Load tree file
    else:
        print(f'\n\nGH32 {GH} genes...', end = '\n') #Print gene type
        treefile =  os.path.expanduser('~') + f'/GH_project/data/fasta/GH32/trees/{GH}_repset.mafft.faa.treefile' #Load tree file
    outplot = f'{GH}_CDS_phylogeny.png' #Specify name of the output plot

    gene_types = {} # Create dictionary to assign gene types
    for gene_dom in GH_doms: #Loop through locus tags
        gene = replace_strain_name(gene_dom) #Change locus tag
        if gene_dom in GS1: #If the locus tag is in GS1
            gene_types[gene] = 'GS1' #Assign it to type GS1
        elif gene_dom in GS2 or (gene_dom[:-2] + '_2' in GS2): #If the locus tag is in GS2 (or part of a GS2_BRS)
            if gene_dom[:-2] + '_1' in BRS: #If the locus tag is also in BRS
                gene_types[gene] = 'GS2_BRS' #Assign the gene to GS2_BRS
            else: gene_types[gene] = 'GS2' #Otherwise, assign the gene to GS2
        elif gene_dom in BRS: #If the locus tag is in BRS
            gene_types[gene] = 'BRS' #Assign the gene to BRS
        elif gene_dom in GS3: #If the locus tag is in GS3
            gene_types[gene] = 'GS3' #Assign the gene to GS3
        elif gene_dom in S1:
            gene_types[gene] = 'S1'
        elif gene_dom in S2a+S2b:
            gene_types[gene] = 'S2'
        elif gene_dom in S3:
            gene_types[gene] = 'S3'
        else: print('Something is VERY wrong here...') #Otherwise, print an error message
            
    tabs = [file for file in os.listdir(indir) if any(list(strain in file for strain in strains)) and file.endswith('.tab')] #snakemake.input.tabs #Maybe I should regenerate the tabs adding locus tag and annotation separately...

# =============================================================================
# Read Interproscan to get the start and end positions of domains within each gene   
# =============================================================================
    domains = [file for file in os.listdir(domain_path) if any(list(strain.replace('Fhon2', 'fhon2') in file for strain in strains))] #Create a list of tab files
    for domain in domains: #Loop through tab files
        domain_file = f'{domain_path}/{domain}' #Get the full path to the file
        with open(domain_file, 'r') as dfile: #Open the file
            freader = csv.reader(dfile, delimiter = '\t') #Read the file contents
            for line in freader: #Loop through the lines in the file
                gene_locus = replace_strain_name(line[0]).replace('-', '').replace('fhon2', 'Fhon2') #Get the locus tag in each line
                if ((line[0] in GH_doms) or (gene_locus.upper() in GH_doms) or (f'{gene_locus.upper()}_1' in GH_doms) or (f'{line[0]}_1' in GH_doms) or ('MP2' in domain_file and int(gene_locus.split('_')[1]) < 14000)) and ('Glycosyl hydrolase family 70' in line[5] or 'glyco_32' in line[5]): #Get the locus tags of interest
                    if 'MP2' in domain_file: #If the strain is MP2
                        gene_locus = gene_locus.replace('_13350', '_03850').replace('_13360', '_03845') #Make extra changes to the locus tags
                    pos = (int(line[6])*3, int(line[7])*3) #Get the position of the domains
                    if gene_locus in domain_dict.keys() and domain_dict[gene_locus] != pos: #If the locus tag is already stored
                        if 'DSM' not in domain_file: #And DSM is not in the file name
                            domain_dict[f'{gene_locus}_1'] = domain_dict[gene_locus] #The position of the BRS domain is the previously stored one
                            del domain_dict[gene_locus] #Remove locus tag from dictionary keys
                            gene_locus = f'{gene_locus}_2' #The new position is the position of GS2
                        else: #If DSM is in the file name
                            domain_dict[f'{gene_locus}_2'] = domain_dict[gene_locus] #The previously stored position corresponds to GS2
                            del domain_dict[gene_locus] #Remove locus tag key
                            gene_locus = f'{gene_locus}_1' #The newly stored position corresponds to BRS
                    domain_dict[gene_locus] = pos #Store the position in the dictionary
    
# =============================================================================
# Create rooted tree from treefile
# =============================================================================
    t = Tree(treefile, format = 0) #Read tree file
    
    #Root trees
    if GH == 'GS1':
        outnode = 'A1001_12310' #Root of GS1
    elif GH == 'GS2':
        outnode = t.get_common_ancestor('H3B104X_13220', 'H4B505J_12900') #Root of GS2
    elif GH == 'BRS':
        outnode = t.get_common_ancestor('LDX55_06330', 'H4B204J_13340') #Root of BRS
    elif GH == 'S1':
        outnode = t.search_nodes(name = 'K2W83_RS06175')[0]
    elif GH == 'S2':
        outnode = t.get_common_ancestor('H4B412M_13220', 'A1805_12800')
    elif GH == 'S3':
        outnode = t.search_nodes(name = 'A1404_13450')[0]
    t.set_outgroup(outnode) #Set the root
    
# =============================================================================
# Set style of tree
# =============================================================================
    ts = TreeStyle() #Create default tree style
    ts.show_branch_length = False # Hide support values
    ts.scale =  2000 #General tree scale
    ts.branch_vertical_margin = 5 # Space between branches
    ts.show_branch_support = False
    ts.show_leaf_name = False # Hide unformatted leaf names
    ts.show_scale = True # Hide tree scale
    
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
                support_face = TextFace(int(n.support), fgcolor = color, fsize = 10) #Format support value text
                n.add_face(support_face, column=0, position='branch-top') #Add the text to the tree
    
# =============================================================================
# Save CDS information for the main plot
# =============================================================================
    #Get strain names
    leaves = t.get_leaves() #Get tree leaves sorted by phylogeny
    lnames = [replace_strain_name(leaf.name) for leaf in leaves] #Get the strain names to be plotted 
    lnames.reverse() #Reverse list with leaf names
    
    gene_dict = {} #Create empty dictionary
    mhpD_dict = {}
    nox_dict = {}
    
    for leaf in lnames: # Loop through leaf names (locus tags)
        for file in sorted(tabs): # Loop through CDS tab files
            strain = leaf.split('_')[0] # Get strain name from locus tag
            if strain in file.replace('-', ''): # If strain name (without -) is in the tab file:
                print(strain, leaf) 
                check = False #Checks to keep the right gene order in the plots 
                oppA_check = False
                border = 'grey' #Color of the border of the plotted CDS
                with open(f'{indir}/{file}', 'r') as gfile: #Open CDS tab file
                    genes = csv.reader(gfile, delimiter='\t') #Read CDS tab file
                    next(genes)
                    
                    for gene in genes: #Loop through genes in CDS tab file
                        gene[1] = gene[1].replace('join(', '')
                        
                        if any(list(gs in gene[0] for gs in GH70_doms + GS2_BRS + GS3 + S1 + S2a + S2b + S3)) or ('mhpD' in gene[0] or 'hpcG' in gene[0] or 'tesE' in gene[0] or 'nox' in gene[0]) or ('oppA' in gene[0] and strain in mhpD_dict.keys() and int(gene[1]) > mhpD_dict[strain][0] and not oppA_check) or (strain != 'MP2'  and strain in nox_dict.keys() and nox_dict[strain][1] + 300 > int(gene[1]) > nox_dict[strain][1]) or (strain == 'DSM' and gene[0] in ['K2W83_RS06170', 'K2W83_RS06165', 'K2W83_RS06190', 'K2W83_RS06195']) or (strain == 'IBH001' and gene[0] in ['LDX55_06315', 'LDX55_06320', 'LDX55_06340', 'LDX55_06345']) or (strain == 'MP2' and gene[0] in ['APS55_RS03835', 'APS55_RS03840', 'APS55_RS03855', 'APS55_RS03860']): #If it is one of the genes of interest
                                
                            if any(list(bs in gene[0] for bs in S2b)): #If it is S2b
                                print(gene[0], 'is S2b') #Print gene type
                                gene[0] = 'S2b' #Change gene name
                                
                            elif any(list(bs in gene[0] for bs in S2a)): #If the gene is S2a
                                print(gene[0], 'is S2a') 
                                gene[0] = 'S2a'
                                
                            elif any(list(bs in gene[0] for bs in S1)): #If the gene is S1
                                print(gene[0], 'is S1')
                                gene[0] = 'S1'
                                
                            elif any(list(gs in gene[0] for gs in GS3)): #If the gene is GS3
                                print(gene[0], 'is GS3')
                                gene[0] = 'GS3'
                                
                            elif any(list(gs in gene[0] for gs in GS1)): #If the gene is GS1
                                print(gene[0], 'is GS1')
                                gene[0] = 'GS1'
                                
                            elif any(list(gs in gene[0] for gs in GS2_BRS)): #If the gene is in GS2_BRS
                                print(gene[0], 'is GS2_BRS') 
                                gene[0] = 'GS2_BRS'
                                    
                            elif any(list(gs in gene[0] for gs in GS2)): #If the gene is in GS2
                                print(gene[0], 'is GS2')
                                gene[0] = 'GS2'
                                        
                            elif any(list(bs in gene[0] for bs in BRS)): #If the gene is BRS
                                print(gene[0], 'is BRS')
                                if 'A1401_12770' in gene[0]: #If the strain has a small extra BRS segment
                                    gene[1] = int(gene[1]) + 11 #Place the start a bit before the segment
                                    
                                gene[0] = 'BRS'
                                
                            elif any(list(bs in gene[0] for bs in S3)): #If the gene is S3
                                print(gene[0], 'is S3')
                                gene[0] = 'S3'
                            
                            else: #If it is other gene of interest
                                print(gene[0], 'is other')
                                if gene[0] == 'hpcG' or gene[0] == 'mhpD' or gene[0] == 'tesE' or gene[0] == 'K2W83_RS06190' or gene[0] == 'LDX55_06340' or gene[0] == 'APS55_RS03840':
                                    gene[0] = 'mhpD'
                                    mhpD_dict[strain] = (int(gene[1]), int(gene[2]))
                                elif gene[0] == 'oppA' or gene[0] == 'K2W83_RS06195' or gene[0] == 'LDX55_06345':
                                    oppA_check = True
                                    gene[0] = 'oppA'
                                elif strain == 'MP2' and check == False and gene[0] == 'APS55_RS03835': #If the strain is MP2 (reversed)
                                    start = int(gene[1]) - padding - gapscale #Set the start position based on oppA
                                    end = start + segment_length #Get the end position
                                    check = True
                                    gene[0] = 'oppA'
                                elif gene[0] == 'nox' or gene[0] == 'K2W83_RS06165' or gene[0] == 'LDX55_06315' or gene[0] == 'APS55_RS03860':
                                    nox_dict[strain] = (int(gene[1]), int(gene[2]))
                                    gene[0] = 'nox'
                                    if strain != 'MP2': #If the strain is not MP2 (not reversed)
                                        start = int(gene[1]) - padding - gapscale #Set the start to the position of nox
                                        end = start + segment_length
                                else:
                                    gene[0] = 'CDS7'


                            format_str = f'Arial|14|black|{gene[0]}' #Text on the CDS to be plotted
                            if gene[0] == 'GS2_BRS': #Gradient fill for multi-GH70 domain proteins
                                col = f'{"rgradient:" + gene_colors["BRS"] + "|" + gene_colors["GS2"]}'
                            else:
                                if gene[0] in gene_colors:
                                    col = gene_colors[gene[0]] #Solid fill for the other CDS
                                else:
                                    col = gene_colors['mhpD']
                                
                            info_list = [int(gene[1])-start, int(gene[2])-start, '[]', 0, 20, border, col, format_str] #Create a list with all the CDS formatting
                            
                            if leaf not in gene_dict.keys(): #If the locus tag is not in the keys of the plotting dictionary
                                gene_dict[leaf] = [info_list] #Add the list with plotting information to the dictionary
                            else: #Otherwise
                                gene_dict[leaf].append(info_list) #Append the information to the dictionary
            
    #Code to adjust the positions of CDS in the plots for each gene
    if GH == 'BRS':
        length_dict = {key:segment_length+1000 for key in gene_dict}
    elif GH == 'GS1':
        length_dict = {key:segment_length+6000 for key in gene_dict}
    elif GH == 'GS2':
        length_dict = {key:segment_length-2000 for key in gene_dict}
    elif GH == 'S3':
        length_dict = {key:segment_length-3000 for key in gene_dict}
    elif GH == 'S2':
        length_dict = {key:segment_length+1000 for key in gene_dict}
    elif GH == 'S1':
        length_dict = {key:segment_length-4000 for key in gene_dict}
    else:
        length_dict = {key:segment_length for key in gene_dict}
    
    #Reverse the plots in all strains except for MP2
    [fix_strand(gene_dict[key]) for key in gene_dict if 'MP2' not in key]
    
    #Create a sequence of desired length to plot the CDS on top    
    seq_dict = {l: f'{"-"*(gapscale)+"A"*(int(length_dict[l])-gapscale+padding-100)}' for l in lnames}
    
# =============================================================================
# Code to add a more intense color to the domains in the plot
# =============================================================================
    gene_domains = {} #Create an empty dictionary
    
    if GH == 'BRS':
        check_dict = gene_dict

    for key in gene_dict.keys(): #Loop through the locus tags that were retrieved
        new_domains = [] #Create empty list to store extra domains
        for element in gene_dict[key]: #Loop through domains to be plotted next to each leaf (locus tag)
            divide = False #Don't divide the gene
            if GH == 'BRS' and f'|{GH}' in element[7] and element[1] - element[0] > 2000 and not key.endswith('_1'): #If it is the BRS plot and the gene is a single-domain BRS
                divide = True #Divide the gene
                print(key, GH, 'Domain BRS')
            elif GH == 'BRS' and '|GS2_BRS' in element[7] and key.endswith('_1'): #If the plot is for BRS and the domain is a BRS in a double-domain GH70
                divide = True
                print(key, GH, 'Domain BRS from GS2_BRS')
            elif GH == 'GS2' and f'|{GH}' in element[7]: #If the plot is for GS2 and the domain is a GS2
                divide = True
                # print(key, GH, 'Domain GS2')
            elif GH == 'GS1' and f'|{GH}' in element[7]: #If the plot is for GS1 and the domain is a GS1
                divide = True
                # print(key, GH, 'Domain GS1')
            elif GH == 'S1' and f'|{GH}' in element[7]:
                divide = True
                # print(key, GH, f'Domain {GH}')
            elif GH == 'S2' and f'|{GH}' in element[7]:
                if (key in S2a and 'S2a' in element[7]) or (key in S2b and 'S2b' in element[7]):
                    divide = True
                # print(key, GH, f'Domain {GH}')
            elif GH == 'S3' and f'|{GH}' in element[7]:
                divide = True
                # print(key, GH, f'Domain {GH}')
                
                
            if divide: #If the gene should be divided
                new_domain = element.copy() #Clone CDS to be plotted
                new_domain[0] = element[0] + domain_dict[key][0] - 1 #Modify start of segment to start of the domain
                new_domain[1] = element[0] + domain_dict[key][1] - 1 #Modify end of segment to end of the domain
                new_domain[7] = new_domain[7].replace('GS2_BRS', GH) #Modify text of segment if it is a GS2_BRS gene
                new_domain[5] = 'black' #Modify edge color
                for GH_typ in GH70_types+GH32_types: #Loop through GH types and change the domain colors to brighter colors
                    if f'|{GH_typ}' in element[7] and 'GS2_BRS' not in element[7]:
                        element[6] = alt_colors[GH_typ]
                        if 'S2b' in element[7]:
                            element[6] = alt_colors['S2b']
                    else:
                        element[6] = element[6].replace(gene_colors[GH_typ], alt_colors[GH_typ])

                if 'GS2_BRS' not in element[7]: #Remove the name of the gene from the segments of the gene beyond the domain if the gene is not a GS2_BRS
                    element[7] = element[7].replace(gene_types[key], '')
                    element[7] = element[7].replace('|a', '|')
                    element[7] = element[7].replace('k|b', 'k|')
                next_domain = element.copy() #Clone CDS to be plotted
                if GH == 'GS2' and 'GS2_BRS' in element[7]: #If the gene is GS2_BRS, remove the gene name from the shorter segment to be plotted
                    element[7] = element[7].replace('GS2_BRS', '')
                elif GH == 'BRS' and 'GS2_BRS' in next_domain[7]:
                    next_domain[7] = next_domain[7].replace('GS2_BRS', '')
                next_domain[0] = new_domain[1] + 1 #Create a new segment to plot that covers from the end of the domain to the end of the CDS
                element[1] = new_domain[0] - 1 #Modify the original segment so that it ends at the start of the domain
                new_domains.append(element) #Add new CDS segments to the list
                new_domains.append(next_domain)
                new_domains.append(new_domain)
            else: #If the gene shouldn't be divided, change the colors to less bright colors
                for GH_typ2 in GH70_types+GH32_types:
                    if f'|{GH_typ2}' in element[7] and 'GS2_BRS' not in element[7]:
                        element[6] = alt_colors[GH_typ2]
                    elif GH_typ2 in element[7]:
                        element[6] = element[6].replace(gene_colors[GH_typ2], alt_colors[GH_typ2])
                    if 'S2b' in element[7]:
                        element[6] = alt_colors['S2b']
                new_domains.append(element)
        gene_dict[key] = new_domains # I need to add some kind of check where I check that the key gene type is the same as the gene type to add domains, and then I need to check that it is the right domain when there are two genes that are classified as the same

# =============================================================================
# Code to center domains
# =============================================================================
    sum_dict = {} #Create an empty dictionary
    for key, value in gene_dict.items(): #Loop through information in the dictionary
        for entry in gene_dict[key]: #Loop through CDS to be plotted for a specific leaf of the phylogeny
            if f'|{GH}' in entry[7] and (any(color in entry[6] for color in gene_colors.values()) or (any(color in entry[6] for color in alt_colors.values())) and GH != 'BRS'): #If the type of a gene matches the gene type that should be highlighted on the tree
                sum_dict[key] = entry[0] #Store the start position of the gene

    max_value = max(sum_dict.values()) #Get the highest starting position
    sub_dict = {k: max_value - v for k, v in sum_dict.items()} #Calculate the difference to the start position for all remaining genes that belong to the category of interest

    
    for key in gene_dict.keys(): #Loop through leaves in the phylogeny
        for element in gene_dict[key]: #Loop through CDS/CDS segments to be plotted
            element[0] += sub_dict[key] #Add the difference to the maximum value to all the genes to be plotted so that the domains are lined up
            element[1] += sub_dict[key]
    

# =============================================================================
# Create the plot
# =============================================================================
    for leaf in leaves: #Loop through leaves in the tree
        lname = replace_strain_name(leaf.name) #Get the right name for each leaf
        strain = lname.split('_')[0] #Get the name of the strain
        color = leaf_color.get(strain, None) #Color leaves according to strain phylogroup
        name_face = TextFace(lname, fgcolor = color, fsize = 18) #Create text with corrected locus tag
        leaf.add_face(name_face, column=0, position='branch-right') #Add formatted leaf names

        seqFace = SeqMotifFace(seq_dict[lname], motifs = gene_dict[lname], 
                                seq_format = 'line', gap_format = 'blank', 
                                scale_factor = 0.03) #Create variable storing the CDS to plot
        (t & f'{leaf.name}').add_face(seqFace, 2, 'aligned') #Add CDS representation to leaf
        
# =============================================================================
# Print tree
# =============================================================================
    t.ladderize(1)
    t.render(f'{outdir}/{outplot}', tree_style = ts) #Render tree in PNG format
    t.render(f'{outdir}/{outplot.replace(".png", ".tiff")}', tree_style = ts) #Render tree in TIFF format
    t.render(f'{outdir}/{outplot.replace(".png", ".svg")}', tree_style = ts) #Render tree in SVG format
    
