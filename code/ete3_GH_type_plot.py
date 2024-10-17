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

color_dict = {'1': 'blue', '2': 'orange', '3': 'green', '4': 'red', '5': 'purple',
              '6': 'magenta', '7': 'brown'}

hpr_colors = ['#ADB8CE', '#8C95A7', '#687287', '#505B6F', '#58506F', '#64506F',
              '#6F5150', '#6F6250', '#696F50']

gene_colors = {'GS1': '#FF7575', 'GS2': '#FF75B6', 'BRS': '#E875FF',
               'S1': '#FFE0C3', 'S2a': '#FFF0CB', 'S2b': '#FFF0CB', 
               'S3': '#FFFCD0', 'mhpD': '#CCCCCC', 'hpcG': '#CCCCCC', 
               'oppA': '#CCCCCC', 'GS3': '#FFC2D0', 'GS4': '#FF8A75'}

alt_colors = {'GS1': '#FFC0C0', 'GS2': '#FFC0DE', 'BRS': '#F4C0FF'}

# =============================================================================
# Set target region and load input files
# =============================================================================
# TO DO: Annotate script and integrate into Snakemake

segment_length = 25000
gapscale = 1000
padding = 500
config_file = os.path.expanduser('~') + '/GH_project/config.yml'
indir =  os.path.expanduser('~') + '/GH_project/plots/tabfiles'
outdir =  os.path.expanduser('~') + '/GH_project/plots/trees'
GH_types = ['BRS', 'GS1', 'GS2']
domain_path =  os.path.expanduser('~') + '/interproscan'
doub_GS1_strains = ['H3B2-06M', 'H4B1-11J', 'H4B5-05J', 'H4B4-06M']

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
        
    for el in reversed(my_info_list):
        my_index = my_info_list.index(el)
        if my_index == len(my_info_list)-1:
            start = padding + gapscale
        else:
            start = pass_next - el[1] + end
        end = start + (el[1] - el[0])
        pass_next = el[0]
        my_info_list[my_index][0] = start
        my_info_list[my_index][1] = end
    return my_info_list

# =============================================================================
# Read types from Snakemake
# =============================================================================
# Read locus tags of different gene subsets
with open(config_file) as conf:
    py_config = yaml.safe_load(conf)
    GS1_repr = py_config['GS1_repr']
    GS2_repr = py_config['GS2_repr']
    BRS_repr = py_config['BRS_repr']
    GH70_repr = list(np.unique(GS1_repr + GS2_repr + BRS_repr))
    
    S1_repr = remove_minus(py_config['S1_repr'])
    S2a_repr = remove_minus(py_config['S2a_repr'])
    S2b_repr = remove_minus(py_config['S2b_repr'])
    S3_repr = remove_minus(py_config['S3_repr'])
    
    
    GS1 = py_config['GS1']
    GS2 = py_config['GS2']
    GS3 = py_config['GS3']
    GS3_repr = [tag.split('_')[0] for tag in GS3]
    S1 = py_config['S1']
    S2a = py_config['S2a']
    S2b = py_config['S2b']
    S3 = py_config['S3']
    BRS = ['A1401_12770'] + py_config['BRS']
    GH70_doms = GS1 + GS2 + BRS
    strains = py_config['representatives']
    s_nominus = [strain.replace('-', '') for strain in strains]
    # Retrieve BRS domains of the GS2_BRS proteins
    GS2_BRS =  [GH_gene.replace('_2', '') for GH_gene in GS2 if GH_gene.replace('_2', '_1') in BRS] #+ [replace_strain_name(GH_gene.replace('_1', '')) for GH_gene in BRS if GH_gene.replace('_1', '') in GS2]
    # Retrieve the most different GS1 group in strains with two GS1
    curated_strains = [strain.replace('-', '').replace('DSMZ12361', 'DSM') for strain in strains]
    GS1_strains = [replace_strain_name(gs).split('_')[0] for gs in GS1 if replace_strain_name(gs).split('_')[0] in curated_strains]
    double_GS1 = [strain for strain, count in Counter(GS1_strains).items() if count > 1]
    GS1_dom = [glu1 for glu1 in GS1 if any(dgs in glu1 for dgs in double_GS1)]
    max_GS1 = {strain: gsd for strain in double_GS1 for gsd in GS1_dom if strain in gsd}

domain_dict = {}
for GH in GH_types:
    print(f'\n\nGH70 {GH} genes...', end = '\n')
    treefile =  os.path.expanduser('~') + f'/GH_project/data/fasta/GH70/trees/{GH}_repset.mafft.faa.treefile'
    # count_file = snakemake.input.counts
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outplot = f'{GH}_phylogeny.png'

# =============================================================================
# Save gene types to dictionary
# =============================================================================
    # Create dictionary assigning gene types
    gene_types = {}
    for gene_dom in GH70_doms:
        gene = replace_strain_name(gene_dom)
        if gene_dom in GS1:
            gene_types[gene] = 'GS1'
        elif gene_dom in GS2 or (gene_dom[:-2] + '_2' in GS2):
            if gene_dom[:-2] + '_1' in BRS:
                gene_types[gene] = 'GS2_BRS'
                # print(gene[:-2])
            else: gene_types[gene] = 'GS2'
        elif gene_dom in BRS:
            gene_types[gene] = 'BRS'
        elif gene_dom in GS3:
            gene_types[gene] = 'GS3'
        else: print('Something is VERY wrong here...')
            
    tabs = [file for file in os.listdir(indir) if any(list(strain in file for strain in strains))] #snakemake.input.tabs #Maybe I should regenerate the tabs adding locus tag and annotation separately...

# =============================================================================
# Read Interproscan to get the start and end positions of domains within each gene   
# =============================================================================
    domains = [file for file in os.listdir(domain_path) if any(list(strain.replace('Fhon2', 'fhon2') in file for strain in strains))]
    for domain in domains:
        domain_file = f'{domain_path}/{domain}'    
        with open(domain_file, 'r') as dfile:
            freader = csv.reader(dfile, delimiter = '\t')
            for line in freader:
                gene_locus = replace_strain_name(line[0]).replace('-', '').replace('fhon2', 'Fhon2')
                if ((line[0] in GH70_doms) or (gene_locus.upper() in GH70_doms) or (f'{gene_locus.upper()}_1' in GH70_doms) or (f'{line[0]}_1' in GH70_doms) or ('MP2' in domain_file and int(gene_locus.split('_')[1]) < 14000)) and 'Glycosyl hydrolase family 70' in line[5]:
                    if 'MP2' in domain_file:
                        gene_locus = gene_locus.replace('_13350', '_03850').replace('_13360', '_03845')
                    pos = (int(line[6])*3, int(line[7])*3)
                    if gene_locus in domain_dict.keys() and domain_dict[gene_locus] != pos:
                        if 'DSM' not in domain_file:
                            domain_dict[f'{gene_locus}_1'] = domain_dict[gene_locus]
                            del domain_dict[gene_locus]
                            gene_locus = f'{gene_locus}_2'
                        else:
                            domain_dict[f'{gene_locus}_2'] = domain_dict[gene_locus]
                            del domain_dict[gene_locus]
                            gene_locus = f'{gene_locus}_1'
                    domain_dict[gene_locus] = pos
    
# =============================================================================
# Create rooted tree from treefile
# =============================================================================
    t = Tree(treefile, format = 0) #Create tree
    
    #Root tree
    if GH == 'GS1':
        outnode = 'A1001_12310'
    elif GH == 'GS2':
        outnode = t.get_common_ancestor('H3B104X_13220', 'H4B504J_13480') # GS2
    else:
        outnode = t.get_common_ancestor('LDX55_06330', 'H4B204J_13340') # BRS
    t.set_outgroup(outnode)
    
# =============================================================================
# Set style of tree
# =============================================================================
    ts = TreeStyle() #Create default tree style
    ts.show_branch_length = False # Hide support values
    ts.scale =  500 #General tree scale
    ts.branch_vertical_margin = 5 # Space between branches
    ts.show_branch_support = False
    ts.show_leaf_name = False # Hide unformatted leaf names
    ts.show_scale = False # Hide tree scale
    
    #Add legend
    ts.legend_position = 3 #Bottom left
    ts.legend.add_face(TextFace(' '*58, fsize = 10), column = 0)
    
    #Formatting of nodes (and branch lines)
    ns = NodeStyle()
    ns['size'] = 0 #Nodes are not visible
    ns['vt_line_width'] = 2
    ns['hz_line_width'] = 2
    ns['hz_line_type'] = 0
    for n in t.traverse():
        n.set_style(ns)
        if n not in t.get_leaves() and n.support > 1:
            if n.support >= 95:
                color = 'black' #'#0E9E00'
            else: color = 'dimgrey'
            if n.support >= 80:
                support_face = TextFace(int(n.support), fgcolor = color, fsize = 10)
                n.add_face(support_face, column=0, position='branch-top')
    
# =============================================================================
# Save CDS information for the main plot
# =============================================================================
    #Get strain names
    leaves = t.get_leaves() #Sorted by phylogeny
    lnames = [replace_strain_name(leaf.name) for leaf in leaves]
    lnames.reverse()
    
    gene_dict = {}
    
    #Add leaf names and plot region
    for leaf in lnames: # Loop through leaf names (locus tags)
        for file in sorted(tabs): # Loop through CDS tab files
            strain = leaf.split('_')[0] # Get strain name from locus tag
            if strain in file.replace('-', ''): # If strain name (without -) is in the tab file:
                print(strain, leaf)
                check = False
                GS1_check = 0
                S3_check = False
                if strain in S3_repr:
                    S3_check = True
                S2a_check = False
                print(strain, S3_check)
                border = 'grey'
                with open(f'{indir}/{file}', 'r') as gfile:
                    genes = csv.reader(gfile, delimiter='\t')
                    
                    for gene in genes:
                        
                        if any(list(gs in gene[0] for gs in GH70_doms + GS2_BRS + GS3 + S1 + S2a + S2b + S3)) or (S3_check == True and S2a_check == True and ('oppA' in gene[0] or 'mhpD' in gene[0] or 'hpcG' in gene[0])): #GS1+GS2+BRS)):
                                
                            if any(list(bs in gene[0] for bs in S2b)):
                                print(gene[0], 'is S2b')
                                start = int(gene[1]) - padding - gapscale
                                end = start + segment_length
                                GS1_check += 1
                                gene[0] = 'S2b'
                                
                            elif any(list(bs in gene[0] for bs in S2a)):
                                print(gene[0], 'is S2a')
                                if GS1_check == 0:
                                    start = int(gene[1]) - padding - gapscale
                                    end = start + segment_length  
                                GS1_check += 1
                                S2a_check = True
                                gene[0] = 'S2a'
                                
                            elif any(list(bs in gene[0] for bs in S1)):
                                print(gene[0], 'is S1')
                                start = int(gene[1]) - padding - gapscale
                                end = start + segment_length
                                GS1_check += 1
                                gene[0] = 'S1'
                                
                            elif any(list(gs in gene[0] for gs in GS3)):
                                print(f'Processing gene {gene[0]}...')
                                print(gene[0], 'is GS3')
                                start = int(gene[1]) - padding - gapscale
                                end = start + segment_length
                                GS1_check += 1
                                gene[0] = 'GS3'
                                
                            elif any(list(gs in gene[0] for gs in GS1)):
                                print(gene[0], 'is GS1')
                                if strain != 'MP2':
                                    if GS1_check == 0:
                                        start = int(gene[1]) - padding - gapscale
                                        end = start + segment_length
                                    check = True
                                    GS1_check += 1
                                gene[0] = 'GS1'
                                
                            elif any(list(gs in gene[0] for gs in GS2_BRS)):
                                print(gene[0], 'is GS2_BRS')
                                gene[0] = 'GS2_BRS'
                                    
                            elif any(list(gs in gene[0] for gs in GS2)):
                                print(gene[0], 'is GS2')
                                gene[0] = 'GS2'
                                if strain == 'MP2' and check == False:
                                    start = int(gene[1]) - padding - gapscale
                                    end = start + segment_length
                                    check = True
                                        
                            elif any(list(bs in gene[0] for bs in BRS)):
                                print(gene[0], 'is BRS')
                                if 'A1401_12770' in gene[0]:
                                    gene[1] = int(gene[1]) + 11
                                    
                                gene[0] = 'BRS'
                                
                            elif any(list(bs in gene[0] for bs in S3)):
                                print(gene[0], 'is S3')
                                gene[0] = 'S3'
                            
                            else:
                                print(gene[0], 'is other')


                            format_str = f'Arial|14|black|{gene[0]}' #Text on the gene
                            if gene[0] == 'GS2_BRS': #I use a gradient for multi-GH70 domain proteins.
                                col = f'{"rgradient:" + gene_colors["BRS"] + "|" + gene_colors["GS2"]}'
                            else:
                                col = gene_colors[gene[0]]
                                
                            info_list = [int(gene[1])-start, int(gene[2])-start, '[]', 0, 20, border, col, format_str]
                            
                            if leaf not in gene_dict.keys():
                                gene_dict[leaf] = [info_list]
                            else:
                                gene_dict[leaf].append(info_list)
            
    if GH == 'BRS':
        length_dict = {key:segment_length for key in gene_dict}
    elif GH == 'GS1':
        length_dict = {key:segment_length+2000 for key in gene_dict}
    elif GH == 'GS2':
        length_dict = {key:segment_length-5500 for key in gene_dict}
    
    [fix_strand(gene_dict[key]) for key in gene_dict if 'MP2' not in key]
    
    seq_dict = {l: f'{"-"*(gapscale)+"A"*(int(length_dict[l])-gapscale+padding-100)}' for l in lnames} #Create a sequence of desired length
    
# =============================================================================
# Code to add a more intense color to the domains in the plot
# =============================================================================
    gene_domains = {}

    for key in gene_dict.keys():
        print(key)
        new_domains = []
        for element in gene_dict[key]:
            divide = False
            if GH == 'BRS' and key in GS2_BRS and GH in element[7] and element == gene_dict[key][-1]:
                divide = True
                print(key, GH, 'True1')
            elif GH == 'BRS' and f'|{GH}' in element[7] and element[1] - element[0] > 2000 and ((key == 'H4B204J_13340' or key == 'IBH001_06330') or ('GS2_BRS' not in gene_dict[key][-1][7])):
                divide = True
                print(key, GH, 'True2')
            elif GH == 'BRS' and '|GS2_BRS' in element[7] and key.endswith('_1'):
                divide = True
                print(key, GH, 'True4')
            elif GH == 'GS2' and f'|{GH}' in element[7]:
                divide = True
            elif GH == 'GS1' and f'|{GH}' in element[7]:
                print(key, element)
                check2 = (key in max_GS1.values()) + (key.split('_')[0] in S2b_repr) + (key.split('_')[0] in S2a_repr) + (key.split('_')[0] in S1_repr) + (key.split('_')[0] in GS3_repr)
                print(key, check2)
                if element == gene_dict[key][check2] or 'MP2' in key:
                    divide = True
                print(key, GH, 'True3')
                
            if divide:
                new_domain = element.copy()
                new_domain[0] = element[0] + domain_dict[key][0] - 1
                new_domain[1] = element[0] + domain_dict[key][1] - 1
                new_domain[7] = new_domain[7].replace('GS2_BRS', GH)
                new_domain[5] = 'black'
                for GH_typ in GH_types:
                    if GH_typ in element[7] and 'GS2_BRS' not in element[7]:
                        element[6] = alt_colors[GH_typ]
                    else:
                        element[6] = element[6].replace(gene_colors[GH_typ], alt_colors[GH_typ])
                if 'GS2_BRS' not in element[7]:
                    element[7] = element[7].replace(gene_types[key], '')
                next_domain = element.copy()
                if GH == 'GS2' and 'GS2_BRS' in element[7]:
                    element[7] = element[7].replace('GS2_BRS', '')
                elif GH == 'BRS' and 'GS2_BRS' in next_domain[7]:
                    next_domain[7] = next_domain[7].replace('GS2_BRS', '')
                next_domain[0] = new_domain[1] + 1
                element[1] = new_domain[0] - 1
                new_domains.append(element)
                new_domains.append(next_domain)
                new_domains.append(new_domain)
            else:
                for GH_typ2 in GH_types:
                    if GH_typ2 in element[7] and 'GS2_BRS' not in element[7]:
                        element[6] = alt_colors[GH_typ2]
                    elif GH_typ2 in element[7]:
                        element[6] = element[6].replace(gene_colors[GH_typ2], alt_colors[GH_typ2])
                new_domains.append(element)
        gene_dict[key] = new_domains # I need to add some kind of check where I check that the key gene type is the same as the gene type to add domains, and then I need to check that it is the right domain when there are two genes that are classified as the same

# =============================================================================
# Code to center domains
# =============================================================================
    sum_dict = {}
    for key,value in gene_dict.items():
        for entry in gene_dict[key]:
            if GH in entry[7] and any(color in entry[6] for color in gene_colors.values()):
                sum_dict[key] = entry[0]
    max_value = max(sum_dict.values())
    sub_dict = {k: max_value - v for k, v in sum_dict.items()}
    
    for key in gene_dict.keys():
        for element in gene_dict[key]:
            element[0] += sub_dict[key]
            element[1] += sub_dict[key]
    

# =============================================================================
# Create the plot
# =============================================================================
    for leaf in leaves:
        lname = replace_strain_name(leaf.name)
        strain = lname.split('_')[0]
        color = leaf_color.get(strain, None) #Color leaves according to strain phylogroup
        name_face = TextFace(lname, fgcolor = color, fsize = 18)
        leaf.add_face(name_face, column=0, position='branch-right') #Add formatted leaf names

        if lname in seq_dict.keys():
            gene_list = gene_dict[lname]
            seqFace = SeqMotifFace(seq_dict[lname], motifs = gene_list, seq_format = 'line', gap_format = 'blank', scale_factor = 0.03) #Original 0.04, Add genomic segment info to node
            (t & f'{leaf.name}').add_face(seqFace, 2, 'aligned')
        
# =============================================================================
# Print tree
# =============================================================================
    t.render(f'{outdir}/{outplot}', tree_style = ts)
    t.render(f'{outdir}/{outplot.replace(".png", ".tiff")}', tree_style = ts)
    t.render(f'{outdir}/{outplot.replace(".png", ".svg")}', tree_style = ts) 
