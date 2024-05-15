#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates a plot of the phylogeny of A. kunkeei strains, with
a representation of the GS1, GS2 and BRS genes that are present. The
phylogenies used include representative strains only, and only the GS2
and BRS phylogenies were included.

@author: Marina Mota Merlo
"""

from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace
import colorsys
import random
#Use PhyloTree
from Bio import SeqIO
import re, csv, os, logging, traceback
from matplotlib.colors import to_hex
from math import sqrt
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

leaf_color = {'A0901': '#2AA380', 'A1001': '#A3C615', 'A1002': '#2AA380',
              'A1003': '#2C85D8', 'A1201': '#2AA380', 'A1202': '#21C1D1',
              'A1401': '#21C1D1', 'A1404': '#7ACA2E', 'A1802': '#2AA380',
              'A1803': '#2AA380', 'A1805': '#21C1D1', 'A2001': '#2AA380',
              'A2002': '#A3C615', 'A2003': '#A3C615', 'A2101': '#21C1D1',
              'A2102': '#2AA380', 'A2103': '#2C85D8', 'G0101': '#2C85D8', 
              'G0102': '#2C85D8', 'G0103': '#2C85D8', 'G0401': '#2C85D8',
              'G0402': '#2C85D8', 'G0403': '#21C1D1', 'G0404': '#2C85D8',
              'G0405': '#2C85D8', 'G0406': '#21C1D1', 'G0407': '#2C85D8',
              'G0408': '#2C85D8', 'G0410': '#2C85D8', 'G0412': '#2C85D8',
              'G0414': '#2C85D8', 'G0415': '#2C85D8', 'G0417': '#2C85D8',
              'G0420': '#21C1D1', 'G0601': '#2C85D8', 'G0602': '#2C85D8',
              'G0702': '#2C85D8', 'G0801': '#2C85D8', 'G0802': '#2C85D8',
              'G0803': '#2C85D8', 'G0804': '#2C85D8', 'Fhon2': '#2C85D8',
              'H1B104J': '#2C85D8', 'H1B105A': '#2C85D8', 'H1B302M': '#2C85D8',
              'H2B105J': '#2AA380', 'H3B101A': '#2C85D8', 'H3B101J': '#2C85D8',
              'H3B101X': '#2AA380', 'H3B102A': '#2C85D8', 'H3B102X': '#2AA380',
              'H3B103J': '#2C85D8', 'H3B103X': '#2AA380','H3B103M': '#2C85D8', 
              'H3B104J': '#2C85D8', 'H3B104X': '#2C85D8', 'H3B107A': '#2C85D8', 
              'H3B109M': '#2C85D8', 'H3B110M': '#2C85D8', 'H3B111A': '#2C85D8', 
              'H3B111M': '#2C85D8', 'H3B202M': '#2AA380', 'H3B202X': '#2C85D8', 
              'H3B203J': '#2C85D8', 'H3B203M': '#2AA380', 'H3B204J': '#2C85D8', 
              'H3B204M': '#2AA380', 'H3B205J': '#2C85D8', 'H3B206M': '#2AA380', 
              'H3B207X': '#2C85D8', 'H3B208X': '#2C85D8', 'H3B209X': '#21C1D1', 
              'H4B101A': '#2AA380', 'H4B102A': '#2AA380', 'H4B103J': '#2AA380', 
              'H4B104A': '#2AA380', 'H4B111J': '#2AA380', 'H4B114J': '#2AA380', 
              'H4B116J': '#2AA380', 'H4B202J': '#2C85D8', 'H4B203M': '#2AA380', 
              'H4B204J': '#2C85D8', 'H4B205J': '#2C85D8', 'H4B206J': '#2AA380', 
              'H4B210M': '#2AA380', 'H4B211M': '#2C85D8', 'H4B303J': '#2AA380', 
              'H4B402J': '#2C85D8', 'H4B403J': '#2AA380', 'H4B404J': '#2AA380', 
              'H4B405J': '#2AA380', 'H4B406M': '#2AA380', 'H4B410M': '#2AA380', 
              'H4B411M': '#2AA380', 'H4B412M': '#2C85D8', 'H4B501J': '#2C85D8', 
              'H4B502X': '#2C85D8', 'H4B503X': '#2C85D8', 'H4B504J': '#21C1D1', #'#2AA380', 
              'H4B505J': '#21C1D1', 'H4B507J': '#2C85D8', 'H4B507X': '#2C85D8', 
              'H4B508X': '#2C85D8', 'MP2': '#21C1D1', 'IBH001': 'black', 
              'DSM': 'black'}

color_dict = {'1': 'blue', '2': 'orange', '3': 'green', '4': 'red', '5': 'purple',
              '6': 'magenta', '7': 'brown'}

hpr_colors = ['#ADB8CE', '#8C95A7', '#687287', '#505B6F', '#58506F', '#64506F',
              '#6F5150', '#6F6250', '#696F50']

gene_colors = {'nox': '#7E9C07', 'GS1': '#FF0000', 'tes': '#BB9900', 'opp': '#2FAD26',
               'sbn': '#6FB039', 'eps': '#39B0AB', 'ada': '#1AA682', 'fab': '#20B5D6',
               'ybi': '#92C62A', 'yhf': '#21E190', 'gal': '#B8C14C', 'acc': '#00E7C8',
               'mhp': '#0FD168', 'mnt': '#0FA2D1', 'umu': '#0F7FD1', 'ywq': '#00AB72',
               'tag': '#0081B2', 'yhd': '#76B200', 'mdt': '#33CD8E', 'maa': '#117EDA',
               'bet': '#69DA11', 'rnh': '#35DA11', 'kef': '#11DA60', 'ydg': '#11DA9D',
               'yic': '#53C1A0', 'cfi': '#53C173', 'dpp': '#6BC153', 'drr': '#53C1A0',
               'cap': '#00AB72', 'hpc': '#5399C1', 'nam': '#C1BB53', 'gsi': '#ADC153',
               'gla': '#BEE800', 'yfk': '#63E800', 'COQ': '#00E85F', 'lem': '#00E8A6',
               'ald': '#00E5E8', 'yvg': '#00BBE8', 'acp': '#0082E8', 'thi': '#00E8D0',
               'ins': '#7FE800', 'adh': '#E8DA00', 'met': '#E8C500', 'dsr': '#E85800',
               'pur': '#C8B817', 'yra': '#17ABC8', 'nch': '#B600FF', 'BRS': '#D930BD',
               'gtf': '#000000', 'GS2': '#FF009B', 'glf': '#54C190', 'mur': '#48CDA7',
               'pgl': '#4ACD48', 'asd': '#48CDC3', 'dap': '#7ACD48', 'lys': '#A8CB3C',
               'ppc': '#D2D518', 'ydh': '#18BED5', 'nai': '#18D581', 'pep': '#18D56B',
               'gly': '#18D548', 'sho': '#AAD518', 'rlm': '#C4BF2A', 'yfe': '#2AC4B8',
               'hom': '#2AAAC4', 'dac': '#98C42A', 'arn': '#67C42A', 'asp': '#2EDDA5',
               'suf': '#95DD2E', 'nrd': '#D1CD00', 'ela': '#00D165', 'nap': '#00D178',
               'ade': '#00D198', 'gar': '#00D1B7', 'thl': '#00C7D1', 'mco': '#D1A800',
               'gnd': '#BED100', 'zwf': '#9ED100', 'lep': '#75D100', 'yda': '#42D100',
               'gmu': '#D7A300', 'add': '#D7C400', 'crc': '#00D72A', 'slm': '#00D775',
               'clp': '#00D793', 'ser': '#00D75B', 'obg': '#00D782', 'uvr': '#00D7A6',
               'mac': '#00D7C0', 'mva': '#00C7D7', 'apt': '#00A0D7', 'rec': '#0079D7',
               'oxy': '#CACD33', 'npr': '#AECD33', 'znu': '#97CD33', 'hba': '#80CD33',
               'hrt': '#68CD33', 'btu': '#33CD89', 'ten': '#33CDAC', 'S1 ': '#FF9600',
               'S2a': '#FFC800', 'S2b': '#FFC800', 'S3 ': '#FFEF00', 'skf': '#0082E8',
               'din': '#00D198'}

alt_colors = {'GS1': '#C7727D', 'GS2': '#C772A8', 'BRS': '#C772C2'}

# =============================================================================
# Set target region and load input files
# =============================================================================

# Idea: Re-do the trees and use A1003 to root all of them
# TO DO: Add S1-3 to the GS1 tree (and the S1 in DSM to the GS2/BRS trees)
# TO DO: Center domains

segment_length = 18000
gapscale = 1000
padding = 500
config_file = '../config.yml'
indir = '../plots/tabfiles'
outdir = '../plots/trees'
GH_types = ['GS1', 'GS2', 'BRS']
domain_path = '../../interproscan'
doub_GS1_strains = ['H3B2-06M', 'H4B1-11J', 'H4B5-05J', 'H4B4-06M']
# comp_dir = os.path.commonpath(tabs)

def replace_strain_name(locus_tag):
    locus_tag = locus_tag.replace('LDX55', 'IBH001').replace('APS55', 'MP2').replace('K2W83', 'DSM').replace('RS', '').replace('FHON', 'Fhon').replace('AKU', '')
    return locus_tag

def fix_strand(my_info_list):  
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

domain_dict = {}
for GH in GH_types:
    print(f'\n\nGH70 {GH} genes...', end = '\n')
    treefile = f'../data/fasta/GH70/trees/{GH}_repset.mafft.faa.treefile'
    # count_file = snakemake.input.counts
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outplot = f'{GH}_phylogeny.png'

# =============================================================================
# Read types from Snakemake
# =============================================================================
    
    # GS2 = snakemake.params.GS2
    # BRS = snakemake.params.BRS
    
    with open(config_file) as conf:
        py_config = yaml.safe_load(conf)
        GS1_repr = py_config['GS1_repr']
        GS2_repr = py_config['GS2_repr']
        BRS_repr = py_config['BRS_repr']
        GH70_repr = list(np.unique(GS1_repr + GS2_repr + BRS_repr))
        GS1 = py_config['GS1']
        GS2 = py_config['GS2']
        BRS = ['A1401_12770'] + py_config['BRS']
        GH70_doms = GS1 + GS2 + BRS
        strains = py_config['representatives']
        # Retrieve BRS domains of the GS2_BRS proteins
        GS2_BRS =  [replace_strain_name(GH_gene.replace('_2', '')) for GH_gene in GS2 if GH_gene.replace('_2', '') in BRS] + [replace_strain_name(GH_gene) for GH_gene in BRS if GH_gene.replace('_2', '') in GS2]
        # Retrieve the most different GS1 group in strains with two GS1
        curated_strains = [strain.replace('-', '').replace('DSMZ12361', 'DSM') for strain in strains]
        GS1_strains = [replace_strain_name(gs).split('_')[0] for gs in GS1 if replace_strain_name(gs).split('_')[0] in curated_strains]
        double_GS1 = [strain for strain, count in Counter(GS1_strains).items() if count > 1]
        GS1_dom = [glu1 for glu1 in GS1 if any(dgs in glu1 for dgs in double_GS1)]
        max_GS1 = {strain: gsd for strain in double_GS1 for gsd in GS1_dom if strain in gsd}
    
    # Create dictionary assigning gene types
    gene_types = {}
    for gene_dom in GH70_doms:
        gene = replace_strain_name(gene_dom)
        if gene_dom in GS1:
            gene_types[gene] = 'GS1'
        elif gene_dom in GS2:
            if gene_dom in BRS:
                gene_types[gene] = 'GS2_BRS'
            else: gene_types[gene] = 'GS2'
        else:
            gene_types[gene] = 'BRS'
            
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
                if ((line[0] in GH70_doms) or (gene_locus.upper() in GH70_doms) or ('MP2' in domain_file and int(gene_locus.split('_')[1]) < 14000)) and 'Glycosyl hydrolase family 70' in line[5]:
                    if 'MP2' in domain_file:
                        gene_locus = gene_locus.replace('_13350', '_03850').replace('_13360', '_03845')
                    pos = (int(line[6])*3, int(line[7])*3)
                    if gene_locus in domain_dict.keys() and domain_dict[gene_locus] != pos:
                        gene_locus = f'{gene_locus}_2'
                    domain_dict[gene_locus] = pos
    
# =============================================================================
# Create rooted tree from treefile
# =============================================================================
    t = Tree(treefile, format = 0) #Create tree
    
    #Root tree
    if GH == 'GS1':
        outnode = 'A1003_12540'
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
        # n = 0
        for file in sorted(tabs): # Loop through CDS tab files
            strain = leaf.split('_')[0] # Get strain name from locus tag
            if strain in file.replace('-', ''): # If strain name (without -) is in the tab file:
                print(strain, leaf)
                check = False
                GS1_check = 0
                
                with open(f'{indir}/{file}', 'r') as gfile:
                    genes = csv.reader(gfile, delimiter='\t')
                    
                    for gene in genes:
                        border = 'white'
                        
                        if any(list(gs in gene[0] for gs in GS1+GS2+BRS)):
                            if any(list(gs in gene[0] for gs in GS1)) or (gene[0] == 'gtfC' and check == False) or (strain in doub_GS1_strains and GS1_check < 2):
                                print(gene[0], 'is GS1')
                                if strain != 'MP2':
                                    if GS1_check == 0:
                                        start = int(gene[1]) - padding - gapscale
                                        end = start + segment_length
                                    check = True
                                    GS1_check += 1
                                gene[0] = 'GS1'
                                
                            elif any(list(gs.replace('_2', '') in gene[0] for gs in GS2)) or (gene[0] == 'gtfC' and check == True):
                                # if GH == 'GS2':
                                #     border = 'black'
                                if any(list(bs.replace('_2', '') in gene[0] for bs in BRS)):
                                    print(gene[0], 'is GS2_BRS')
                                    # if GH == 'BRS' and replace_strain_name(gene[0]) == leaf.replace('_2', ''):
                                    #     border = 'black'
                                    gene[0] = 'GS2_BRS'
                                    
                                else:
                                    print(gene[0], 'is GS2')
                                    gene[0] = 'GS2'
                                    if strain == 'MP2' and check == False:
                                        start = int(gene[1]) - padding - gapscale
                                        end = start + segment_length
                                        check = True
                                        
                            elif any(list(bs in gene[0] for bs in BRS)):
                                print(gene[0], 'is BRS')
                                # if GH == 'BRS' and replace_strain_name(gene[0]) == leaf.replace('_2', ''):
                                #     border = 'black'
                                if 'A1401_12770' in gene[0]:
                                    gene[1] = int(gene[1]) + 11
                                    
                                gene[0] = 'BRS'
                                    
                            format_str = f'Arial|14|white|{gene[0]}' #Text on the gene
                            if gene[0] == 'GS2_BRS': #I use a gradient for multi-GH70 domain proteins.
                                col = f'{"rgradient:" + gene_colors["BRS"] + "|" + gene_colors["GS2"]}'
                            else:
                                col = gene_colors[gene[0]]
                            if border == 'white':
                                border = col[-7:]
                                
                            info_list = [int(gene[1])-start, int(gene[2])-start, '[]', 0, 20, border, col, format_str]
                            
                            if leaf not in gene_dict.keys():
                                gene_dict[leaf] = [info_list]
                            else:
                                gene_dict[leaf].append(info_list)
        
    length_dict = {key:20000 for key in gene_dict} #gene_dict[key][-1][1]+100 for key in gene_dict}
    
    [fix_strand(gene_dict[key]) for key in gene_dict if 'MP2' not in key]
    
    seq_dict = {l: f'{"-"*(gapscale)+"A"*(int(length_dict[l])-gapscale+padding-100)}' for l in lnames} #Create a sequence of desired length
    
    gene_domains = {}

    for key in gene_dict.keys():
        print(key)
        new_domains = []
        # add_BRS = True
        for element in gene_dict[key]:
            divide = False
        #     if (GH == 'BRS' and key in GS2_BRS and 'GS2_BRS' in element[7]) or (GH == 'BRS' and key not in GS2_BRS):
        #         add_BRS = True
        #         print(element[7], GH)
        #     else:
        #         add_BRS = False
            if GH == 'BRS' and key in GS2_BRS and GH in element[7] and element == gene_dict[key][-1]:
                divide = True
                print(key, GH, 'True1')
            elif GH == 'BRS' and f'|{GH}' in element[7] and 'GS2_BRS' not in gene_dict[key][-1][7] and element[1] - element[0] > 2000:
                divide = True
                print(key, GH, 'True2')
            elif GH == 'GS2' and f'|{GH}' in element[7]:
                divide = True
            elif GH == 'GS1' and f'|{GH}' in element[7]:
                check = key in max_GS1.values()
                if element == gene_dict[key][check] or 'MP2' in key:
                    divide = True
                print(key, GH, 'True3')
                
            if divide:
                new_domain = element.copy()
                new_domain[0] = element[0] + domain_dict[key][0] - 1
                new_domain[1] = element[0] + domain_dict[key][1] - 1
                new_domain[7] = new_domain[7].replace('GS2_BRS', GH)
                # if f'{key}_2' in GS2_BRS or key.replace('_2', '') in GS2_BRS or key in GS2_BRS:
                #     new_domain[7] = new_domain[7].replace('GS2_BRS', GH)
                #     print(GH)
                # element[2] = '[]'
                for GH_typ in GH_types:
                    if GH_typ in element[7] and 'GS2_BRS' not in element[7]:
                        element[5] = alt_colors[GH_typ]
                        element[6] = alt_colors[GH_typ]
                    else:
                        element[5] = element[5].replace(gene_colors[GH_typ], alt_colors[GH_typ])
                        element[6] = element[6].replace(gene_colors[GH_typ], alt_colors[GH_typ])
                element[7] = element[7].replace('GS2_BRS', '').replace(gene_types[key], '')
                next_domain = element.copy()
                next_domain[0] = new_domain[1] + 1
                element[1] = new_domain[0] - 1
                new_domains.append(element)
                new_domains.append(new_domain)
                new_domains.append(next_domain)
            else:
                # element[2] = '[]'
                for GH_typ2 in GH_types:
                    if GH_typ2 in element[7] and 'GS2_BRS' not in element[7]:
                        element[5] = alt_colors[GH_typ2]
                        element[6] = alt_colors[GH_typ2]
                    elif GH_typ2 in element[7]:
                        element[5] = element[5].replace(gene_colors[GH_typ2], alt_colors[GH_typ2])
                        element[6] = element[6].replace(gene_colors[GH_typ2], alt_colors[GH_typ2])
                new_domains.append(element)
        gene_dict[key] = new_domains # I need to add some kind of check where I check that the key gene type is the same as the gene type to add domains, and then I need to check that it is the right domain when there are two genes that are classified as the same
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
    #t.convert_to_ultrametric() #For nicer visualization
    t.render(f'{outdir}/{outplot}', tree_style = ts) 
