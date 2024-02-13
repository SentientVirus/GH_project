#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

Another version of the script to plot the deleted region. This version uses
the versions of the files from NCBI, and is therefore more updated.

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

# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

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
              'G0803': '#2C85D8', 'G0804': '#2C85D8', 'fhon2': '#2C85D8',
              'H1B1-04J': '#2C85D8', 'H1B1-05A': '#2C85D8', 'H1B3-02M': '#2C85D8',
              'H2B1-05J': '#2AA380', 'H3B1-01A': '#2C85D8', 'H3B1-01J': '#2C85D8',
              'H3B1-01X': '#2AA380', 'H3B1-02A': '#2C85D8', 'H3B1-02X': '#2AA380',
              'H3B1-03J': '#2C85D8', 'H3B1-03X': '#2AA380','H3B1-03M': '#2C85D8', 
              'H3B1-04J': '#2C85D8', 'H3B1-04X': '#2C85D8', 'H3B1-07A': '#2C85D8', 
              'H3B1-09M': '#2C85D8', 'H3B1-10M': '#2C85D8', 'H3B1-11A': '#2C85D8', 
              'H3B1-11M': '#2C85D8', 'H3B2-02M': '#2AA380', 'H3B2-02X': '#2C85D8', 
              'H3B2-03J': '#2C85D8', 'H3B2-03M': '#2AA380', 'H3B2-04J': '#2C85D8', 
              'H3B2-04M': '#2AA380', 'H3B2-05J': '#2C85D8', 'H3B2-06M': '#2AA380', 
              'H3B2-07X': '#2C85D8', 'H3B2-08X': '#2C85D8', 'H3B2-09X': '#21C1D1', 
              'H4B1-01A': '#2AA380', 'H4B1-02A': '#2AA380', 'H4B1-03J': '#2AA380', 
              'H4B1-04A': '#2AA380', 'H4B1-11J': '#2AA380', 'H4B1-14J': '#2AA380', 
              'H4B1-16J': '#2AA380', 'H4B2-02J': '#2C85D8', 'H4B2-03M': '#2AA380', 
              'H4B2-04J': '#2C85D8', 'H4B2-05J': '#2C85D8', 'H4B2-06J': '#2AA380', 
              'H4B2-10M': '#2AA380', 'H4B2-11M': '#2C85D8', 'H4B3-03J': '#2AA380', 
              'H4B4-02J': '#2C85D8', 'H4B4-03J': '#2AA380', 'H4B4-04J': '#2AA380', 
              'H4B4-05J': '#2AA380', 'H4B4-06M': '#2AA380', 'H4B4-10M': '#2AA380', 
              'H4B4-11M': '#2AA380', 'H4B4-12M': '#2C85D8', 'H4B5-01J': '#2C85D8', 
              'H4B5-02X': '#2C85D8', 'H4B5-03X': '#2C85D8', 'H4B5-04J': '#21C1D1', #'#2AA380', 
              'H4B5-05J': '#21C1D1', 'H4B5-07J': '#2C85D8', 'H4B5-07X': '#2C85D8', 
              'H4B5-08X': '#2C85D8', 'MP2': '#21C1D1'}

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


# =============================================================================
# Set target region and load input files
# =============================================================================

start = 1270000 #Chosen so that reference gene is included and distance is 40000
end = 1301000 #Original 1301000
scale = 200 #Width for the presence/absence plot
tabs = snakemake.input.tabs
comp_dir = os.path.commonpath(tabs)
treefile = snakemake.input.tree #'kunkeei_nonclosed.tree'
count_file = snakemake.input.counts
outdir = 'plots/trees'
if not os.path.exists(outdir):
    os.makedirs(outdir)
outplot = 'phylogeny.png'

# =============================================================================
# Read types from Snakemake
# =============================================================================

GS1 = snakemake.params.GS1
GS2 = snakemake.params.GS2
BRS = snakemake.params.BRS
short = snakemake.params.short
NGB = snakemake.params.NGB
S1 = snakemake.params.S1
S2a = snakemake.params.S2a
S2b = snakemake.params.S2b
S3 = snakemake.params.S3

# =============================================================================
# Create rooted tree from treefile
# =============================================================================
t = Tree(treefile, format = 1) #Create tree

#Root tree
# outnode = t.get_common_ancestor('A2002', 'A1001')
# t.set_outgroup(outnode)

# =============================================================================
# Set style of tree
# =============================================================================
ts = TreeStyle() #Create default tree style
ts.show_branch_length = False #Hide support values
ts.scale =  5000 #General tree scale
#ts.tree_width = 200
#ts.force_topology = True
ts.branch_vertical_margin = 5 #Space between branches
ts.show_leaf_name = False #Hide unformatted leaf names
ts.show_scale = False #Hide tree scale

#Add legend
ts.legend_position = 3 #Bottom left
ts.legend.add_face(TextFace(' '*58, fsize = 10), column = 0)

#Add gene name text at the botton of the tree
names = ['GH70', 'GS1', 'GS2', 'BRS', 'NCB', 'Short', ' ', 'GH32', 'S1', 'S2', 
         'S3']
for n in range(len(names)):
    if names[n] in names[8:10]:
        text_face = TextFace(f'{names[n]}  ', fsize=15, tight_text = True)
    else:
        text_face = TextFace(f'{names[n]}', fsize=15, tight_text = True)
    text_face.rotation = 290
    text_face.margin_bottom = -9 #Reduce space between labels
    ts.legend.add_face(text_face, column = n+1)


#Formatting of nodes (and branch lines)
ns = NodeStyle()
ns['size'] = 0 #Nodes are not visible
for n in t.traverse():
   n.set_style(ns)

# =============================================================================
# Save CDS information for the main plot
# =============================================================================
#Get strain names
leaves = t.get_leaves() #Sorted by phylogeny
lnames = [leaf.name for leaf in leaves]
lnames.reverse()

gene_dict = {}
seq_dict = {l: f'{"A"*int(end-start)}' for l in lnames} #Create a sequence of desired length

#Add leaf names and plot region
for leaf in lnames:
    n = 0
    for filename in sorted(os.listdir(comp_dir)):
        if leaf == 'fhon2':
            leaf = 'Fhon2'
        if leaf in filename and filename.endswith('gpr.tab'):
            if leaf == 'MP2':
                leaf = '55'
            print(leaf)
            check = False
            with open(f'{comp_dir}/{filename}', 'r') as gfile:
                genes = csv.reader(gfile, delimiter='\t')
                for gene in genes:
                    if 'name' not in gene: #Skip first line
                        if leaf != '55' and gene[0] == 'ohrR' and not check: #Gene used to center the plotted fragments in the right position
                            start = int(gene[1])+5000
                            if leaf == 'H3B2-03M':
                                start += 1400
                            elif leaf == 'H1B1-05A':
                                start += 1300
                            elif leaf == 'H1B3-02M':
                                start += 1450
                            end = start + 30500 #Original 30500
                            check = True
                        if leaf == '55':
                            start = 846400
                            end = start + 30500 #Original 30500
                            check = True
                            #if leaf == 'G0407':
                            #    print((start, end))
                    if check and 'name' not in gene and end >= int(gene[1]) >= start and ((leaf != '55' and gene[3] == '-') or (leaf == '55' and gene[3] == '+')): #If the gene is in the right range and strand.
                        if gene[0][3:] in GS1 or (gene[0] == 'gtfC' and leaf != 'G0417'):
                            gene[0] = 'GS1'
                        elif gene[0][3:] in GS2:
                            gene[0] = 'GS2'
                        elif f'{gene[0][3:]}_2' in GS2 and gene[0][3:] in BRS:
                            gene[0] = 'GS2_BRS'
                        elif gene[0][3:] in BRS:
                            gene[0] = 'BRS'
                        elif gene[0][3:] in S1:
                            gene[0] = 'S1'
                        elif gene[0][3:] in S2a:
                            gene[0] = 'S2a'
                        elif gene[0][3:] in S2b:
                            gene[0] = 'S2b'
                        elif gene[0][3:] in S3:
                            gene[0] = 'S3 '
                        elif gene[0] == 'gtfC' and leaf == 'G0417':
                            gene[0] = 'S1'
                        elif gene[0] == 'GH39' or 'AKU' in gene[0] or 'APS55' in gene[0]:
                            gene[0] = 'hpr'
                        elif gene[0] == 'tesE': #Wrongly annotated by Prokka, here I correct that.
                            gene[0] = 'mhpD'
                        print(f'{leaf}_leaf', gene[0])
                        format_str = f'Arial|14|white|{gene[0]}' #Text on the gene
                        if gene[0][0] == 'S':
                            format_str = f'Arial|14|black|{gene[0]}'
                        if gene[0] != 'hpr' and gene[0][0:2] != 'S1':
                            col = gene_colors[gene[0][0:3]]
                            if gene[0] == 'GS2_BRS': #I use a gradient for multi-GH70 domain proteins.
                                col = f'{"rgradient:" + gene_colors["BRS"] + "|" + gene_colors["GS2"]}'
                        elif gene[0][0:2] == 'S1':
                            col = gene_colors['S1 ']
                        else:
                            if leaf != '55':
                                col = hpr_colors[n % 9] #Hypothetical proteins in greyish colors.
                            else: col = hpr_colors[::-1][n % 9]
                            n += 1
                            
                        if leaf != '55':
                            info_list = [int(gene[1])-start, int(gene[2])-start, '[]', 0, 20, col[-7:], col, format_str]
                        else:
                            info_list = [abs(int(gene[2])-end), abs(int(gene[1])-end), '[]', 0, 20, col[-7:], col, format_str]
                        
                        new_name0 = leaf.replace('-', '')
                        if new_name0 not in gene_dict.keys():
                            gene_dict[new_name0] = [info_list]
                        else:
                            gene_dict[new_name0].append(info_list)

                #gene_dict[new_name0].sort(key=lambda x: x[0])

# =============================================================================
# Info for the presence/absence plot           
# =============================================================================
comp_dir = '/'.join(count_file.split('/')[:-1])

presence_dict = {}
val_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

#Here I get the count info for GH70 domains and save it into a dictionary
with open(count_file) as presence:
    gtf_pres = csv.reader(presence, delimiter = '\t')
    for line in gtf_pres:
        if line[0] != 'strain':
            if line[0] == 'MP2':
                line[0] = '55'
            presence_dict[line[0]] = [line[6], line[1], line[2], line[3], 
                                      line[4], line[5], line[10], line[7], 
                                      line[8], line[9]]
            val_list = [max(val_list[i], int(presence_dict[line[0]][i])) for i in range(len(val_list))]

# =============================================================================
# Create the plot
# =============================================================================
for leaf in leaves:
    new_name = leaf.name.replace('-', '')
    if new_name == 'fhon2':
        new_name = 'Fhon2'
    elif new_name == 'MP2':
        new_name = '55'
    color = leaf_color.get(leaf.name, None) #Color leaves according to strain phylogroup
    name_face = TextFace(leaf.name, fgcolor = color, fsize = 15)
    leaf.add_face(name_face, column=0, position='branch-right') #Add formatted leaf names
    if new_name in presence_dict.keys():
        seq = '-'*scale + 140*'-' #Create a gap-only sequence that can be hidden
        num_list = presence_dict[new_name]
        motif = [int(1/20*scale), scale, 'o', 0, 21, '#DADADA', '#DADADA', None] #Define position of grey circles
        motifs = [motif.copy() for _ in range(10)]
        #motifs.append([161, 180, 'o', 0, 20, 'white', 'white', None])
        colors = [leaf_color[leaf.name], '#FF0000', '#FF009B', '#D930BD', 
                  '#B600FF', '#7000FF', leaf_color[leaf.name], '#FF9600', 
                  '#FFC800', '#FFEF00']
        for i in range(len(num_list)): #Color circles depending on presence/absence of GH70 type + add count
            motifs[i][0] += int(3/20*scale)*i
            motifs[i][1] = motifs[i][0] + int(2/20*scale)
            if i >= 6:
                motifs[i][0] += int(2/20*scale)
                motifs[i][1] += int(2/20*scale)
            if int(num_list[i]) > 0:
                motifs[i][5] = colors[i]
                motifs[i][6] = colors[i]
                motifs[i][7] = f'Arial|12|white|{num_list[i]}'
                if 10 >= i > 6:
                    motifs[i][7] = f'Arial|12|black|{num_list[i]}'

    if leaf.name in seq_dict.keys():
        gene_list = gene_dict[new_name]
        seqFace = SeqMotifFace(seq, motifs = motifs, seq_format = 'line', gap_format = 'blank') #Add presence/absence info to node
        (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #The number represents the column
        seqFace = SeqMotifFace(seq_dict[leaf.name], motifs = gene_list, seq_format = 'line', gap_format = 'blank', scale_factor = 0.05) #Original 0.04, Add genomic segment info to node
        (t & f'{leaf.name}').add_face(seqFace, 2, 'aligned')

    
# =============================================================================
# Print tree
# =============================================================================
t.ladderize(1)
#node_A = t.get_common_ancestor('H4B5-03X', 'H1B3-02M')
#node_A.ladderize(1)
node_A = t.get_common_ancestor('A1802', 'H3B1-02X')
node_A.ladderize(0)
node_B = t.get_common_ancestor('H4B5-01J', 'H4B4-12M')
node_B.ladderize(1)
node_C = t.get_common_ancestor('A1802', 'H4B4-06M')
node_C.ladderize(1)
node_D = t.get_common_ancestor('H3B1-02X', 'H2B1-05J')
node_D.ladderize(1)
#node_C = t.get_common_ancestor('H4B5-04J', 'H4B4-04J')
#node_C.ladderize(0)
#node_D = t.get_common_ancestor('H4B2-06J', 'H3B1-01X')
#node_D.ladderize(1)
t.convert_to_ultrametric() #For nicer visualization
t.render(f'{outdir}/{outplot}', tree_style = ts) 
