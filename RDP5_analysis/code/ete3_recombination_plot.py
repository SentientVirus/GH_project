#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

@author: marmo435
"""

from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace
import colorsys
import random
#Use PhyloTree
from Bio import SeqIO
import re, csv, os
from matplotlib.colors import to_hex
from math import sqrt
import pandas as pd

class recomb_obj:
    def __init__(self, recombinant, minor, major, start, end, RDPvalue):
        self.recombinant = recombinant
        self.minor = minor
        self.major = major
        self.start = start
        self.end = end
        self.pvalue = RDPvalue
    def __str__(self):
        return f'<{self.recombinant}: {self.major} acquired a segment ({self.start}:{self.end}), p-value {self.pvalue} from {self.minor}>'

leaf_color = {'A0901': '#D55E00', 'A1001': 'white', 'A1002': '#D55E00',
              'A1003': '#0072B2', 'A1201': '#D55E00', 'A1202': '#33B18F',
              'A1401': '#33B18F', 'A1404': '#DA73B3', 'A1802': '#D55E00',
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
              'G0803': '#0072B2', 'G0804': '#0072B2', 'FHON2': '#0072B2',
              'H1B104J': '#0072B2', 'H1B105A': '#0072B2', 'H1B302M': '#0072B2',
              'H2B105J': '#D55E00', 'H3B101A': '#0072B2', 'H3B101J': '#0072B2',
              'H3B101X': '#D55E00', 'H3B102A': '#0072B2', 'H3B102X': '#D55E00',
              'H3B103J': '#0072B2', 'H3B103X': '#D55E00','H3B103M': '#0072B2', 
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
              'DSM': 'black'}

# leaf_color = {'A0901': '#00ABB0', 'A1001': '#E7957A', 'A1002': '#00ABB0',
#               'A1003': '#C50060', 'A1201': '#00ABB0', 'A1202': '#392A98',
#               'A1401': '#392A98', 'A1404': '#2085FF', 'A1802': '#00ABB0',
#               'A1803': '#00ABB0', 'A1805': '#9EB539', 'A2001': '#00ABB0',
#               'A2002': '#E7957A', 'A2003': '#E7957A', 'A2101': '#9EB539',
#               'A2102': '#00ABB0', 'A2103': '#C50060', 'G0101': '#C50060', 
#               'G0102': '#C50060', 'G0103': '#C50060', 'G0401': '#C50060',
#               'G0402': '#C50060', 'G0403': '#810E21', 'G0404': '#C50060',
#               'G0405': '#C50060', 'G0406': '#810E21', 'G0407': '#C50060',
#               'G0408': '#C50060', 'G0410': '#C50060', 'G0412': '#C50060',
#               'G0414': '#C50060', 'G0415': '#C50060', 'G0417': '#C50060',
#               'G0420': '#810E21', 'G0601': '#C50060', 'G0602': '#C50060',
#               'G0702': '#C50060', 'G0801': '#C50060', 'G0802': '#C50060',
#               'G0803': '#C50060', 'G0804': '#C50060', 'fhon2': '#C50060',
#               'H1B1-04J': '#C50060', 'H1B1-05A': '#C50060', 'H1B3-02M': '#C50060',
#               'H2B1-05J': '#FF5733', 'H3B1-01A': '#C50060', 'H3B1-01J': '#C50060',
#               'H3B1-01X': '#FF5733', 'H3B1-02A': '#C50060', 'H3B1-02X': '#FF5733',
#               'H3B1-03J': '#C50060', 'H3B1-03X': '#FF5733','H3B1-03M': '#C50060', 
#               'H3B1-04J': '#C50060', 'H3B1-04X': '#C50060', 'H3B1-07A': '#C50060', 
#               'H3B1-09M': '#C50060', 'H3B1-10M': '#C50060', 'H3B1-11A': '#C50060', 
#               'H3B1-11M': '#C50060', 'H3B2-02M': '#FF5733', 'H3B2-02X': '#C50060', 
#               'H3B2-03J': '#C50060', 'H3B2-03M': '#FF5733', 'H3B2-04J': '#C50060', 
#               'H3B2-04M': '#00ABB0', 'H3B2-05J': '#C50060', 'H3B2-06M': '#FF5733', 
#               'H3B2-07X': '#C50060', 'H3B2-08X': '#C50060', 'H3B2-09X': '#392A98', 
#               'H4B1-01A': '#FF5733', 'H4B1-02A': '#FF5733', 'H4B1-03J': '#FF5733', 
#               'H4B1-04A': '#FF5733', 'H4B1-11J': '#FF5733', 'H4B1-14J': '#FF5733', 
#               'H4B1-16J': '#FF5733', 'H4B2-02J': '#C50060', 'H4B2-03M': '#FF5733', 
#               'H4B2-04J': '#C50060', 'H4B2-05J': '#C50060', 'H4B2-06J': '#00ABB0', 
#               'H4B2-10M': '#FF5733', 'H4B2-11M': '#C50060', 'H4B3-03J': '#FF5733', 
#               'H4B4-02J': '#C50060', 'H4B4-03J': '#FF5733', 'H4B4-04J': '#FF5733', 
#               'H4B4-05J': '#FF5733', 'H4B4-06M': '#FF5733', 'H4B4-10M': '#FF5733', 
#               'H4B4-11M': '#FF5733', 'H4B4-12M': '#C50060', 'H4B5-01J': '#C50060', 
#               'H4B5-02X': '#C50060', 'H4B5-03X': '#C50060', 'H4B5-04J': '#9EB539', 
#               'H4B5-05J': '#392A98', 'H4B5-07J': '#C50060', 'H4B5-07X': '#C50060', 
#               'H4B5-08X': '#C50060', 'MP2': '#157C00'}

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

aln_len = 31420
prefix = 'subset'
start = 1 #Chosen so that reference gene is included and distance is 40000
end = start + aln_len
treefile = os.path.expanduser('~') + '/GH_project/trees/own_representatives.tree'
outplot = os.path.expanduser('~') + '/GH_project/plots/trees/RDP5_tree.png'
recomb = os.path.expanduser('~') + '/GH_project/RDP5_analysis/RDP_output/recomb_toplot.csv'
in_tab = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab/all_subsets_positions_aln.tab'
t = Tree(treefile, format = 9)

ts = TreeStyle()

t.set_outgroup('A1001')

ts.show_branch_length = False
#ts.show_branch_support = True
ts.scale =  10 #General tree scale, original 4000
#ts.tree_width = 0.001
#ts.force_topology = True
ts.branch_vertical_margin = 10 #Space between branches
ts.show_leaf_name = False #Hide unformatted leaf names
ts.show_scale = False #Hide tree scale

# ts.legend_position = 3 #Bottom left
# ts.legend.add_face(TextFace(' '*96, fsize = 10), column = 0)
# names = ['GH70', 'GS1', 'GS2', 'BRS', 'NCB', 'Short', ' ', 'GH32', 'S1', 'S2', 'S3']
# for n in range(len(names)):
#     if names[n] in names[8:11]:
#         text_face = TextFace(f'{names[n]}  ', fsize=20, tight_text = True)
#     else:
#         text_face = TextFace(f'{names[n]}', fsize=20, tight_text = True)
#     text_face.rotation = 290
#     text_face.margin_bottom = -10 #Reduce space between labels
#     ts.legend.add_face(text_face, column = n+1)


#Formatting of nodes (and branch lines)
ns = NodeStyle()
ns['size'] = 0 #Nodes are not visible
#ns['hz_line_type'] = 1
ns['hz_line_color'] = 'lightgrey'
ns['vt_line_color'] = 'black'

ns2 = NodeStyle()
ns2['size'] = 0
ns2['hz_line_color'] = 'white'
ns2['vt_line_color'] = 'white'

invisible_nodes = [t.get_common_ancestor('A1001'), 
                   t.get_common_ancestor('IBH001', 'DSMZ12361'), 
                   t.get_common_ancestor('IBH001', 'A1003'), 
                   t.search_nodes(name = 'IBH001')[0], 
                   t.search_nodes(name = 'DSMZ12361')[0], 
                   t.search_nodes(name = 'A1001')[0]]

for n in t.traverse():
    if n not in invisible_nodes:
        n.set_style(ns)
    else: 
        n.set_style(ns2)
        n.img_style['hz_line_color'] = 'white'
        n.img_style['vt_line_color'] = 'white'

#Get strain names
leaves = t.get_leaves() #Sorted by phylogeny
lnames = [leaf.name for leaf in leaves]
lnames_sort = sorted(lnames)
lno = {lname: lnames_sort.index(lname)+1 for lname in lnames_sort}
lnames.reverse()

recomb_tracts = {leaf.name: [] for leaf in leaves}
#with open('GS1_S2-3_70kb_poster.csv') as rectracts:
with open(recomb) as rectracts:
    tracts = (line for line in rectracts)
    for i in range(16):
        next(tracts)
    for l in tracts:
        if '=' not in l and len(l.split(',')) > 1:
            l = l.split(',')
            if (l[2] and l[3]) != '':
                recomb = l[8].strip('^')
                start = ''.join(c for c in l[2] if c.isdigit())
                end = ''.join(c for c in l[3] if c.isdigit())
                print(start, end, int(start), int(end))
                new_tract = [l[8], l[9], l[10], start, end, l[11]]
                recomb_tracts[recomb].append(new_tract)
                
gene_positions = []

loci = []
with open(in_tab) as tabfile:
    tab_df = pd.read_csv(tabfile, sep = '\t')
    for index, row in tab_df.iterrows():
        if row['strain'] == 'A0901':
            locus = [row['start'], row['end'], '[]', 0, 20, 'slategrey', 'white', f'Arial|12|black|{row["gene_name"]}']
            if locus[0] == 915:
                locus[0] += 2
            loci.append(locus)
    
                
# gene_dict = {}
# seq_dict = {l: f'{"A"*int(end-start)}' for l in lnames} #Create a sequence of desired length

#In this section, I add leaf names
# for leaf in lnames:
#     rectract = recomb_tracts[leaf]
    # n = 0
    # for filename in sorted(os.listdir(comp_dir)):
    #     if leaf in filename and filename.endswith('gpr.tab'):
    #         check = False
    #         with open(comp_dir + filename, 'r') as gfile:
    #             genes = csv.reader(gfile, delimiter='\t')
    #             for gene in genes:
    #                 if 'name' not in gene: #Skip first line
    #                     if gene[0] == 'ohrR' and not check: #Gene used to center the plotted fragments in the right position
    #                         start = int(gene[1])+5000
    #                         if leaf == 'MP2' or leaf == 'H3B2-03M':
    #                             start += 1400
    #                         elif leaf == 'H1B1-05A':
    #                             start += 1300
    #                         elif leaf == 'H1B3-02M':
    #                             start += 1450
    #                         end = start + 30500 #Original 30500
    #                         check = True
    #                         #if leaf == 'G0407':
    #                         #    print((start, end))
    #                 if check and 'name' not in gene and end >= int(gene[1]) >= start and gene[3] == '-': #If the gene is in the right range and strand.
    #                     if gene[0] == 'out':
    #                         gene[0] = 'hpr'
    #                     elif gene[0] == 'gtf2a':
    #                         gene[0] = 'S2a'
    #                     elif gene[0] == 'gtf2b':
    #                         gene[0] = 'S2b'
    #                     elif gene[0] == 'gtf3':
    #                         gene[0] = 'S3 '
    #                     elif gene[0] == 'tesE': #Wrongly annotated by Prokka, here I correct that.
    #                         gene[0] = 'mhpD'
    #                     elif gene[0] == 'cap8A':
    #                         gene[0] = 'ywqC'
    #                     elif gene[0][0:3] == 'glu': #I change the annotation to make it consistent
    #                         gene[0] = gene[0].replace(gene[0][0:3], 'GS')
    #                         if len(gene[0]) == 8:
    #                             gene[0] = gene[0].replace(gene[0][4:8], 'BRS')
    #                     elif gene[0][0:3] == 'brs':
    #                         gene[0] = gene[0].replace(gene[0][0:4], 'BRS')
    #                     elif gene[0][0:4] == 'gtf1':
    #                         gene[0] = gene[0].replace('gtf', 'S')
    #                         gene[0] = gene[0] + ' '
    #                     format_str = f'Arial|14|white|{gene[0]}' #Text on the gene
    #                     if gene[0][0] == 'S':
    #                         format_str = f'Arial|14|black|{gene[0]}'
    #                     if gene[0] != 'hpr' and gene[0][0:2] != 'S1':
    #                         col = gene_colors[gene[0][0:3]]
    #                         if gene[0] == 'GS2_BRS': #I use a gradient for multi-GH70 domain proteins.
    #                             col = f'{"rgradient:" + gene_colors["BRS"] + "|" + gene_colors["GS2"]}'
    #                     elif gene[0][0:2] == 'S1':
    #                         col = gene_colors['S1 ']
    #                     else:
    #                         col = hpr_colors[n % 9] #Hypothetical proteins in greyish colors.
    #                         n += 1
    #                     info_list = [int(gene[1])-start, int(gene[2])-start, '[]', 0, 20, col[-7:], col, format_str]
    #                     if leaf not in gene_dict.keys():
    #                         gene_dict[leaf] = [info_list]
    #                     else:
    #                         gene_dict[leaf].append(info_list)

#I add presence/absence info for each gene type            
# comp_dir = '../../PhD/precipitation'

# presence_dict = {}
# val_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# #Here I get the count info for GH70 domains and save it into a dictionary
# with open(f'{comp_dir}/strain_GH70_32.tab') as presence: #f'{comp_dir}/strain_gtfs_choline_weighted.tab') as presence:
#     gtf_pres = csv.reader(presence, delimiter = '\t')
#     for line in gtf_pres:
#         if line[0] != 'strain':
#             presence_dict[line[0]] = [line[6], line[1], line[2], line[3], 
#                                       line[4], line[5], line[10], line[7], 
#                                       line[8], line[9]]
#             val_list = [max(val_list[i], int(presence_dict[line[0]][i])) for i in range(len(val_list))]

#Here I add the count info + plot genome fragment
motif_dict = {}
seq = 'A'*aln_len 
miniseq = 'A'*10
for leaf in leaves:
    color = leaf_color.get(leaf.name.replace('-', '').upper(), None) #Color leaves according to strain phylogroup
    if 'DSM' in leaf.name:
        leaf_name = '*DSM'
    elif 'IBH001' in leaf.name:
        leaf_name = '**' + leaf.name
    else: leaf_name = leaf.name
    name_face = TextFace(leaf_name, fgcolor = color, fsize = 16)
    leaf.add_face(name_face, column=0, position='branch-right') #Add formatted leaf names
    if leaf.name in recomb_tracts.keys():
        recomb_tract = recomb_tracts[leaf.name]
        motifs = []
        for tract in recomb_tract:
            check = False
            if tract[5] != 'NS':
                if float(tract[5]) < 1E-200:
                    color = 'white'
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 1E-150:
                    color = 'yellow'
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 1E-100:
                    color = 'orange'
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 1E-50:
                    color = 'red'
                    color2 = 'white'
                    check = True
                elif float(tract[5]) < 1E-25:
                    color = 'darkred'
                    color2 = 'white'
                    check = False
                else: color, color2 = 'black', 'white'
            #Find the way to avoid superposition of tracts
            #tract[1] = tract[1].replace('Z12361', '').replace('-', '').replace('IBH001', 'IBH')
            if 'Unknown' in tract[1]:
                tract[1] = tract[1].replace('Unknown (', '').replace(')', '?')
            tract[1] = tract[1].replace('Z12361', '')
                
            # for lname in lnames:
            #     if lname in tract[1]:
            #         tract[1] = tract[1].replace(lname, str(lno[lname]))
                    
            if check:
                info_list = [int(tract[3]), int(tract[4]), '[]', 0, 20, 'slategrey', color, f'Arial|12|{color2}|{tract[1]}']
                motifs.append(info_list)
            
        tract_list = recomb_tracts[leaf.name]
        
        # leaf_annot = [[1, 100, '[]', 0, 20, 'white', 'white', f'Arial|14|black|{lno[leaf.name]}']]
        
        to_add = []
        for i in range(len(motifs)):
            if motifs[i][1] - motifs[i][0] > 1000 or motifs[i][0] == 15453 or motifs[i][1] == 15452:
                if leaf.name == 'H4B2-06J':
                    seqFace = SeqMotifFace(seq, motifs = [motifs[i]], height = 20, seq_format = '[]', gap_format = 'blank', scale_factor = 0.05) #Add presence/absence info to node
                    (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #The number represents the column
                else: to_add.append(motifs[i])
            
        if leaf.name != 'H4B2-06J' and to_add != []:
            seqFace = SeqMotifFace(seq, motifs = to_add, height = 20, seq_format = '[]', gap_format = 'blank', scale_factor = 0.05) #Add presence/absence info to node
            (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #The number represents the column
        elif leaf.name == 'A1001':
            seqFace = SeqMotifFace(seq, motifs = loci, height = 20, seq_format = 'blank', gap_format = 'blank', scale_factor = 0.05)
            (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #The number represents the column
            
        seqFace = SeqMotifFace(miniseq, height = 20, seq_format = 'blank', gap_format = 'blank') #Add presence/absence info to node
        (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned') #The number represents the column
    motif_dict[leaf.name] = motifs
    
#Here I print the tree
t.ladderize(1)
t.convert_to_ultrametric() #For nicer visualization
t.render(outplot, tree_style = ts)
t.render(outplot.replace('png', 'svg'), tree_style = ts) 
t.render(outplot.replace('png', 'pdf'), tree_style = ts) 
