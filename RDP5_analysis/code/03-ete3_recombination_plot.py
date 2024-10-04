#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

@author: marmo435
"""

from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace, RectFace
# from TreeStyle import RectFace
# import colorsys
# import random
# from Bio import SeqIO
import os
# import re, csv
# from matplotlib.colors import to_hex
# from math import sqrt
import pandas as pd
from matplotlib import pyplot as plt

#Update to show also major parent and fix support values
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

to_include = ['A0901', 'A1001', 'A1003', 'A1202', 'A1401', 'A1805',
              'DSMZ12361', 'fhon2', 'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A', 
              'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 'H3B2-02X',
              'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 'H4B1-11J', 
              'H4B2-02J', 'H4B2-04J', 'H4B2-05J', 'H4B2-06J', 'H4B2-11M', 
              'H4B4-02J', 'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 
              'H4B5-03X', 'H4B5-04J', 'H4B5-05J', 'IBH001', 'MP2']

# OBS! Prune the non-closed phylogeny with ete3!
aln_len = 31420
prefix = 'subset'
start = 1 #Chosen so that reference gene is included and distance is 40000
end = start + aln_len
treefile = os.path.expanduser('~') + '/GH_project/trees/kunkeei_nonclosed.tree'
outplot = os.path.expanduser('~') + '/GH_project/plots/trees/RDP5_tree.png'
recomb = os.path.expanduser('~') + '/GH_project/RDP5_analysis/RDP_output/all_subsets_7methods_5+_filtered.csv'
in_tab = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab/all_subsets_positions_aln.tab'


t = Tree(treefile, format = 1)

t.set_outgroup('A1001')

for node in t.traverse():
    if node.is_leaf() and node.name not in to_include:
        node.delete()
    elif node.name == 'A1001':
        root_node = node
    elif node.name == 'fhon2':
        node.name = 'Fhon2'
    elif node.name == 'H4B2-06J':
        IBH_sis = node
    elif node.name == 'H3B1-04J':
        DSM_sis = node

        
# This works now, but I should also add sister nodes where DSM and IBH001 should be placed.
root_node.add_sister(name = 'IBH001', dist = 0)
root_node.add_sister(name = 'DSMZ12361', dist = 0)
IBH_sis.add_child(name = '**', dist = 0)
IBH_sis.add_child(name = IBH_sis.name, dist = 0)
IBH_sis.name = ''
DSM_sis.add_child(name = '*', dist = 0)
DSM_sis.add_child(name = DSM_sis.name, dist = 0)
DSM_sis.name = ''

ts = TreeStyle()

ts.show_branch_length = False
#ts.show_branch_support = True
ts.scale = 4000 #General tree scale, original 4000
#ts.tree_width = 0.001
#ts.force_topology = True
ts.branch_vertical_margin = 10 #Space between branches
ts.show_leaf_name = False #Hide unformatted leaf names
ts.show_scale = False #Hide tree scale

ts.legend_position = 2 # Upper right
ts.legend.margin_right = -100
ts.legend.add_face(TextFace('p-value', fsize = 15), column = 0)
ts.legend.add_face(TextFace(' ', fsize = 12), column = 1)
ts.legend.add_face(TextFace('<1E-150', fsize = 12), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#944654', 
                            label = None), column = 1)
ts.legend.add_face(TextFace('<1E-100', fsize = 12), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#C46877', 
                            label = None), column = 1)
ts.legend.add_face(TextFace('<1E-50', fsize = 12), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#DB8484', 
                            label = None), column = 1)
ts.legend.add_face(TextFace('<1E-25', fsize = 12), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#F7B9B0', 
                            label = None), column = 1)
ts.legend.add_face(TextFace('<1E-10', fsize = 12), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#F7CFC1', 
                            label = None), column = 1)
ts.legend.add_face(TextFace('<0.05', fsize = 12), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = 'linen', 
                            label = None), column = 1)


#Formatting of nodes (and branch lines)
ns = NodeStyle()
ns['size'] = 0 #Nodes are not visible
#ns['hz_line_type'] = 1
ns['hz_line_color'] = 'lightgrey'
ns['hz_line_width'] = 4
ns['vt_line_color'] = 'slategrey'
ns['vt_line_width'] = 4

ns2 = NodeStyle()
ns2['size'] = 0
ns2['hz_line_color'] = 'white'
ns2['vt_line_color'] = 'white'

ns3 = NodeStyle()
ns3['size'] = 0
ns3['hz_line_color'] = 'lightgrey'
ns3['hz_line_width'] = 2
ns3['vt_line_color'] = 'slategrey'
ns3['vt_line_width'] = 2
ns3['hz_line_type'] = 2

invisible_nodes = [t.get_common_ancestor('A1001'), 
                   t.get_common_ancestor('IBH001', 'DSMZ12361'), 
                   t.get_common_ancestor('IBH001', 'A1003'), 
                   t.search_nodes(name = 'IBH001')[0], 
                   t.search_nodes(name = 'DSMZ12361')[0], 
                   t.search_nodes(name = 'A1001')[0]]

for n in t.traverse():
    if n not in invisible_nodes and '*' not in n.name:
        n.set_style(ns)
    elif '*' in n.name:
        n.set_style(ns3)
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
                support = [l[11], l[12], l[13], l[14], l[15], l[16], l[17]]
                while 'NS' in support:
                    support.remove('NS')
                support = [float(sup) for sup in support]
                if support != []:
                    support = sorted(support, reverse = False)[:3]
                    support = sum(support)/len(support)
                else: support = 1
                # print(start, end, int(start), int(end))
                new_tract = [l[8], l[9], l[10], start, end, support]
                recomb_tracts[recomb].append(new_tract)
                
gene_positions = []

loci = []
labels = {}
countvar = 0
with open(in_tab) as tabfile:
    tab_df = pd.read_csv(tabfile, sep = '\t')
    tab_df = tab_df.reindex(index = tab_df.index[::-1])
    for index, row in tab_df.iterrows():
        if row['strain'] == 'A0901':
            if row['strand'] == -1:
                shape = '>'
            else:
                shape = '<'
            labels[countvar] = row['gene_name']
            countvar += 1
            locus = [row['start']-1, row['end']-1, shape, 0, 40, 'slategrey', '#9DAAB7', f'Arial|12|black|{row["gene_name"]}']
            if row['gene_name'] == 'bcrA':
                locus[0] += 1
            loci.append(locus)

#Here I add the count info + plot genome fragment
segment_count = 0
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
                if float(tract[5]) < 1E-150:
                    color = '#944654' #'#f08080'
                    color2 = 'white'
                    check = True
                elif float(tract[5]) < 1E-100:
                    color = '#C46877' #'#f4978e'
                    color2 = 'white'
                    check = True
                elif float(tract[5]) < 1E-50:
                    color = '#DB8484' #'#f8ad9d'
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 1E-25:
                    color = '#F7B9B0' #'#fbc4ab'
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 1E-10:
                    color = '#F7CFC1' #'#ffdab9'
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 0.05:
                    color, color2 = 'linen', 'black'
                    check = True
                else: check = False
            if 'Unknown' in tract[1]:
                tract[1] = tract[1].replace('Unknown (', '').replace(')', '?')
            tract[1] = tract[1].replace('Z12361', '')
                    
            if check:
                info_list = [int(tract[3]), int(tract[4]), '[]', 0, 20, 'lightgrey', color, f'Arial|12|{color2}|{tract[1]}']
                motifs.append(info_list)
            
        tract_list = recomb_tracts[leaf.name]
        
        to_add = []
        for i in range(len(motifs)):
            if motifs[i][1] - motifs[i][0] > 500: #or motifs[i][0] == 15453 or motifs[i][1] == 15452:
                if len(to_add) > 0:
                    check = []
                    for motif in to_add:
                        if (motifs[i][1] < motif[0]) or (motifs[i][0] > motif[1]):
                            check.append(False)
                        else: check.append(True)
                    if sum(check) != 0:
                        seqFace = SeqMotifFace(seq, motifs = to_add, height = 20, seq_format = '[]', gap_format = 'blank', scale_factor = 0.04, bgcolor = 'lightgrey', fgcolor = 'lightgrey')
                        (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned')
                        text_face = TextFace('|', fsize = 15, fgcolor = 'yellow')
                        text_face.margin_top = -25.6
                        text_face.margin_left = 634
                        (t & f'{leaf.name}').add_face(text_face, 1, 'aligned')
                        segment_count += 1
                        to_add = [motifs[i]]
                    else: to_add.append(motifs[i])     
                else: to_add.append(motifs[i])
                
        if to_add != []:
            seqFace = SeqMotifFace(seq, motifs = to_add, height = 20, seq_format = '[]', gap_format = 'blank', scale_factor = 0.04, bgcolor = 'lightgrey', fgcolor = 'lightgrey') #Add presence/absence info to node
            (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #The number represents the column
            text_face = TextFace('|', fsize = 15, fgcolor = 'yellow')
            text_face.margin_top = -25.6
            text_face.margin_left = 634
            (t & f'{leaf.name}').add_face(text_face, 1, 'aligned')
            segment_count += 1
        elif leaf.name == 'A1001':
            seqFace = SeqMotifFace(seq, motifs = loci, height = 20, seq_format = 'blank', gap_format = 'blank', scale_factor = 0.04, bgcolor = 'lightgrey', fgcolor = 'lightgrey')
            (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #The number represents the column
            
        seqFace = SeqMotifFace(miniseq, height = 20, seq_format = 'blank', gap_format = 'blank') #Add strain number
        (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned') #The number represents the column
    motif_dict[leaf.name] = motifs
    
    
#Here I print the tree
t.ladderize(1)
t.convert_to_ultrametric() #For nicer visualization
t.render(outplot, tree_style = ts)
t.render(outplot.replace('png', 'svg'), tree_style = ts) 
t.render(outplot.replace('png', 'pdf'), tree_style = ts) 
