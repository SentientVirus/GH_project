#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates plots of the domains of GH70 and GH32 proteins
next to a strain phylogeny. To plot a specific gene type, the script has
to be manually updated.

@author: Marina Mota Merlo
"""

from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace
import re

leaf_color = {'A0901': '#A027FF', 'A1001': '#E5BA60', 'A1002': '#A027FF',
              'A1003': '#1E55F6', 'A1201': '#A027FF', 'A1202': '#00AEFF',
              'A1401': '#00AEFF', 'A1404': '#FF74D6', 'A1802': '#A027FF',
              'A1803': '#A027FF', 'A1805': '#00AEFF', 'A2001': '#A027FF',
              'A2002': '#E5BA60', 'A2003': '#E5BA60', 'A2101': '#00AEFF',
              'A2102': '#A027FF', 'A2103': '#1E55F6', 'G0101': '#1E55F6', 
              'G0102': '#1E55F6', 'G0103': '#1E55F6', 'G0401': '#1E55F6',
              'G0402': '#1E55F6', 'G0403': '#00AEFF', 'G0404': '#1E55F6',
              'G0405': '#1E55F6', 'G0406': '#00AEFF', 'G0407': '#1E55F6',
              'G0408': '#1E55F6', 'G0410': '#1E55F6', 'G0412': '#1E55F6',
              'G0414': '#1E55F6', 'G0415': '#1E55F6', 'G0417': '#1E55F6',
              'G0420': '#00AEFF', 'G0601': '#1E55F6', 'G0602': '#1E55F6',
              'G0702': '#1E55F6', 'G0801': '#1E55F6', 'G0802': '#1E55F6',
              'G0803': '#1E55F6', 'G0804': '#1E55F6', 'FHON2': '#1E55F6',
              'H1B104J': '#1E55F6', 'H1B105A': '#1E55F6', 'H1B302M': '#1E55F6',
              'H2B105J': '#A027FF', 'H3B101A': '#1E55F6', 'H3B101J': '#1E55F6',
              'H3B101X': '#A027FF', 'H3B102A': '#1E55F6', 'H3B102X': '#A027FF',
              'H3B103J': '#1E55F6', 'H3B103X': '#A027FF','H3B103M': '#1E55F6', 
              'H3B104J': '#1E55F6', 'H3B104X': '#1E55F6', 'H3B107A': '#1E55F6', 
              'H3B109M': '#1E55F6', 'H3B110M': '#1E55F6', 'H3B111A': '#1E55F6', 
              'H3B111M': '#1E55F6', 'H3B202M': '#A027FF', 'H3B202X': '#1E55F6', 
              'H3B203J': '#1E55F6', 'H3B203M': '#A027FF', 'H3B204J': '#1E55F6', 
              'H3B204M': '#A027FF', 'H3B205J': '#1E55F6', 'H3B206M': '#A027FF', 
              'H3B207X': '#1E55F6', 'H3B208X': '#1E55F6', 'H3B209X': '#00AEFF', 
              'H4B101A': '#A027FF', 'H4B102A': '#A027FF', 'H4B103J': '#A027FF', 
              'H4B104A': '#A027FF', 'H4B111J': '#A027FF', 'H4B114J': '#A027FF', 
              'H4B116J': '#A027FF', 'H4B202J': '#1E55F6', 'H4B203M': '#A027FF', 
              'H4B204J': '#1E55F6', 'H4B205J': '#1E55F6', 'H4B206J': '#A027FF', 
              'H4B210M': '#A027FF', 'H4B211M': '#1E55F6', 'H4B303J': '#A027FF', 
              'H4B402J': '#1E55F6', 'H4B403J': '#A027FF', 'H4B404J': '#A027FF', 
              'H4B405J': '#A027FF', 'H4B406M': '#A027FF', 'H4B410M': '#A027FF', 
              'H4B411M': '#A027FF', 'H4B412M': '#1E55F6', 'H4B501J': '#1E55F6', 
              'H4B502X': '#1E55F6', 'H4B503X': '#1E55F6', 'H4B504J': '#00AEFF',
              'H4B505J': '#00AEFF', 'H4B507J': '#1E55F6', 'H4B507X': '#1E55F6', 
              'H4B508X': '#1E55F6', 'MP2': '#00AEFF', 'IBH001': 'black', 
              'DSM': 'black'}



# shapes = {'RE': '[]', 'EL': '()', 'DI': '<>', 'TR': 'o'}

# domain_file = '../../PhD/Alignment/domains/domain_file.txt'
treefile = 'data/fasta/GH32/trees/GH32_repset.mafft.faa.treefile' #'data/fasta/GH32/trees/GH32_repset.mafft.faa.treefile' #'data/fasta/GH70/trees/GH70_functional_outgroup_repset.mafft.faa.treefile'
outfile = 'figure2.svg' #'figure2.png' #'figure1.png'
t = Tree(treefile, format = 0)

#outnode = t.get_common_ancestor('AAU08014.2_L_reuteri', 'AOR73699.1_L_fermentum') #GH70
outnode = t.get_common_ancestor('H3B203M_12520', 'A1404_13450') #GH32
#outnode = 'A1003_12540' #GS1
#outnode = t.get_common_ancestor('H4B5-05J_12900', 'H3B1-03M_13220') #GS2
#outnode = 'H3B1-01A_13250' #BRS
#outnode = t.get_common_ancestor('H3B1-02X_04940', 'H3B2-03M_04710') #NCB
#outnode = 'H3B1-01J_13900' #Short
t.set_outgroup(outnode)

ts = TreeStyle()
#ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = False
ts.scale =  500 #1500
#ts.tree_width = 200
#ts.force_topology = True
ts.show_leaf_name = False

ns = NodeStyle()
ns['size'] = 0
ns['vt_line_width'] = 5
ns['hz_line_width'] = 5
ns['hz_line_type'] = 0
for n in t.traverse():
   n.set_style(ns)
   if n not in t.get_leaves() and n.support > 1:
       if n.support >= 95:
           color = 'black' #'#0E9E00'
       else: color = 'dimgrey'
       # elif 90 <= n.support < 95:
       #     color = 'dimgrey'
       # elif 80 <= n.support < 90:
       #     color = 'grey' #'#5D9E00'
       # elif 70 <= n.support < 80:
       #     color = 'darkgrey' #'#809E00'
       # elif 60 <= n.support < 70:
       #     color = 'silver' #'#9E9E00'
       # elif 40 <= n.support < 60:
       #     color = 'lightgrey' #'#9E6D00'
       # elif 20 <= n.support < 40:
       #     color = '#9E4F00'
       # elif n.support < 20:
       #     color = '#9E0000'
       if n.support >= 50:
           support_face = TextFace(int(n.support), fgcolor = color, fsize = 24)
           n.add_face(support_face, column=0, position='branch-top')
   
leaves = t.get_leaves() #Sort by phylogeny
lnames = [leaf.name for leaf in leaves]

# with open(domain_file, 'r') as dfile:
#     domains = (d for d in dfile)
#     dom_dict = {}
#     seq_dict = {}
#     for domain in domains:
#         domain = domain.replace('\n', '')
#         domain = re.split(r'[|,]', domain)
#         if len(domain) >= 3:
#             seq_dict[domain[0]] = 'A'*int(domain[1])
#             dom_dict[domain[0]] = []
#             n = (len(domain)-1)
#             i = 2
#             while i < n: #Change signal peptide to blue circle
#                 if domain[i+4] == 'SP' or domain[i+4] == 'CWB':
#                     size = 15
#                     if domain[i+4] == 'CWB':
#                         domain[i+3] = '#585858'
#                         domain[i+4] = 'CB'
#                     else:
#                         domain[i] = 'TR'
#                         domain[i+3] = '#00AAFF'
#                 else:
#                     size = 24
#                     if domain[i+4] == 'GH70':
#                         domain[i+3] = '#FF0000' # #C70039 (all) #B600FF (NCB), #D930BD (BRS), #FF0000 (GS1), #FF009B (GS2), #7000FF (short)
#                     elif domain[i+4] == 'CB':
#                         domain[i+4] = 'GB'
#                         domain[i+3] = '#222222'
#                 if domain[i+4] != 'GH70':
#                     dom_dict[domain[0]].append([int(domain[i+1]), int(domain[i+2]), 
#                     shapes[domain[i]], None, 40, domain[i+3], domain[i+3], None])
#                 else:
#                     dom_dict[domain[0]].append([int(domain[i+1]), int(domain[i+2]), 
#                     shapes[domain[i]], None, 40, domain[i+3], domain[i+3],
#                     f'arial|{size}|white|{domain[i+4]}'])
                    
#                 i += 5
                
            
for leaf in leaves:
    # if leaf.name[-2:] == '_2':
    #     nleaf = leaf.name[:-2]
    # else:
    nleaf = leaf.name
    if 'LDX55' in nleaf:
        nleaf = nleaf.replace('LDX55', 'IBH001')
    elif 'APS55' in nleaf:
        nleaf = nleaf.replace('APS55_RS', 'MP2_')
    elif 'K2W83' in leaf.name:
        nleaf = nleaf.replace('K2W83_RS', 'DSM_')
    if nleaf[-2:] == '_2':
        mleaf = nleaf[:-2]
    else:
        mleaf = nleaf
    color = leaf_color.get(mleaf.split('_')[0], None)
    name_face = TextFace(nleaf, fgcolor = color, fsize = 40)
    leaf.add_face(name_face, column=0, position='branch-right')
    # if nleaf in seq_dict.keys():
    #     seqFace = SeqMotifFace(seq_dict[nleaf], motifs = dom_dict[nleaf], seq_format = 'line', scale_factor = 0.8)
    #     (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned')
    
    
#print(t)
t.ladderize(1)
# t.convert_to_ultrametric()
t.render(outfile, tree_style = ts)
