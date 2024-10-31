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
import os

leaf_color = {'A0901': '#D55E00', 'A1001': '#771853', 'A1002': '#D55E00',
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



# Maybe I should re-do the treefile getting domains only after running Interproscan on all the genes in the outgroup file
treefile = os.path.expanduser('~') + '/GH_project/add_species/results/alignment/GH70_species.mafft.faa.treefile' 
outfile = os.path.expanduser('~') + '/GH_project/add_species/plots/GH70_multispecies_tree.png'
outdir = os.path.dirname(outfile)

if not os.path.exists(outdir):
   os.makedirs(outdir)


t = Tree(treefile, format = 0)
outnode = t.get_common_ancestor('AOR73699.1_4_3_gtf_L_fermentum', 'AAU08014.2_4_6_gtf_L_reuteri')
t.set_outgroup(outnode)

ts = TreeStyle()
#ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = False
ts.scale =  500 #1500
#ts.tree_width = 200
#ts.force_topology = True
ts.show_leaf_name = False
ts.scale_length = 0.2

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

                
# Fix this part for MP2. To fix MP2 successfully, I have to add the old locus tag from the Interproscan in the previous script
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
    if '4_' in nleaf:
        nleaf = nleaf.replace('4_', '4,')
    if nleaf[-2:] == '_2':
        mleaf = nleaf[:-2]
    else:
        mleaf = nleaf
    color = leaf_color.get(mleaf.split('_')[0], None)
    name_face = TextFace(nleaf, fgcolor = color, fsize = 40)
    leaf.add_face(name_face, column=0, position='branch-right')
    
    
#print(t)
t.ladderize(1)
# t.convert_to_ultrametric()
t.render(outfile, tree_style = ts)
t.render(outfile.replace('png', 'svg'), tree_style = ts)
t.render(outfile.replace('png', 'pdf'), tree_style = ts)
