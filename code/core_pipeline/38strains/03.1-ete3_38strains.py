#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates a plot of the GH70 phylogeny, including
outgroups.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 
# =============================================================================
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
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

to_include = ['A0901', 'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805',
             'DSM', 'fhon2', 'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A',
             'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 'H3B2-02X',
             'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 'H4B1-11J',
             'H4B2-02J', 'H4B2-04J', 'H4B2-05J', 'H4B2-06J', 'H4B2-11M',
             'H4B4-02J', 'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J',
             'H4B5-03X', 'H4B5-04J', 'H4B5-05J', 'IBH001', 'MP2']

# Maybe I should re-do the treefile getting domains only after running Interproscan on all the genes in the outgroup file
treefile1 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.treefile' 
treefile2 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.nofilter.treefile' 
treefile3 = os.path.expanduser('~') + '/GH_project/all_core/38_strains/trees/38strains.treefile'
treefiles = [treefile1, treefile2, treefile3]

outdir = os.path.expanduser('~') + '/GH_project/plots/trees'

if not os.path.exists(outdir):
   os.makedirs(outdir)

for treefile in treefiles:
    outfile = f'{outdir}/{treefile.split("/")[-1].split(".treefile")[0]}.png'
    print(outfile)

    t = Tree(treefile, format = 0)
    if treefile != treefile3:
        outnode = t.search_nodes(name = 'fhon13')[0] #t.get_common_ancestor('A1001', ['A2002', 'A2003']) #GH70
    else: outnode = t.search_nodes(name = 'A1001')[0]
    t.set_outgroup(outnode)
    to_keep = [node for node in t.traverse() if node.name in to_include] #Get all t>
    t.prune(to_keep) #Prune the tree
    
    ts = TreeStyle()
    #ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.scale =  2000
    #ts.tree_width = 200
    #ts.force_topology = True
    ts.show_leaf_name = False
    ts.scale_length = 0.05
    
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
           if n.support >= 50:
               support_face = TextFace(int(n.support), fgcolor = color, 
                                       ftype = 'Arial', fsize = 30)
               n.add_face(support_face, column=0, position='branch-top')
       
    leaves = t.get_leaves() #Sort by phylogeny
    lnames = [leaf.name for leaf in leaves]
    
                    
    # Fix this part for MP2. To fix MP2 successfully, I have to add the old locus tag from the Interproscan in the previous script
    for leaf in leaves:
        nleaf = leaf.name.replace('-', '').upper()
        color = leaf_color.get(nleaf, None)
        name_face = TextFace(leaf.name, ftype = 'Arial', fgcolor = color, fsize = 40)
        leaf.add_face(name_face, column=0, position='branch-right')
        
    #print(t)
    t.ladderize(1)
    # t.convert_to_ultrametric()
    t.render(outfile, tree_style = ts)
    t.render(outfile.replace('png', 'tiff'), tree_style = ts)
    t.render(outfile.replace('png', 'svg'), tree_style = ts)
