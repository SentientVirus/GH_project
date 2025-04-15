#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script plots the strain phylogeny considering only the representative 38
strains, either because they are not in the original tree or because they have
been pruned out of the final tree.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================

from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import os

# =============================================================================
# 1. Define input paths and formatting variables
# =============================================================================

#Dictionary to color strains by phylogroup
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

#Input treefiles to be plotted
treefile1 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.treefile' 
treefile2 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.nofilter.treefile' 
treefile3 = os.path.expanduser('~') + '/GH_project/all_core/trees/all_strains.treefile'
treefiles = [treefile1, treefile2, treefile3]

#Output directory
outdir = os.path.expanduser('~') + '/GH_project/plots/trees'

#Create output directory if it does not exist
if not os.path.exists(outdir):
   os.makedirs(outdir)
   
# =============================================================================
# 2. Loop through tree files and generate output plot
# =============================================================================

for treefile in treefiles: #Loop through tree files
    outfile = f'{outdir}/{treefile.split("/")[-1].split(".treefile")[0]}_fullset.png' #Define path to the output plot
    print(outfile) #Print path to output

    t = Tree(treefile, format = 0) #Create a tree object from the tree file
    
    if treefile != treefile3:
        outnode = t.search_nodes(name = 'fhon13')[0] #Set the outgroup to the phylogroup F strain A1001
    else:
        outnode = t.search_nodes(name = 'A1001')[0]
    t.set_outgroup(outnode) #Root the tree
    
    ts = TreeStyle() #Create tree style object
    ts.show_branch_length = False #Hide branch lengths
    ts.show_branch_support = False #Hide branch supports to add formatted text
    ts.show_leaf_name = False #Hide leaf names to add formatted text
    ts.scale =  30000 #Scale of the tree
    ts.scale_length = 0.01 #Scale legend bar
    
    ns = NodeStyle() #Create node style
    ns['size'] = 0 #Hide nodes
    ns['vt_line_width'] = 5 #Set width of vertical lines
    ns['hz_line_width'] = 5 #Set width of horizontal lines
    ns['hz_line_type'] = 0 #Horizontal lines will be solid lines
    
    for n in t.traverse(): #Loop through nodes in the tree
       n.set_style(ns) #Apply the style to each node
       if n not in t.get_leaves() and n.support >= 50: #If the node is not a leaf and support bigger or equal than 50%
           if n.support >= 95: #If the node support is bigger or equal than 95%
               color = 'black' #Color the support in black
           else: color = 'dimgrey' #Otherwise, color the support in grey
           support_face = TextFace(int(n.support), fgcolor = color, #Create text for support values and set its color
                                   ftype = 'Arial', fsize = 20) #Sent font type and size
           n.add_face(support_face, column = 0, position = 'branch-top') #Add the text to the node
       
    for leaf in t.get_leaves(): #Loop through the leaves sorted by phylogeny 
        nleaf = leaf.name.replace('-', '').upper()  #Create leaf names without the '-' symbol
        color = leaf_color.get(nleaf, None) #Use that to retrieve the leaf color
        name_face = TextFace(leaf.name, ftype = 'Arial', fgcolor = color, 
                             fsize = 40) #Create the text with leaf names
        leaf.add_face(name_face, column = 0, position = 'branch-right') #Add the text to the leaf
        
    t.ladderize(1) #Reorder the nodes to keep the order consistent in all plots
    nodeA = t.get_common_ancestor('H3B1-01A', 'G0101')
    nodeA.ladderize(0)
    nodeB = t.get_common_ancestor('H1B1-04J', 'G0101')
    nodeB.ladderize(1)
    nodeC = t.get_common_ancestor('G0101', 'H4B2-05J')
    nodeC.ladderize(0)
    nodeD = t.get_common_ancestor('fhon2', 'H3B2-02X')
    nodeD.ladderize(1)
    nodeE = t.get_common_ancestor('H3B2-02X', 'H1B3-02M')
    nodeE.ladderize(1)
    nodeF = t.get_common_ancestor('DSM', 'H3B1-04J')
    nodeF.ladderize(1)
    nodeG = t.get_common_ancestor('A1202', 'A1401')
    nodeG.ladderize(1)
    nodeH = t.get_common_ancestor('H4B2-06J', 'H3B2-06M')
    nodeH.ladderize(0)
    nodeI = t.get_common_ancestor('H3B2-03M', 'IBH001')
    nodeI.ladderize(1)
    nodeJ = t.get_common_ancestor('A0901', 'IBH001')
    nodeJ.ladderize(0)
    nodeK = t.get_common_ancestor('H4B4-12M', 'H4B5-01J')
    nodeK.ladderize(1)
    
    t.render(outfile, tree_style = ts) #Save the tree plot to a PNG file
    t.render(outfile.replace('png', 'tiff'), tree_style = ts) #Save plot to TIFF
    t.render(outfile.replace('png', 'svg'), tree_style = ts) #Save plot to SVG
