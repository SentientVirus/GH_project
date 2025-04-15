#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 18:37:42 2025

@author: marina
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates plots of the domains of GH70 and GH32 proteins
next to a strain phylogeny. To plot a specific gene type, the script has
to be manually updated.

@author: Marina Mota-Merlo
"""
##IN PROGRESS, now increase the scale of the trees, remove support values (but keep copies),
#change direction of the 38 strain tree and put everything together

# =============================================================================
# 0. Import required packages or modules
# =============================================================================
from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace
import re
import os

# =============================================================================
# 1. Define input variables and formatting variables
# =============================================================================

#Dictionary to assign colors to each locus tag depending on strain phylogroup
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
              'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': '#D55E00', 
              'DSM': '#0072B2'}

workdir = os.path.expanduser('~') + '/GH_project'
GH_dir = f'{workdir}/data/fasta/GH70/trees'
core_dir = f'{workdir}/all_core/38_strains/trees'
outdir = f'{workdir}/plots/trees/tanglegrams'

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it


for GH_type in ['GS1', 'NGB']:
    treefile1 = f'{GH_dir}/{GH_type}_repset.mafft.faa.treefile'
    treefile2 = f'{core_dir}/{GH_type}/38strains_{GH_type}.treefile'
    
    for treefile in [treefile1, treefile2]:
        if '38_strains' not in treefile:
            check = False
            outfile = f'{outdir}/{GH_type}_tree.png' #Output image with the tree representation
        else:
            check = True
            outfile = f'{outdir}/{GH_type}_38_strains_tree.png' #Output image with the tree representation
    
        # =============================================================================
        # 1. Load the tree file and set the outgroup
        # =============================================================================
        
        t = Tree(treefile, format = 0)  #Load the tree into a Tree object
        if GH_type == 'GS1' and not check:
            outnode = 'A1001_12310'
        elif GH_type == 'GS1' and check:
            outnode = 'A1001'
        elif GH_type == 'NGB' and not check:
            outnode = t.get_midpoint_outgroup()
        elif GH_type == 'NGB' and check:
            outnode = t.get_midpoint_outgroup()
         
        # if GH_type != 'NGB' or check:
        t.set_outgroup(outnode) #Root the tree on the outgroup node
        
        # =============================================================================
        # 2. Create and modify a tree style
        # =============================================================================
        
        ts = TreeStyle() #Create a tree style 
        ts.show_branch_length = False #Hide branch lengths
        ts.show_branch_support = False #Hide support values
        ts.show_leaf_name = False #Hide leaf names
        ts.scale = 5000 #Set the scale of the tree
        ts.scale_length = 0.01 #Set the length of the tree
        # ts.tree_width = 1 #Set the width of the branches, overwritten by the scale
        # ts.force_topology = True #Forces all the branches to be of equal length
        
        # =============================================================================
        # 3. Create and apply a node style
        # =============================================================================
        
        ns = NodeStyle() #Create a node style
        ns['size'] = 0 #Remove node representations as circles
        ns['vt_line_width'] = 1 #Set the width of vertical lines
        ns['hz_line_width'] = 1 #Set the width of horizontal lines
        ns['hz_line_type'] = 0 #Make horizontal lines solid
        
        for n in t.traverse(): #Loop through the nodes in the tree
           n.set_style(ns) #Apply the style to the nodes
           if n not in t.get_leaves() and n.support >= 50: #If the node is not a leaf and the support value is higher of equal to 50
               if n.support >= 95: #If the support value is > 95
                   color = 'black' #Color the support value in black
               else: color = 'dimgrey' #Otherwise, color it in grey
        
               # support_face = TextFace(int(n.support), fgcolor = color, fsize = 24,
               #                         ftype = 'Arial') #Create a text with the support value
               # n.add_face(support_face, column = 0, position='branch-top') #Add the text to the corresponding node in the tree
               
        # =============================================================================
        # 4. Modify the style of the leaf names in the tree
        # =============================================================================
        
        leaves = t.get_leaves() #Sort the leaves (locus tags) by phylogeny
        
        for leaf in leaves: #Loop through the leaves of the tree
            nleaf = leaf.name #Create an additional variable to store the leaf name
            if 'LDX55' in nleaf: #If LDX55 is in the locus tag
                nleaf = nleaf.replace('LDX55', 'IBH001') #Change it to IBH001
            elif 'APS55' in nleaf: #If APS55 is in the locus tag
                nleaf = nleaf.replace('APS55', 'MP2') #Change it to MP2
            elif 'K2W83' in leaf.name: #If K2W83 is in the locus tag
                nleaf = nleaf.replace('K2W83', 'DSM') #Change it to DSM
            
            mleaf = nleaf.split('_')[0]
            color = leaf_color.get(mleaf.replace('-', '').upper(), None) #Get strain names from the leaf name and use them to get the leaf color 
            if mleaf.startswith('H') and '-' not in mleaf:
                mleaf = mleaf[:4] + '-' + mleaf[4:]
            elif mleaf == 'FHON2':
                mleaf = 'fhon2'
            name_face = TextFace(mleaf, fgcolor = color, fsize = 40, ftype = 'Arial') #Create a text with locus tags
            leaf.add_face(name_face, column = 0, position = 'branch-right') #Add the text to the right leaf in the tree
            
        # =============================================================================
        # 5. Change final details of tree formatting and save the tree to a file
        # =============================================================================
        
        t.ladderize(1) #Reverse the order in which the leaves of the tree are plotted
        if '38_strains' in treefile:
            ts.orientation = 1
        # t.convert_to_ultrametric() #Convert the tree to a cladogram
        t.render(outfile, tree_style = ts) #Save the tree to a PNG image
        t.render(outfile.replace('png', 'svg'), tree_style = ts) #Save the tree to an SVG image
        t.render(outfile.replace('png', 'pdf'), tree_style = ts) #Save the tree to a PDF image
