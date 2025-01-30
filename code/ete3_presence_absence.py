#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates a plot of the GH70 phylogeny, next to which it plots the
presence/absence of GH70 and GH32 genes.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Load required packages
# =============================================================================
from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace
import os

# =============================================================================
# 1. Define variables to do the formatting
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

#This can be removed after I generate a tree of the 38 strains
to_include = ['A0901', 'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805',
             'DSM', 'fhon2', 'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A',
             'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 'H3B2-02X',
             'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 'H4B1-11J',
             'H4B2-02J', 'H4B2-04J', 'H4B2-05J', 'H4B2-06J', 'H4B2-11M',
             'H4B4-02J', 'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J',
             'H4B5-03X', 'H4B5-04J', 'H4B5-05J', 'IBH001', 'MP2']

#Dictionary of GH genes
domain_no = {0: 'NGB', 1: 'S3', 2: 'GS2', 3: 'GS2-BRS', 4: 'BRS', 5: 'GS1', 
             6: 'GS4', 7: 'GS3', 8: 'S2a', 9: 'S2b', 10: 'S1'}

#List of GH genes
names = list(domain_no.values())

#Matrix with gene counts per strain (0 = absent, 1 = present, -1 = partial/non-functional)
                   #Strain      N  3 G2 G2B B G1 G4 G3 2a 2b  1
strain_presence = {'A0901':    [1, 1, 0, 0, 0,-1, 0, 0, 1, 0, 0],
                   'A1001':    [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0],
                   'A1003':    [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                   'A1202':    [1, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0],
                   'A1401':    [1, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0],
                   'A1404':    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   'A1805':    [1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0],
                   'DSM':      [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1],
                   'fhon2':    [1, 0,-1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'G0101':    [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
                   'G0403':    [0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'H1B1-04J': [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0],
                   'H1B1-05A': [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                   'H1B3-02M': [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'H3B1-01A': [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'H3B1-04J': [0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0],
                   'H3B1-04X': [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'H3B2-02X': [1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],
                   'H3B2-03J': [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'H3B2-03M': [1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0],
                   'H3B2-06M': [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
                   'H3B2-09X': [1, 0,-1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'H4B1-11J': [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
                   'H4B2-02J': [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
                   'H4B2-04J': [1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
                   'H4B2-05J': [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
                   'H4B2-06J': [0, 1, 0, 0, 0,-1, 0, 0, 1, 0, 0],
                   'H4B2-11M': [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
                   'H4B4-02J': [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0],
                   'H4B4-05J': [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
                   'H4B4-06M': [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
                   'H4B4-12M': [0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0],
                   'H4B5-01J': [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
                   'H4B5-03X': [1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0],
                   'H4B5-04J': [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'H4B5-05J': [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                   'IBH001':   [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
                   'MP2':      [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0]
                   }

scale = 60 #To set the size of the dots
#To color the dots
colors = ['#336E74', '#84CDE8', '#FFE570', '#BDE384', '#84E8A7', '#FF707C', 
          '#FFC36F', '#FF986F', '#C870FF', '#F87BFF', '#656ED4']

# =============================================================================
# 2. Set path to input files
# =============================================================================
treefile1 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.treefile' #Tree with recombination filters
treefile2 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.nofilter.treefile' #Tree without recombination filter
treefile3 = os.path.expanduser('~') + '/GH_project/all_core/38_strains/trees/38strains.treefile' #Tree without recombination filter
treefiles = [treefile1, treefile2, treefile3] #List with both of the tree files

outdir = os.path.expanduser('~') +  '/GH_project/plots/trees'

if not os.path.exists(outdir):
   os.makedirs(outdir)

for treefile in treefiles:
    outfile = f'{outdir}/{treefile.split("/")[-1].split(".treefile")[0]}.png'
    print(outfile)

    t = Tree(treefile, format = 0)
    if treefile != treefile3:
        outnode = t.search_nodes(name = 'fhon13')[0] #t.get_common_ancestor('A1001', ['A2002', 'A2003']) #GH70
        t.set_outgroup(outnode)
        to_keep = [node for node in t.traverse() if node.name in to_include] #Get all t>
        t.prune(to_keep) #Prune the tree
    else:
        outnode = t.search_nodes(name = 'A1001')[0]
        t.set_outgroup(outnode)
        outfile = f'{outdir}/representative_38.png'
    
    ts = TreeStyle()
    #ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.scale = 30000
    ts.show_leaf_name = False
    ts.scale_length = 0.01
    ts.legend_position = 4 #Bottom left
    ts.draw_guiding_lines = True
    ts.guiding_lines_color = 'lightgrey'
    
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
               support_face = TextFace(int(n.support), ftype = 'Arial', 
                                       fgcolor = color, fsize = 24)
               n.add_face(support_face, column=0, position='branch-top')
       
    leaves = t.get_leaves() #Sort by phylogeny
    lnames = [leaf.name for leaf in leaves]
    
                    
    # Fix this part for MP2. To fix MP2 successfully, I have to add the old locus tag from the Interproscan in the previous script
    for leaf in leaves:
        nleaf = leaf.name.replace('-', '').upper()
        color = leaf_color.get(nleaf, None)
        name_face = TextFace(leaf.name, ftype = 'Arial', fgcolor = color, 
                             fsize = 40)
        leaf.add_face(name_face, column=0, position='aligned')
        
        if leaf.name == 'A1001':
            name_string = ''
            for i in range(len(names)):
                gene_text = f'{names[i]}'
                text_face = TextFace(gene_text, fsize=28, tight_text = True)
                text_face.rotation = 290
                text_face.margin_bottom = -9 #Reduce space between labels
                if i == 0:
                    ts.legend.add_face(text_face, column = i)
                    empty_face = TextFace(' ', ftype = 'Arial', fsize=28, 
                                          tight_text = True)
                    ts.legend.add_face(empty_face, column = i+1)
                    ts.legend.add_face(empty_face, column = i+2)
                else:
                    ts.legend.add_face(text_face, column = i+2)
                if i == len(names)-1:
                    empty_face = TextFace(' '*2, ftype = 'Arial', fsize=28, 
                                          tight_text = True)
                    ts.legend.add_face(empty_face, column = len(names)+2)
        
        motifs = ['']*len(names)
        seqFace = SeqMotifFace('A'*300, motifs = '', seq_format = 'blank', 
                               gap_format = 'blank') #Add presence/absence info to node
        (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned') #The number represents the column
        for n in range(len(names)):
            motifs[n] = [5, 55, 'o', None, 50, '', '', '']
            motifs[n][5] = colors[n]
            motifs[n][6] = colors[n]
            if n == 0:
                motifs[n][7] = 'Arial|18|white|' #'{domain_no[n]}'
            else:
                motifs[n][7] = 'Arial|18|black|' #'{domain_no[n]}'
                
            if strain_presence[leaf.name][n] == 0:
                motifs[n][0] = 15
                motifs[n][1] = 45
                motifs[n][4] = 30
                motifs[n][5] = '#E7E7E7'
                motifs[n][6] = '#E7E7E7'
                motifs[n][7] = 'Arial|18|black|'
            elif strain_presence[leaf.name][n] == -1:
                motifs[n][0] = 10
                motifs[n][1] = 50
                motifs[n][4] = 40
                motifs[n][5] = 'darkgrey'
                motifs[n][6] = 'darkgrey'
                motifs[n][7] = 'Arial|18|black|'
            
            seq = 'A'*scale
            if n == 0:
                seq += 'A'*40
            seqFace = SeqMotifFace(seq, motifs = [motifs[n]], seq_format = 'blank', 
                                   gap_format = 'blank') #Add presence/absence info to node
            (t & f'{leaf.name}').add_face(seqFace, n+1, 'aligned') #The number represents the column

        
    #print(t)
    t.ladderize(1)
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
    nodeJ.ladderize(1)
    nodeK = t.get_common_ancestor('H4B4-12M', 'H4B5-01J')
    nodeK.ladderize(1)
    nodeL = t.get_common_ancestor('IBH001', 'H4B2-06J')
    nodeL.ladderize(1)
    # t.convert_to_ultrametric()
    t.render(outfile, tree_style = ts)
    t.render(outfile.replace('png', 'tiff'), tree_style = ts)
    t.render(outfile.replace('png', 'svg'), tree_style = ts)