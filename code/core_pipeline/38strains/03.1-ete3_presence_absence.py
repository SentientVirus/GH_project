#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates a plot of the 38 representative strains, next to which it 
plots the presence/absence of GH70 and GH32 genes.

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

#Representative set of strains
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

scale = 60 #Variable to scale the presence/absence circles

#To color the circles
colors = ['#336E74', '#84CDE8', '#FFE570', '#BDE384', '#84E8A7', '#FF707C', 
          '#FFC36F', '#FF986F', '#C870FF', '#F87BFF', '#656ED4']

# =============================================================================
# 2. Set path to input files
# =============================================================================

treefile1 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.treefile' #Tree with recombination filters
treefile2 = os.path.expanduser('~') + '/GH_project/trees/phylogeny_bootstrap.nofilter.treefile' #Tree without recombination filter
treefile3 = os.path.expanduser('~') + '/GH_project/all_core/38_strains/trees/38strains.treefile' #Tree of the 38 representative strains
treefiles = [treefile1, treefile2, treefile3] #List with all of the tree files

outdir = os.path.expanduser('~') +  '/GH_project/plots/trees' #Directory to save the plots

if not os.path.exists(outdir): #Create the output directory if it does not exist
   os.makedirs(outdir)

# =============================================================================
# 3. Loop through input files and generate the plots
# =============================================================================

for treefile in treefiles: #Loop through tree files
    outfile = f'{outdir}/{treefile.split("/")[-1].split(".treefile")[0]}_GH.png'
    print(outfile)

    t = Tree(treefile, format = 0) #Create a tree object from the tree file
    
    if treefile != treefile3: #If the tree file has more than just the representative strains
        to_keep = [node for node in t.traverse() if node.name in to_include] #Set the nodes to be kept (representative strains)
        t.prune(to_keep) #Prune the tree
    outnode = t.search_nodes(name = 'A1001')[0] #Set the outgroup to the phylogroup F strain A1001
    t.set_outgroup(outnode) #Root the tree
    
    ts = TreeStyle() #Create tree style object
    ts.show_branch_length = False #Hide branch lengths
    ts.show_branch_support = False #Hide branch supports to add formatted text
    ts.show_leaf_name = False #Hide leaf names to add formatted text
    ts.scale =  30000 #Scale of the tree
    ts.scale_length = 0.01 #Scale legend bar
    ts.legend_position = 4 #Legend placement at bottom left
    ts.draw_guiding_lines = True #Draw guiding lines
    ts.guiding_lines_color = 'lightgrey' #Set the color of guiding lines
    
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
                                   ftype = 'Arial', fsize = 30) #Sent font type and size
           n.add_face(support_face, column = 0, position = 'branch-top') #Add the text to the node
       
    for leaf in t.get_leaves(): #Loop through the leaves sorted by phylogeny 
        nleaf = leaf.name.replace('-', '').upper()  #Create leaf names without the '-' symbol
        color = leaf_color.get(nleaf, None) #Use that to retrieve the leaf color
        name_face = TextFace(leaf.name, ftype = 'Arial', fgcolor = color, 
                             fsize = 40) #Create the text with leaf names
        leaf.add_face(name_face, column = 0, position = 'aligned') #Add the text to the leaf
        
        if leaf.name == 'A1001': #If the node is the outgroup
            name_string = '' #Create an empty string to add gene names
            for i in range(len(names)): #Loop through gene names
                gene_text = f'{names[i]}' #Create text with each gene name
                text_face = TextFace(gene_text, fsize = 28, tight_text = True) #Format text
                text_face.rotation = 290 #Rotate text
                text_face.margin_bottom = -9 #Reduce space between labels
                if i == 0: #If it's the first gene name
                    ts.legend.add_face(text_face, column = i) #Add it to the plot
                    empty_face = TextFace(' ', ftype = 'Arial', fsize=28, #Create text with empty spaces
                                          tight_text = True)
                    ts.legend.add_face(empty_face, column = i+1) #Add these spaces to two columns (to separate NGB from the other GHs)
                    ts.legend.add_face(empty_face, column = i+2)
                else: #Otherwise
                    ts.legend.add_face(text_face, column = i+2) #Add the gene name text
                if i == len(names)-1: #If it is the last gene name
                    empty_face = TextFace(' '*2, ftype = 'Arial', fsize=28, 
                                          tight_text = True) #Create text with spaces
                    ts.legend.add_face(empty_face, column = len(names)+2) #Add it after the last gene (to center gene names better)
        
        motifs = ['']*len(names) #Create an empty list of motifs to plot
        seqFace = SeqMotifFace('A'*300, motifs = '', seq_format = 'blank', #Create an empty motif to add space between leaf names and presence/absence 
                               gap_format = 'blank')
        (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned') #Add empty motif to the first column
        
        for n in range(len(names)): #Loop through gene names
            motifs[n] = [5, 55, 'o', None, 50, '', '', ''] #Start adding info to plot a motif
            motifs[n][5] = colors[n] #Update edge colors according to gene type
            motifs[n][6] = colors[n] #Update fill colors according to gene type
            
            # if n == 0: #Uncomment to add gene type on circles
            #     motifs[n][7] = f'Arial|18|white|{names[n]}' 
            # else:
            #     motifs[n][7] = f'Arial|18|black|{names[n]}'
                
            if strain_presence[leaf.name][n] == 0: #If the gene is absent
                motifs[n][0] = 15 #Set start position of circle
                motifs[n][1] = 45 #Set end position of circle (equal to end - start)
                motifs[n][4] = 30 #Set height of circle
                motifs[n][5] = '#E7E7E7' #Set edge color to light grey
                motifs[n][6] = '#E7E7E7' #Set fill color to light grey
                # motifs[n][7] = 'Arial|18|black|'
            elif strain_presence[leaf.name][n] == -1: #If the gene is partially present, perhaps lost
                motifs[n][0] = 10 #Set start position of circle
                motifs[n][1] = 50 #Set end position of circle
                motifs[n][4] = 40 #Set height of circle (equal to end - start)
                motifs[n][5] = 'darkgrey' #Set edge color to dark grey
                motifs[n][6] = 'darkgrey' #Set fill color to dark grey
                # motifs[n][7] = 'Arial|18|black|'
            
            seq = 'A'*scale #Create sequence object on top to which circles are plotted
            if n == 0: #If it's the first gene (NGB)
                seq += 'A'*40 #Increase the length of the sequence object to add space between NGB and the other genes
            seqFace = SeqMotifFace(seq, motifs = [motifs[n]], seq_format = 'blank', #Create object with the motifs, make the sequence tract invisible
                                   gap_format = 'blank') 
            (t & f'{leaf.name}').add_face(seqFace, n+1, 'aligned') #Add each gene to a different column

        
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