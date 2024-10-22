#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates plots of the domains of GH70 and GH32 proteins
next to a strain phylogeny. To plot a specific gene type, the script has
to be manually updated.

@author: Marina Mota Merlo
"""

from ete3 import Tree, TreeStyle, NodeStyle
from ete3 import SeqMotifFace, TextFace, RectFace, CircleFace
import re, os

# =============================================================================
# 1. Define inputs and plot formatting variables
# =============================================================================
#Dictionary to color strains by phylogroup
leaf_color = {'A0901': '#D55E00', 'A1001': '#771853', 'A1002': '#D55E00',
              'A1003': '#0072B2', 'A1201': '#D55E00', 'A1202': '#33B18F',
              'A1401': '#33B18F', 'A1404': '#FF74D6', 'A1802': '#D55E00',
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
              'G0803': '#0072B2', 'G0804': '#0072B2', 'Fhon2': '#0072B2',
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
              'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': 'black', 
              'DSMZ12361': 'black'}

#Read text as a particular shape
shapes = {'RE': '[]', 'EL': '()', 'DI': '<>', 'TR': 'o'}

#Glycosyl hydrolase types to make graphical representations from
GH_types = ['GH32', 'GH70']

workdir = os.path.expanduser('~') + '/GH_project' #Work directory

# =============================================================================
# 2. Create the plot
# =============================================================================
for GH_type in GH_types: #Loop through GH types

    domain_file = f'{workdir}/data/tabs/{GH_type}_domain_file.txt' #File with domain information
    treefile = f'{workdir}/data/fasta/{GH_type}/trees/{GH_type}_functional_outgroup_repset.mafft.faa.treefile' #Tree file
    outfile = f'{workdir}/plots/trees/{GH_type}_domains.png' #Output file
    
    if GH_type == 'GH70': #If the gene type is GH70
        t = Tree(treefile, format = 0) #Read the tree file
        outnode = t.get_common_ancestor('AAU08014.2_L_reuteri', 'AOR73699.1_L_fermentum') #Get the root of the tree
    elif GH_type == 'GH32': #If the gene type is GH32
        t = Tree(treefile.replace('functional_outgroup_', ''), format = 0) #Read the tree file
        outnode = t.get_common_ancestor('H3B203M_12520', 'A1404_13450') #Get the root of the tree

    t.set_outgroup(outnode) #Set the root of the tree
    
    ts = TreeStyle() #Create a tree style
    ts.scale = 500 #Set scale of the tree
    ts.show_branch_length = False #Remove branch length numbers
    ts.show_branch_support = False #Remove support values
    ts.show_leaf_name = False #Remove leaf names
    
    ns = NodeStyle() #Create node style
    ns['size'] = 0 #Remove circles over the nodes
    ns['vt_line_width'] = 5 #Set width of vertical lines
    ns['hz_line_width'] = 5 #Set width of horizontal lines
    ns['hz_line_type'] = 0 #Use solid horizontal lines
    
    for n in t.traverse(): #Loop through nodes in the tree
       n.set_style(ns) #Assign style to the nodes
       if n not in t.get_leaves() and n.support > 80: #If the node is not a leaf and the support is over 80
           if n.support >= 95: #If the support is over 95
               color = 'black' #Color support values in black
           else: color = 'dimgrey' #Otherwise, color support values in grey
           support_face = TextFace(int(n.support), fgcolor = color, fsize = 24)
           n.add_face(support_face, column=0, position='branch-top')
       
    leaves = t.get_leaves() #Get tree leaves
    lnames = [leaf.name for leaf in leaves] #Create a list of leaf names following the phylogeny 
    
    with open(domain_file, 'r') as dfile: #Open file with domain information
        domains = (d for d in dfile) #Get all the information in the file
        dom_dict = {} #Create empty dictionary with domain information
        seq_dict = {} #Create empty dictionary with sequences
        for domain in domains: #Loop through domains
            domain = domain.replace('\n', '') #Remove special characters from line
            domain = re.split(r'[|,]', domain) #Separate elements in line
            if len(domain) >= 3: #If the line contains at least 3 elements
                seq_dict[domain[0]] = 'A'*int(domain[1]) #Use the first element (locus tag) as a key and create a sequence of a length indicated by the second element
                dom_dict[domain[0]] = [] #Create an empty list to later store domain information
                n = (len(domain)-1) #Get length of line
                i = 2 #Loop starting at the third element
                while i < n: #Loop until reaching the end of the list
                    #Assign colors and shapes to each domain depending on the type
                    if domain[i+4] == 'CWB':
                        domain[i+3] = '#cce769'
                        domain[i+4] = 'CB'
                    elif domain[i+4] == 'SP':
                        domain[i] = 'TR'
                        domain[i+3] = '#9975ff'
                    elif domain[i+4] == 'GH70':
                        domain[i+3] = '#FF7575'
                    elif domain[i+4] == 'GH32':
                        domain[i+3] = '#FFB875'
                    elif domain[i+4] == 'CB':
                        domain[i+4] = 'GB'
                        domain[i+3] = '#75cd5e'
                    elif domain[i+4] == 'DUF':
                        domain[i+3] = '#e875ff'
                        
                    #Add the shape to a dictionary, don't add text on the shape unless the domain is a GH
                    if domain[i+4] != 'GH70' and domain[i+4] != 'GH32':
                        dom_dict[domain[0]].append([int(domain[i+1]), int(domain[i+2]), 
                        shapes[domain[i]], None, 40, domain[i+3], domain[i+3], None])
                    else:
                        dom_dict[domain[0]].append([int(domain[i+1]), int(domain[i+2]), 
                        shapes[domain[i]], None, 40, domain[i+3], domain[i+3],
                        f'arial|24|black|{domain[i+4]}'])
                        
                    i += 5 #Jump five positions in the list (next domain)
                    
    for leaf in leaves: #Loop through leaves in the tree
        nleaf = leaf.name #Get name of the leaf
        if leaf.name[-2] == '_': #If the leaf corresponds to a double-domain gene
            gene_name = leaf.name[:-2] #Use the locus tag of the complete gene
        else: gene_name = '' #Otherwise, don't get an extra locus tag
        if nleaf in seq_dict.keys(): #If the leaf name is in the dictionary with sequences
        #Create a graphical display of all the domains in the gene
            seqFace = SeqMotifFace(seq_dict[nleaf], motifs = dom_dict[nleaf], 
                                   seq_format = 'line', scale_factor = 0.8)
            (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned') #Add display to the leaf
        elif gene_name in seq_dict.keys(): #If the full gene is in the sequence dictionary
            seqFace = SeqMotifFace(seq_dict[gene_name], motifs = dom_dict[gene_name], seq_format = 'line', scale_factor = 0.8)
            (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned')
        #Change the locus tags of certain strains to make them more recognisable
        if 'LDX55' in nleaf:
            nleaf = nleaf.replace('LDX55', 'IBH001')
        elif 'APS55' in nleaf:
            nleaf = nleaf.replace('APS55_RS', 'MP2_')
        elif 'K2W83' in leaf.name:
            nleaf = nleaf.replace('K2W83_RS', 'DSM_')
        
        #Use a different leaf name to retrieve Fhon2 from the color dictionary
        mleaf = nleaf
        if 'FHON2' in leaf.name:
            mleaf = mleaf.replace('FHON', 'Fhon')

        color = leaf_color.get(mleaf.split('_')[0], None) #Get leaf text color
        name_face = TextFace(nleaf, fgcolor = color, fsize = 40) #Create leaf name text
        leaf.add_face(name_face, column=0, position='branch-right') #Add leaf name text to the plot
        
    #Set length of the legend padding to the right
    if GH_type == 'GH70':
        no_blank = 20
    elif GH_type == 'GH32':
        no_blank = 10
    
    #Add legend patches (rectangles or circles) in the same color as domain representations
    ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#75cd5e', label = ''), column=0)
    #Add text with domain name next to the legend patch
    ts.legend.add_face(TextFace(' Glucan-binding domain', fsize = 32), column = 1)
    #Add padding at the right margin
    ts.legend.add_face(TextFace(' '*no_blank, fsize = 32), column = 2)
    #Add different domaint types to the legend depending on gene type
    if GH_type == 'GH70':
        ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#cce769', label = ''), column=0)
        ts.legend.add_face(TextFace(' Cell wall-binding domain', fsize = 32), column = 1)
        ts.legend.add_face(TextFace(' '*no_blank, fsize = 32), column = 2)
        ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#FF7575', label = ''), column=0)
        ts.legend.add_face(TextFace(' GH70 domain', fsize = 32), column = 1)
        ts.legend.add_face(TextFace(' '*no_blank, fsize = 32), column = 2)
        ts.legend.add_face(CircleFace(30, color = '#e875ff', label = ''), column=0)
        ts.legend.add_face(TextFace(' DUF5776', fsize = 32), column = 1)
        ts.legend.add_face(TextFace(' '*no_blank, fsize = 32), column = 2)
    elif GH_type == 'GH32':
        ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#FFB875', label = ''), column=0)
        ts.legend.add_face(TextFace(' GH32 domain', fsize = 32), column = 1)
        ts.legend.add_face(TextFace(' '*no_blank, fsize = 32), column = 2)
    ts.legend.add_face(CircleFace(30, color = '#9975ff', label = ''), column=0)
    ts.legend.add_face(TextFace(' Signal peptide', fsize = 32), column = 1)
    ts.legend.add_face(TextFace(' '*no_blank, fsize = 32), column = 2)
    
    ts.legend_position = 2 #Position the legend at the top right corner
    t.ladderize(1) #Reverse the order of the leaves in the tree so that the outgroup appears at the bottom
    t.render(outfile, tree_style = ts) #Save the tree to a PNG file
    t.render(outfile.replace('png', 'tiff'), tree_style = ts) #Save the tree to a TIFF file
    t.render(outfile.replace('png', 'svg'), tree_style = ts) #Save the tree to an SVG file
