#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

This script creates plots of the domains of GH70 and GH32 proteins
next to a strain phylogeny. To plot a specific gene type, the script has
to be manually updated.

@author: Marina Mota Merlo
"""

from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace, RectFace, CircleFace
import re, os

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


shapes = {'RE': '[]', 'EL': '()', 'DI': '<>', 'TR': 'o'}

GH_types = ['GH32', 'GH70']

workdir = os.path.expanduser('~') + '/GH_project'

for GH_type in GH_types:

    domain_file = f'{workdir}/data/tabs/{GH_type}_domain_file.txt'
    treefile = f'{workdir}/data/fasta/{GH_type}/trees/{GH_type}_functional_outgroup_repset.mafft.faa.treefile' #'data/fasta/GH32/trees/GH32_repset.mafft.faa.treefile' #'data/fasta/GH70/trees/GH70_functional_outgroup_repset.mafft.faa.treefile'
    outfile = f'{workdir}/plots/trees/{GH_type}_domains.png'
    
    if GH_type == 'GH70':
        t = Tree(treefile, format = 0)
        outnode = t.get_common_ancestor('AAU08014.2_L_reuteri', 'AOR73699.1_L_fermentum') #GH70
    elif GH_type == 'GH32':
        t = Tree(treefile.replace('functional_outgroup_', ''), format = 0)
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
           if n.support >= 80:
               support_face = TextFace(int(n.support), fgcolor = color, fsize = 24)
               n.add_face(support_face, column=0, position='branch-top')
       
    leaves = t.get_leaves() #Sort by phylogeny
    lnames = [leaf.name for leaf in leaves]
    
    with open(domain_file, 'r') as dfile:
        domains = (d for d in dfile)
        dom_dict = {}
        seq_dict = {}
        for domain in domains:
            domain = domain.replace('\n', '')
            domain = re.split(r'[|,]', domain)
            if len(domain) >= 3:
                seq_dict[domain[0]] = 'A'*int(domain[1])
                dom_dict[domain[0]] = []
                n = (len(domain)-1)
                i = 2
                while i < n: #Change signal peptide to blue circle
                    if domain[i+4] == 'SP' or domain[i+4] == 'CWB':
                        size = 15
                        if domain[i+4] == 'CWB':
                            domain[i+3] = '#cce769'
                            domain[i+4] = 'CB'
                        else:
                            domain[i] = 'TR'
                            domain[i+3] = '#9975ff'
                    else:
                        size = 24
                        if domain[i+4] == 'GH70':
                            domain[i+3] = '#FF7575' #'#FFB875' # #C70039 (all) #B600FF (NCB), #D930BD (BRS), #FF0000 (GS1), #FF009B (GS2), #7000FF (short)
                        elif domain[i+4] == 'GH32':
                            domain[i+3] = '#FFB875'
                        elif domain[i+4] == 'CB':
                            domain[i+4] = 'GB'
                            domain[i+3] = '#75cd5e'
                        elif domain[i+4] == 'DUF':
                            domain[i+3] = '#e875ff'
                    if domain[i+4] != 'GH70' and domain[i+4] != 'GH32':
                        dom_dict[domain[0]].append([int(domain[i+1]), int(domain[i+2]), 
                        shapes[domain[i]], None, 40, domain[i+3], domain[i+3], None])
                    else:
                        dom_dict[domain[0]].append([int(domain[i+1]), int(domain[i+2]), 
                        shapes[domain[i]], None, 40, domain[i+3], domain[i+3],
                        f'arial|{size}|black|{domain[i+4]}'])
                        
                    i += 5
                    
    # Fix this part for MP2. To fix MP2 successfully, I have to add the old locus tag from the Interproscan in the previous script
    for leaf in leaves:
        # if leaf.name[-2:] == '_2':
        #     nleaf = leaf.name[:-2]
        # else:
        nleaf = leaf.name
        if leaf.name[-2] == '_':
            gene_name = leaf.name[:-2]
        else: gene_name = ''
        if nleaf in seq_dict.keys():
            seqFace = SeqMotifFace(seq_dict[nleaf], motifs = dom_dict[nleaf], seq_format = 'line', scale_factor = 0.8)
            (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned')
        elif gene_name in seq_dict.keys():
            seqFace = SeqMotifFace(seq_dict[gene_name], motifs = dom_dict[gene_name], seq_format = 'line', scale_factor = 0.8)
            (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned')
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
        
        
    #print(t)
    ts.layout_fn = lambda node: True
    ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#75cd5e', label = ''), column=0)
    ts.legend.add_face(TextFace(' Glucan-binding domain', fsize = 32), column = 1)
    if GH_type == 'GH70':
        ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#cce769', label = ''), column=0)
        ts.legend.add_face(TextFace(' Cell wall-binding domain', fsize = 32), column = 1)
        ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#FF7575', label = ''), column=0)
        ts.legend.add_face(TextFace(' GH70 domain', fsize = 32), column = 1)
        ts.legend.add_face(CircleFace(30, color = '#e875ff', label = ''), column=0)
        ts.legend.add_face(TextFace(' DUF5776', fsize = 32), column = 1)
    elif GH_type == 'GH32':
        ts.legend.add_face(RectFace(60, 40, fgcolor = None, bgcolor = '#FFB875', label = ''), column=0)
        ts.legend.add_face(TextFace(' GH32 domain', fsize = 32), column = 1)
    ts.legend.add_face(CircleFace(30, color = '#9975ff', label = ''), column=0)
    ts.legend.add_face(TextFace(' Signal peptide', fsize = 32), column = 1)
    ts.legend_position = 2
    t.ladderize(1)
    # t.convert_to_ultrametric()
    t.render(outfile, tree_style = ts)
