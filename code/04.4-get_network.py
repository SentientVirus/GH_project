#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:44:38 2023

Script to create a phylogenetic network.

@author: Marina Mota Merlo
"""

from matplotlib import pyplot as plt

# A = 'H3B2-03M_GS1'
# B = 'A1001_GS1'
# C = 'MP2_GS1'

# dsAB = 2
# dsBC = 5
# dsAC = 0.85

import networkx as nx
import os
import pandas as pd

GS1 = ['A1001_12310', 'A1003_12540', 'A1202_13520', 'A1401_12750', 
       'A1805_12820', 'FHON2_13540', 'G0101_12800', 'G0403_13100', 
       'H1B104J_13010', 'H1B105A_12300', 'H1B302M_12880', 'H3B101A_13240', 
       'H3B104J_12990', 'H3B104X_13200', 'H3B111M_12560', 'H3B202X_12850', 
       'H3B203J_13370', 'H3B203M_12480', 'H3B206M_12830', 'H3B206M_12840',
       'H3B209X_13340', 'H4B111J_13560', 'H4B111J_13570', 'H4B202J_12880', 
       'H4B204J_13330', 'H4B205J_12990', 'H4B211M_13000', 'H4B402J_12600', 
       'H4B405J_13350', 'H4B405J_13360', 'H4B406M_13450', 'H4B406M_13460',
       'H4B412M_13240', 'H4B501J_12890', 'H4B503X_12670', 'H4B504J_13460',
       'H4B505J_12880', 'APS55_RS03850', 'LDX55_06325', 'K2W83_RS06180']

GS2 = ['G0403_13120', 'H1B302M_12900', 'H3B101A_13260', 'H3B104J_13020_2',
       'H3B104X_13220', 'H3B111M_12590', 'H3B202X_12860', 'H3B203J_13390', 
       'H4B202J_12890_2','H4B204J_13350_2', 'H4B504J_13480', 'H4B505J_12900', 
       'APS55_RS03845', 'LDX55_06335_2', 'K2W83_RS06185']

BRS = ['A1401_12760', 'FHON2_13550', 'G0403_13110', 'H1B302M_12890',
       'H3B101A_13250', 'H3B104J_13020', 'H3B104X_13210', 'H3B203J_13380',
       'H3B209X_13350', 'H4B202J_12890', 'H4B204J_13340', 'H4B204J_13350',
       'H4B504J_13470', 'H4B505J_12890', 'LDX55_06330', 'LDX55_06335',
       'K2W83_RS06185_2']

NGB = ['A0901_05380', 'A1003_04750', 'A1202_05530', 'A1401_04720', 
       'A1805_04920', 'FHON2_04830', 'G0101_04800', 'G0102_04760', 
       'H1B105A_04750', 'H1B302M_04960', 'H3B101A_04720', 'H3B104X_04750', 
       'H3B202X_04770', 'H3B203J_04720', 'H3B203M_04710', 'H3B206M_05060', 
       'H3B209X_04660', 'H4B111J_04920', 'H4B202J_04590', 'H4B204J_04710', 
       'H4B205J_04990', 'H4B211M_04810', 'H4B405J_04950', 'H4B406M_05270', 
       'H4B503X_04680', 'H4B504J_05510', 'H4B505J_04980', 'K2W83_RS02365']

S1 = ['G0101_12790', 'G0102_12700', 'H4B205J_12980', 'H4B211M_12990', 
      'H4B501J_12880', 'K2W83_RS06175']

S2a = ['A0901_13270', 'A1001_12300', 'A1805_12810', 'H1B104J_13000', 
       'H3B203M_12470', 'H4B206J_13400', 'H4B402J_12590', 'H4B412M_13230', 
       'H4B503X_12660']

S2b = ['A1805_12800', 'H1B104J_12990', 'H4B412M_13220']

S3 = ['A0901_13360', 'A1001_12360', 'A1404_13450', 'H3B203M_12520',
      'H4B206J_13490', 'H4B402J_12660', 'H4B412M_13310', 'H4B503X_12760']

# =============================================================================
# 1. Here I define the paths to input and output files
# =============================================================================
home = os.path.expanduser('~')

outdir = f'{home}/GH_project/plots/network'

if not os.path.exists(outdir):
    os.makedirs(outdir)

GH_types = ['GH70', 'GH32']

# OBS! I want to switch to sequence identities instead of dS values, use the percentage identity repset file in tables/

infile = f'{home}/GH_project/tables/percentage_identity_repset.tab' #f'{home}/GH_project/results/GH70/dNdS.tsv'

# outfile = f'{home}/GH_project/plots/network/GH70_identity_network.svg'
ident = 70 #(2/3)*100

# =============================================================================
# 2. Here I open the file with dS values and store each locus tag as a node
# label (perhaps assigning random node colors or coloring by phylogroup), and 
# each pairwise dS value as an edge value.
# =============================================================================

for GH_type in GH_types:
    outfile = f'{home}/GH_project/plots/network/{GH_type}_{int(ident)}identity_network.svg'
    if GH_type == 'GH70':
        GH_list = ['/GS1', '/GS2', '/BRS', '/NGB']
    else:
        GH_list = ['/S1', '/S2', '/S3']

    with open(infile) as reader:
        df = pd.read_csv(reader, sep = '\t')
        node_names = []
        edge_list = []
        tag_dict = {}
        # node_names = list(set(list(df['locus1'])+list(df['locus2'])))
        # edge_list = [(row['locus1'], row['locus2'], row['dS']) for index, row in df.iterrows()]
        for index, row in df.iterrows():
            if any(GH in row['GH_type'] for GH in GH_list):
                node_names += [row[0], row[1]]
                edge_list.append((row[0], row[1], row[3]))
                tag_dict[row[0]] = row[2].split('/')[0]
                tag_dict[row[1]] = row[2].split('/')[1]
    
    node_names = set(node_names)
    color_dict = {'GS1': '#ff0025', 'BRS': '#da00ff', 'GS2': '#ff00a4', 
                  'NGB': '#5a00ff', 'S1': '#ff6d00', 'S2a': '#ffb900', 
                  'S2b': '#ecce23', 'S3': '#d5dd2d'}
    color_edge_dict = {'GS1': '#da707e', 'GS2': '#da70b3', 'BRS': '#cc70da',
                       'NGB': '#9770da', 'S1': '#b59274', 'S2a': '#b5a274', 
                       'S2b': '#b5ac74', 'S3': '#b2b574'}
    G = nx.Graph()
    G.add_nodes_from(node_names)
    G.add_weighted_edges_from(edge_list)
    
    # subset_edges = [edge for edge in edge_list if edge[2] <= 2]
    
    colors = [color_dict[tag_dict[loctag]] for loctag in node_names]
    
    edge_colors = ['#d1d1d1' if tag_dict[edge[0]] != tag_dict[edge[1]] else color_edge_dict[tag_dict[edge[0]]] for edge in list(G.edges())]
    
    # =============================================================================
    # 3. Here I plot everything as a network (using the code above)
    # =============================================================================
    
    fig, ax = plt.subplots()
    ax.margins(0.1)
    ax.axis('off')
    fig.set_size_inches(48, 32)
    #colors = [value['color'] for value in list(G.nodes.values())]
    weights = [pos_w[2]['weight']/10 if pos_w[2]['weight'] > ident else 0 for pos_w in list(G.edges.data())]
    #final_weights = [21/(w*10+1) if w <= 2 else 0 for w in weights]
    
    # This part is needed to specify that I want sequences that are have a high
    # percentage of identity to cluster together, not the opposite
    distances = dict(nx.shortest_path_length(G, weight='weight'))
    for k, v in distances.items():
        for k2, v2 in v.items():
            print(k, k2)
            if distances[k][k2] != 0:
                distances[k][k2] = 1/distances[k][k2]
    
    net_pos = nx.kamada_kawai_layout(G, dist = distances)#, dist = {('A1202_05350', 'FHON2_04830'): 1})
    
    # Divide this into a command to plot edges with lower alpha and a command to plot
    # nodes with higher alpha
    
    nx.draw_networkx_edges(G, pos = net_pos, edgelist = G.edges(), alpha = 0.5,
                            edge_color = edge_colors, width = weights, ax = ax)
    
    nx.draw_networkx_nodes(G, pos = net_pos, nodelist = G.nodes(), alpha = 1, node_color = colors,
                            node_shape = '*', node_size = 1200, ax = ax)
    # nx.draw(G, pos=net_pos, with_labels = False, alpha = 0.8, font_size = 20,
    #         edge_color = edge_colors, width = final_weights, node_color = colors, 
    #         node_shape = '*', node_size = 1200, ax = ax)
    
    nx.draw_networkx_labels(G, pos=net_pos, font_size = 20, font_weight='bold',
                            verticalalignment = 'baseline', ax = ax,
                            horizontalalignment = 'center')
    
    
    # =============================================================================
    # 4. Here I save the network to .svg and .tiff or similar
    # =============================================================================
    
    fig.savefig(outfile, format='svg', dpi=800, pad_inches = 0)
    fig.savefig(outfile.replace('.svg', '.tiff'), format='tiff', dpi=800, pad_inches = 0)
