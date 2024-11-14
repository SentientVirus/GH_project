#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:44:38 2023

Script to create a phylogenetic network.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required libraries and modules
# =============================================================================

from matplotlib import pyplot as plt
import networkx as nx
import os
import pandas as pd

# =============================================================================
# 1. Define gene types
# =============================================================================

GS1 = ['A1001_12310', 'A1202_13520', 'A1401_12750', 'A1805_12820', 
       'FHON2_13540', 'G0101_12800', 'G0403_13100', 'H1B104J_13010', 
       'H1B105A_12300', 'H1B302M_12880', 'H3B101A_13240', 'H3B104J_12990',
       'H3B104X_13200', 'H3B111M_12560', 'H3B202X_12850', 'H3B203J_13370', 
       'H3B203M_12480', 'H3B206M_12840', 'H3B209X_13340', 'H4B111J_13570', 
       'H4B202J_12880', 'H4B204J_13330', 'H4B205J_12990', 'H4B211M_13000', 
       'H4B402J_12600', 'H4B405J_13360', 'H4B406M_13460', 'H4B412M_13240', 
       'H4B501J_12890', 'H4B503X_12670', 'H4B504J_13460', 'H4B505J_12880', 
       'APS55_RS03850', 'LDX55_06325', 'K2W83_RS06180']

GS2 = ['G0403_13120', 'H1B302M_12900', 'H3B101A_13260', 'H3B104J_13020_2',
       'H3B104X_13220', 'H3B111M_12590', 'H3B202X_12860', 'H3B203J_13390', 
       'H4B202J_12890_2','H4B204J_13350_2', 'H4B504J_13480', 'H4B505J_12900', 
       'APS55_RS03845', 'LDX55_06335_2', 'K2W83_RS06185']

GS3 = ['H3B206M_12830', 'H4B111J_13560', 'H4B405J_13350', 'H4B406M_13450']

GS4 = ['A1003_12540']

BRS = ['A1401_12760', 'FHON2_13550', 'G0403_13110', 'H1B302M_12890',
       'H3B101A_13250', 'H3B104J_13020_1', 'H3B104X_13210', 'H3B203J_13380',
       'H3B209X_13350', 'H4B202J_12890_1', 'H4B204J_13340', 'H4B204J_13350_1',
       'H4B504J_13470', 'H4B505J_12890', 'LDX55_06330', 'LDX55_06335_1',
       'K2W83_RS06185_1']

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

GH_types = ['GH70', 'GH32']

# =============================================================================
# 2. Define the paths to input and output files
# =============================================================================
workdir = os.path.expanduser('~') + '/GH_project' #Working directory
outdir = f'{workdir}/plots/network' #Output directory

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

infile = f'{workdir}/tables/percentage_identity_repset.tab' #Input file with % of identity

ident = 70 #Threshold percentage of identity

# =============================================================================
# 3. Open the file with dS values and store each locus tag as a node
# label and each pairwise dS value as an edge value.
# =============================================================================

color_dict = {'GS3': '#FF7594', 'GS4': '#FF8A75', 'GS1': '#FF7575', 
              'BRS': '#E875FF', 'GS2': '#FF75B6', 'NGB': '#C475FF', 
              'S1': '#FFA175', 'S2a': '#FFB875', 'S2b': '#FFD775', 
              'S3': '#FFF575'} #Dictionary to color nodes by gene type

color_edge_dict = {'GS3': '#C96279', 'GS4': '#C96D5C', 'GS1': '#C95E5E',  
                   'GS2': '#C95F91', 'BRS': '#B95ECB', 'NGB': '#995CC6', 
                   'S1': '#C97F5D', 'S2a': '#C9925E', 'S2b': '#C9A95C', 
                   'S3': '#C9C15C'} #Dictionary to color edges by gene type

for GH_type in GH_types: #Loop through GH families (GH70/32)
    outfile = f'{workdir}/plots/network/{GH_type}_{int(ident)}identity_network.svg' #Set the name of the output file
    
    if GH_type == 'GH70': #Set the subtypes of GH70 genes
        GH_list = ['/GS1', '/GS2', '/GS3', '/GS4', '/BRS', '/NGB']
    else: #Set the subtypes of GH32 genes
        GH_list = ['/S1', '/S2', '/S3']

    with open(infile) as reader: #Open input file
        df = pd.read_csv(reader, sep = '\t') #Load file contents as a dataframe
        node_names = [] #Empty list to store node information (pairs of nodes)
        edge_list = [] #Empty list to store edge information (pairs of nodes and identity)
        tag_dict = {} #Dictionary to assign a gene type to each node

        for index, row in df.iterrows(): #Loop trough rows in the dataframe
            if any(GH in row['GH_type'] for GH in GH_list): #If the GH belongs to the right family
                node_names += [row[0], row[1]] #Add node information
                edge_list.append((row[0], row[1], row[3])) #Add edge information
                tag_dict[row[0]] = row[2].split('/')[0] #Add gene type to first locus tag
                tag_dict[row[1]] = row[2].split('/')[1] #Add gene type to second locus tag
    
    node_names = set(node_names) #Convert list to set

    G = nx.Graph() #Create the network graph
    G.add_nodes_from(node_names) #Add nodes to graph
    G.add_weighted_edges_from(edge_list) #Add edges to graph, considering weights (% identity)
    
    colors = [color_dict[tag_dict[loctag]] for loctag in node_names] #Get the colors of each node
    
    edge_colors = ['#d1d1d1' if tag_dict[edge[0]] != tag_dict[edge[1]] else color_edge_dict[tag_dict[edge[0]]] for edge in list(G.edges())] #Get edge colors (grey if the genes belong to different subtypes)
    
    # =============================================================================
    # 3.Plot everything as a network (using the code above)
    # =============================================================================
    
    fig, ax = plt.subplots() #Create plot
    ax.margins(0.1) #Set plot margins
    ax.axis('off') #Remove the plot axis
    fig.set_size_inches(48, 32) #Set figure size

    weights = [pos_w[2]['weight']/10 if pos_w[2]['weight'] > ident else 0 for pos_w in list(G.edges.data())] #Adjust the weights of the edges and remove edges with identity < 70%
    distances = dict(nx.shortest_path_length(G, weight='weight')) #Convert the weights to distances
    for k, v in distances.items(): #Loop through node names (k = node 1, v = node 2 and distance)
        for k2, v2 in v.items(): #Loop through node names and distances (k2 = node 2, v2 = distance)
            print(k, k2) #Print the name of both nodes
            if distances[k][k2] != 0: #If the distance is not 0
                distances[k][k2] = 1/distances[k][k2] #Reverse the distance (so that the higher the % identity, the lower the distance)
    
    net_pos = nx.kamada_kawai_layout(G, dist = distances) #Generate network
    
    # Divide this into a command to plot edges with lower alpha and a command to plot
    # nodes with higher alpha
    
    nx.draw_networkx_edges(G, pos = net_pos, edgelist = G.edges(), alpha = 0.5,
                            edge_color = edge_colors, width = weights, ax = ax) #Plot edges
    
    nx.draw_networkx_nodes(G, pos = net_pos, nodelist = G.nodes(), alpha = 1, 
                           node_color = colors, node_shape = '*', 
                           node_size = 1200, ax = ax) #Plot nodes
    
    nx.draw_networkx_labels(G, pos=net_pos, font_size = 20, font_weight='bold',
                            verticalalignment = 'baseline', ax = ax,
                            horizontalalignment = 'center') #Add labels to nodes
    
    
    # =============================================================================
    # 4. Here I save the network to files
    # =============================================================================
    
    fig.savefig(outfile, format='svg', dpi=800, pad_inches = 0) #Save network to SVG
    fig.savefig(outfile.replace('.svg', '.png'), format='png', dpi=800, 
                pad_inches = 0) #Save network to PNG
    fig.savefig(outfile.replace('.svg', '.tiff'), format='tiff', dpi=800, 
                pad_inches = 0) #Save netwoek to TIFF
