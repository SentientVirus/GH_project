#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:54 2021

Code to plot a representation of the predicted recombination tracts in RDP5
next to the strain phylogeny, using a subset of 36 representative strains,
from which basal strains have been excluded.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Load the packages that are needed
# =============================================================================

from ete3 import Tree, TreeStyle, NodeStyle, SeqMotifFace, TextFace, RectFace
import os
import pandas as pd

# =============================================================================
# 1. Create a cladd to store recombination segment information
# =============================================================================

#OBS! Update to show also major parent 
class recomb_obj: #Object with recombination tract information
    def __init__(self, recombinant, minor, major, start, end, RDPvalue): #Input values
        self.recombinant = recombinant #Recombinant strain
        self.minor = minor #Minor parent
        self.major = major #Major parent
        self.start = start #Start of the recombination tract
        self.end = end #End of the recombination tract
        self.pvalue = RDPvalue #p-value of the recombination tract
    def __str__(self): #Function to print the object
        return f'<{self.recombinant}: {self.major} acquired a segment ({self.start}:{self.end}), p-value {self.pvalue} from {self.minor}>'

# =============================================================================
# 2. Define input and output paths and formatting variables
# =============================================================================

#Dictionary of leaf colors
leaf_color = {'A0901': '#D55E00', 'A1001': 'white', 'A1002': '#D55E00',
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

#Strains to include in the tree
to_include = ['A0901', 'A1001', 'A1003', 'A1202', 'A1401', 'A1805',
              'DSMZ12361', 'fhon2', 'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A', 
              'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 'H3B2-02X',
              'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 'H4B1-11J', 
              'H4B2-02J', 'H4B2-04J', 'H4B2-05J', 'H4B2-06J', 'H4B2-11M', 
              'H4B4-02J', 'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 
              'H4B5-03X', 'H4B5-04J', 'H4B5-05J', 'IBH001', 'MP2']

colors = ['#9DAAB7', '#6586A3', '#647D88', '#626D75', '#B1BFD4'] #Colors for the CDS
colors = colors[::-1] #Reverse the color list to make the colors consistent with the breakpoint plot

aln_len = 31420 #Total length of the alignment
middle_point = 15970 #Point where the two consecutive segments are joined
treefile = os.path.expanduser('~') + '/GH_project/trees/kunkeei_nonclosed.tree' #Path to the treefile
outplot = os.path.expanduser('~') + '/GH_project/plots/trees/RDP5_tree.png' #Path to the output plot
recomb = os.path.expanduser('~') + '/GH_project/RDP5_analysis/RDP_output/all_subsets_7methods_5+_filtered.csv' #Path to the recombination results
in_tab = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab/all_subsets_positions_aln.tab' #Path to the tabs with gene position information

# =============================================================================
# 2. Format the tree
# =============================================================================
#OBS! Try to prune the nodes, rather than deleting them!
t = Tree(treefile, format = 1) #Load treefile

t.set_outgroup('A1001') #Set the root of the tree

to_keep = [node for node in t.traverse() if node.name in to_include] #Get all the nodes to include in the final tree
t.prune(to_keep) #Prune the tree

for node in t.traverse(): #Loop through nodes in the tree
    if node.name == 'A1001': #If the node is A1001
        root_node = node #Save it to a variable
    elif node.name == 'fhon2': #If the node is fhon2
        node.name = 'Fhon2' #Change the node name to Fhon2
    elif node.name == 'H4B2-06J': #If the node name is H4B2-06J
        IBH_sis = node #Create a sister node to print IBH001
    elif node.name == 'H3B1-04J': #If the node is H3B1-04J
        DSM_sis = node #Create a sister node to print DSM

        
# This works now, but I should also add sister nodes where DSM and IBH001 should be placed.
root_node.add_sister(name = 'IBH001', dist = 0) #Add sister node to include IBH001 at the root of the phylogeny
root_node.add_sister(name = 'DSMZ12361', dist = 0) #Add sister node to include DSM at the root of the phylogeny
IBH_sis.add_child(name = '**', dist = 0) #Add sister node to H4B2-06J to point out to IBH001
IBH_sis.add_child(name = IBH_sis.name, dist = 0) #Re-add H4B2-06J as a leaf (it became an internal node when adding IBH001)
IBH_sis.name = '' #Remove the name of the former H4B2-06J node
DSM_sis.add_child(name = '*', dist = 0) #Add sister node to H3B1-04J to point out to IBH001
DSM_sis.add_child(name = DSM_sis.name, dist = 0) #Re-add strain H3B1-04J as a leaf
DSM_sis.name = '' #Remove name from the former H3B1-04J node

ts = TreeStyle() #Create a tree style

ts.show_branch_length = False #Hide branch lengths
ts.scale = 12000 #General tree scale, original 4000
ts.branch_vertical_margin = 10 #Space between branches
ts.show_leaf_name = False #Hide unformatted leaf names
ts.show_scale = False #Hide tree scale
segment_scale = 0.04 #Scale of the segments compared to the tree
segment_height = 20 #Height of the segments and tracts

# =============================================================================
# 3. Add legend to the tree
# =============================================================================

ts.legend_position = 2 #Place the legend at the upper right corner of the plot
# ts.legend.margin_right = 0 #-100

#Add color legend for the p-values to the tree
ts.legend.add_face(TextFace('p-value', fsize = 18, fstyle = 'italic', 
                            bold = True), column = 0) #Legend title
ts.legend.add_face(TextFace('', fsize = 13), column = 1) #Empty legend next to title
ts.legend.add_face(TextFace('<1E-150', fsize = 16), column = 0) #Legend label
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#944654', 
                            label = None), column = 1) #Legend rectangle
ts.legend.add_face(TextFace('<1E-100', fsize = 16), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#C46877', 
                            label = None), column = 1)
ts.legend.add_face(TextFace('<1E-50', fsize = 16), column = 0)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#DB8484', 
                            label = None), column = 1)
#Plot half of the p-values on separate columns
ts.legend.add_face(TextFace('', fsize = 18), column = 2) #Empty label column
ts.legend.add_face(TextFace('', fsize = 18), column = 3) #Empty legend rectangle column
ts.legend.add_face(TextFace('  <1E-25  ', fsize = 16), column = 2) #Legend label
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#F7B9B0', 
                            label = None), column = 3) #Legend rectangle
ts.legend.add_face(TextFace('  <1E-10  ', fsize = 16), column = 2)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#F7CFC1', 
                            label = None), column = 3)
# ts.legend.add_face(TextFace('  <0.05  ', fsize = 12), column = 2)
# ts.legend.add_face(RectFace(45, 35, fgcolor = 'lightgrey', bgcolor = 'linen', 
#                             label = None), column = 3)

ts.legend.add_face(TextFace('\t'*5, fsize = 18), column = 4) #Add gap between the two legends

#Add color legend for the phylogroups
ts.legend.add_face(TextFace('Phylogroup', fsize = 18, bold = True), 
                    column = 5) #Title of the phylogroup legend
ts.legend.add_face(TextFace(' ', fsize = 16), column = 6) #Empty legend next to title
ts.legend.add_face(TextFace('  Phylogroup A', fsize = 16), column = 5) #Legend label
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#0072B2',
                            label = None), column = 6) #Legend rectangle
ts.legend.add_face(TextFace('  Phylogroup B', fsize = 16), column = 5)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#33B18F', 
                            label = None), column = 6)
ts.legend.add_face(TextFace('  Phylogroup C', fsize = 16), column = 5)
ts.legend.add_face(RectFace(45, 15, fgcolor = 'lightgrey', bgcolor = '#D55E00', 
                            label = None), column = 6)

ts.legend.add_face(TextFace('\t'*5, fsize = 15), column = 7) #Add gap between the two legends

# #Add text legend for major and minor parents
# ts.legend.add_face(TextFace('Parents', fsize = 15, bold = True), 
#                     column = 8) #Title of the parent legend
# ts.legend.add_face(TextFace(' ', fsize = 13), column = 9) #Empty legend next to title
# ts.legend.add_face(TextFace('     Minor', fsize = 13), column = 8) #Legend label
# ts.legend.add_face(TextFace('m', fsize = 12, ftype = 'UbuntuCondensed'), column = 9) #Legend text
# ts.legend.add_face(TextFace('     Major', fsize = 13), column = 8) #Legend label
# ts.legend.add_face(TextFace('M', fsize = 12, ftype = 'UbuntuCondensed'), column = 9) #Legend text

# =============================================================================
# 4. Format the nodes (and branch lines)
# =============================================================================

ns = NodeStyle() #Create a node style
ns['size'] = 0 #Hide the nodes
ns['hz_line_color'] = 'lightgrey' #Set the color of horizontal lines
ns['hz_line_width'] = 4 #Set the width of horizontal lines
ns['vt_line_color'] = 'slategrey' #Set the color of vertical lines
ns['vt_line_width'] = 4 #Set the width of vertical lines

ns2 = NodeStyle() #Create another node style to hide branches
ns2['size'] = 0 #Hide the node
ns2['hz_line_color'] = 'white' #Set the width of horizontal lines
ns2['vt_line_width'] = 0 #Set the width of vertical lines

ns3 = NodeStyle() #Create another node style for the branches that were added manually
ns3['size'] = 0 #Hide the nodes
ns3['hz_line_color'] = 'lightgrey' #Set the color of horizontal lines
ns3['hz_line_width'] = 2 #Set the width of horizontal lines
ns3['vt_line_color'] = 'slategrey' #Set the color of vertical lines
ns3['vt_line_width'] = 2 #Set the width of vertical lines
ns3['hz_line_type'] = 2 #Make horizontal lines dashed

#List of nodes that should be invisible
invisible_nodes = [t.get_common_ancestor('A1001'), 
                   t.get_common_ancestor('IBH001', 'DSMZ12361'), 
                   t.get_common_ancestor('IBH001', 'A1003'), 
                   t.search_nodes(name = 'IBH001')[0], 
                   t.search_nodes(name = 'DSMZ12361')[0], 
                   t.search_nodes(name = 'A1001')[0]]

for n in t.traverse(): #Loop through nodes
    if n not in invisible_nodes and '*' not in n.name: #If the nodes are not invisible/manually added
        n.set_style(ns) #Apply the first style
    elif '*' in n.name: #If the nodes are manually added
        n.set_style(ns3) #Apply the dashed line style
    else: #Otherwise
        n.set_style(ns2) #Apply the style to remove the lines

# =============================================================================
# 5. Retrieve recombination tracts
# =============================================================================

#Get strain names
leaves = t.get_leaves() #Get leaf nodes sorted by phylogeny
lnames = [leaf.name for leaf in leaves] #Store leaf names in a list
lnames_sort = sorted(lnames) #Sort leaf names alphabetically
lno = {lname: lnames_sort.index(lname)+1 for lname in lnames_sort} #Assign a number to each leaf
lnames.reverse() #Reverse the list

recomb_tracts = {leaf.name: [] for leaf in leaves} #Create empty lists to store recombination tracts for each leaf
with open(recomb) as rectracts: #Open RDP5 output
    tracts = (line for line in rectracts) #Create a generator to loop through lines in the file
    for i in range(16): #Skip the first 16 lines
        next(tracts)
    for l in tracts: #Loop through the remaining lines
        if '=' not in l and len(l.split(',')) > 1: #If the line does not contain = and its length is > 1
            l = l.split(',') #Split the line by commas
            if (l[2] and l[3]) != '': #If the third and fourth element of the line are not empty
                recomb = l[8].strip('^') #Get the recombinant strain
                start = ''.join(c for c in l[2] if c.isdigit()) #Get the start position
                end = ''.join(c for c in l[3] if c.isdigit()) #Get the end position
                support = [l[11], l[12], l[13], l[14], l[15], l[16], l[17]] #Get the p-values of all methods
                while 'NS' in support: #If the event is not supported by any method
                    support.remove('NS') #Remove the method from the list
                support = [float(sup) for sup in support] #Convert the strings to floating point numbers
                if support != []: #If the list of support values is not empty
                    support = sorted(support, reverse = False)[:3] #Get the three lowest p-values
                    support = sum(support)/len(support) #Calculate the support as the average
                else: support = 1 #Otherwise, set the support to one
                new_tract = [l[8], l[9], l[10], start, end, support] #Store all the tract information
                recomb_tracts[recomb].append(new_tract) #Append the information to the dictionary
           
# =============================================================================
# 6. Retrieve gene information to plot at the bottom of the figure
# =============================================================================

loci = [] #Create an empty list to append gene information
labels = {} #Create an empty dictionary with gene labels
countvar = 0 #Create a count variable
with open(in_tab) as tabfile: #Open the file with positional information
    tab_df = pd.read_csv(tabfile, sep = '\t') #Load the file content into a dataframe
    tab_df = tab_df.reindex(index = tab_df.index[::-1]) #Reverse the index of the dataframe
    for index, row in tab_df.iterrows(): #Loop through rows in the dataframe
        if row['strain'] == 'A0901': #If the strain is A0901 (it could have been almost any other)
            if row['strand'] == -1: #If the strand is the reverse strand
                shape = '>' #Make the gene point to the right
            else: #Otherwise
                shape = '<' #Make the gene point to the left
            labels[countvar] = row['gene_name'] #Save the gene name as the label for the gene
            countvar += 1 #Increase the count
            locus = [row['start']-1, row['end']-1, shape, 0, segment_height*2, 
                     'black', colors[countvar%5],
                     f'Timmana|16|black|{row["gene_name"]}'] #Format the graphic representation of the gene
            if row['gene_name'] == 'bcrA': #If the gene is bcrA
                locus[0] += 1 #Increase the start position by 1 (so that it doesn't overlap with the previous gene)
            loci.append(locus) #Append the locus to the list

#Annotate the CDS of hypothetical proteins manually
loci[7][7] = loci[7][7].replace('unk.', 'CDS1')
loci[9][7] = loci[9][7].replace('unk.', 'CDS2')
loci[10][7] = loci[10][7].replace('unk.', 'CDS3')
loci[11][7] = loci[11][7].replace('unk.', 'CDS4')
loci[13][7] = loci[13][7].replace('unk.', 'CDS5')
loci[14][7] = loci[14][7].replace('unk.', 'CDS7')
loci[16][7] = loci[16][7].replace('unk.', 'CDS8')

# =============================================================================
# 7. Add the recombination tract and gene information to the tree
# =============================================================================

seq = 'A'*middle_point #Create the first half of the segment to be plotted
seq2 = 'A'*(aln_len-middle_point) #Create the second half of the segment
miniseq = 'A'*20 #Create a smaller sequence
for leaf in leaves: #Loop through the leaves in the tree
    color = leaf_color.get(leaf.name.replace('-', '').upper(), None) #Color leaves according to strain phylogroup
    if 'DSM' in leaf.name: #If the leaf name is DSM
        leaf_name = '*DSM' #Add an asterisk before
    elif 'IBH001' in leaf.name: #If the leaf name is IBH001
        leaf_name = '**' + leaf.name #Add two asterisks before
    else: leaf_name = leaf.name #Otherwise, store the leaf name as it is
    
    name_face = TextFace(leaf_name, fgcolor = color, fsize = 16) #Create a TextFace with the gene name
    leaf.add_face(name_face, column = 0, position='branch-right') #Add formatted leaf names to the tree
    
    seqFace = SeqMotifFace(miniseq, height = segment_height, seq_format = 'blank', gap_format = 'blank') #Add separator between segments
    (t & f'{leaf.name}').add_face(seqFace, 0, 'aligned') #Add the SeqMotifFace to column 0 of the plot
    (t & f'{leaf.name}').add_face(seqFace, 2, 'aligned') #Add the SeqMotifFace to column 2 of the plot

    
    if leaf.name in recomb_tracts.keys(): #If there are recombination tracts for that leaf
        recomb_tract = recomb_tracts[leaf.name] #Retrieve the recombination tracts
        motifs = [] #Create an empty list to store the motifs
    
        for tract in recomb_tract: #Loop through recombination tracts for the leaf
            check = False #Set the boolean to include the event false
            if tract[5] != 'NS': #If the p-value is determined
                if float(tract[5]) < 1E-150: #Color differently depending on the p-value
                    color = '#944654' 
                    color2 = 'white'
                    check = True
                elif float(tract[5]) < 1E-100:
                    color = '#C46877' 
                    color2 = 'white'
                    check = True
                elif float(tract[5]) < 1E-50:
                    color = '#DB8484' 
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 1E-25:
                    color = '#F7B9B0'
                    color2 = 'black'
                    check = True
                elif float(tract[5]) < 1E-10:
                    color = '#F7CFC1' 
                    color2 = 'black'
                    check = True
                # elif float(tract[5]) < 0.05:
                #     color, color2 = 'linen', 'black'
                #     check = False
                else: 
                    # color, color2 = 'white', 'black'
                    check = False #Filter out tracts with p-values > 0.05
            if 'Unknown' in tract[1]: #If the minor parent includes "Unknown"
                tract[1] = tract[1].replace('Unknown (', '').replace(')', '?') #Replace it with a question mark
            # if 'Unknown' in tract[2]:
            #     tract[2] = tract[2].replace('Unknown (', '').replace(')', '?') #Replace it with a question mark
            tract[1] = tract[1].replace('Z12361', '') #Remove part of the name of DSM
            # tract[2] = tract[2].replace('Z12361', '') #Remove part of the name of DSM
                    
            if check: #If the tract should be included
                info_list = [int(tract[3]), int(tract[4]), '[]', 0, segment_height, 
                             'lightgrey', color, 
                             f'Verdana|16|{color2}|{tract[1]}'] #Create a list with tract information
                if info_list[1] > middle_point: #If the tract goes on the second segment
                    info_list[0] += 1 #Add one to the start (so that it doesn't overlap with the end of the previous segment)
                    info_list[1] += 1 #Add one to the end
                    
                motifs.append(info_list) #Add the list to a list of motifs to be shown in the figure
            
        tract_list = recomb_tracts[leaf.name] #Load the list of recombination tracts associated to the leaf
        
        to_add = [] #Create an empty list
        for i in range(len(motifs)): #Loop through the motifs to plot
            if motifs[i][1] - motifs[i][0] > 500: #If the distance between the motifs is at least 500 nt long
                if len(to_add) > 0: #If the to_add list is not empty
                    check = [] #Create an empty list to store booleans
                    for motif in to_add: #Loop through motifs that have been included in to_add
                        if (motifs[i][1] < motif[0]) or (motifs[i][0] > motif[1]): #If the motifs don't overlap
                            check.append(False) #Add a False value to the check list
                        else: check.append(True) #Otherwise, add a True check to the list
                        
                    if sum(check) != 0: #If there are overlapping segments
                        to_add1 = [mot for mot in to_add if mot[1] <= middle_point] #Retrieve the motifs on the first half of the segment
                        to_add2 = [] #Create an empty list to retrieve the motifs in the second half of the segment
                        for mot in to_add: #Loop through all the tracts in to_add
                            if mot[1] > middle_point: #If the end of the tract is after the middle point
                                mot[0] -= middle_point #Substract the middle point from the start position
                                mot[1] -= middle_point #Substract the middle point from the end position
                                to_add2.append(mot) #Append the motif to the list
                        seqFace = SeqMotifFace(seq, motifs = to_add1, #Create a SeqMotifFace to plot the tracts
                                               height = segment_height, seq_format = '[]', #Represent the segment as a horizontal bar of height 20
                                               gap_format = 'blank', #Show gaps as blank spaces
                                               scale_factor = segment_scale, #Scale the segment
                                               bgcolor = 'lightgrey', #Background color of the segment
                                               fgcolor = 'lightgrey') #Color of the border of the segment
                        (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #Add the SeqMotifFace to column 1 of the plot
                        
                        seqFace2 = SeqMotifFace(seq2, motifs = to_add2, #Create a SeqMotifFace to plot the tracts
                                               height = segment_height, seq_format = '[]', 
                                               gap_format = 'blank', 
                                               scale_factor = segment_scale, 
                                               bgcolor = 'lightgrey', 
                                               fgcolor = 'lightgrey')
                        (t & f'{leaf.name}').add_face(seqFace2, 3, 'aligned') #Add the SeqMotifFace to column 3 of the plot
                        
                        to_add = [motifs[i]] #Re-start the to_add list
                    else: to_add.append(motifs[i]) #Else, append the motif to the list
                else: to_add.append(motifs[i]) #Else (if the list is empty), append the motif to the list
                
        if to_add != []: #If the motif list is not empty
            to_add1 = [mot for mot in to_add if mot[1] <= middle_point] #Retrieve the motifs on the first half of the segment
            to_add2 = [] #Create an empty list to retrieve the motifs in the second half of the segment
            for mot in to_add: #Loop through all the tracts in to_add
                if mot[1] > middle_point: #If the end of the tract is after the middle point
                    mot[0] -= middle_point #Substract the middle point from the start position
                    mot[1] -= middle_point #Substract the middle point from the end position
                    to_add2.append(mot) #Append the motif to the list
                    
            seqFace = SeqMotifFace(seq, motifs = to_add1, height = segment_height, #Create a SeqMotifFace to plot the tract
                                   seq_format = '[]', gap_format = 'blank', 
                                   scale_factor = segment_scale, 
                                   bgcolor = 'lightgrey', 
                                   fgcolor = 'lightgrey')
            (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #Add the SeqMotifFace to column 1 of the plot
            
            seqFace2 = SeqMotifFace(seq2, motifs = to_add2, height = segment_height, #Create a SeqMotifFace to plot the tract
                                   seq_format = '[]', gap_format = 'blank', 
                                   scale_factor = segment_scale, 
                                   bgcolor = 'lightgrey', 
                                   fgcolor = 'lightgrey')
            (t & f'{leaf.name}').add_face(seqFace2, 3, 'aligned') #Add the SeqMotifFace to column 3 of the plot
            
        elif leaf.name == 'A1001': #If the leaf is the root leaf
            loci1 = [locus for locus in loci if locus[1] <= middle_point] #Retrieve the CDS on the first half of the segment
            loci2 = [] #Create an empty list to retrieve the CDS in the second half of the segment
            for locus in loci: #Loop through all the CDS
                if locus[1] > middle_point: #If the end of the CDS is after the middle point
                    locus[0] -= middle_point - len(miniseq) #Substract the middle point and the separator length from the start position
                    locus[1] -= middle_point - len(miniseq) #Substract the middle point and the separator length from the start position
                    loci2.append(locus)
                    
            seqFace = SeqMotifFace(seq, motifs = loci1, height = segment_height*2, 
                                    seq_format = 'blank', gap_format = 'blank', 
                                    scale_factor = segment_scale, 
                                    bgcolor = 'lightgrey', 
                                    fgcolor = 'lightgrey')
            (t & f'{leaf.name}').add_face(seqFace, 1, 'aligned') #Add the CDS to column 1 of the plot
            
            seqFace2 = SeqMotifFace(seq2, motifs = loci2, height = segment_height*2, 
                                    seq_format = 'blank', gap_format = 'blank', 
                                    scale_factor = segment_scale, 
                                    bgcolor = 'lightgrey', 
                                    fgcolor = 'lightgrey')
            (t & f'{leaf.name}').add_face(seqFace2, 3, 'aligned') #Add the CDS to column 3 of the plot
    
# =============================================================================
# 8. Apply the style to the tree and save it to a file
# =============================================================================

t.ladderize(1) #Reverse the order of the nodes
t.convert_to_ultrametric() #Make the tree ultrametric for nicer visualization
t.render(outplot, tree_style = ts) #Save the styled tree to a PNG file
t.render(outplot.replace('png', 'svg'), tree_style = ts) #Save the styled tree to a SVG file
t.render(outplot.replace('png', 'pdf'), tree_style = ts) #Save the styled tree to a PDF file
