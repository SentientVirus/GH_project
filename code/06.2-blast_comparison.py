#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:28:42 2023

This is the code to generate a Blast comparison plot in PygenomeViz focusing
on the region where the GH70 genes GS1, GS2 and BRS are, together with the
GH32 S1-3 genes.

@author: Marina Mota Merlo
"""

from Bio.Blast.Applications import NcbiblastpCommandline as cline_blast
from matplotlib.lines import Line2D
from Bio.SeqFeature import SimpleLocation
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank as gbk_read
from Bio import GenBank as gbk
import os
import pandas as pd

# =============================================================================
# Color dictionary for strain names
# =============================================================================
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
              'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': '#D55E00', 
              'DSMZ12361': '#0072B2'}

# =============================================================================
# Locus tags
# =============================================================================
GS1 = ['A1001_12310', 'A1202_13520', 'A1401_12750', 'A1805_12820',
       'FHON2_13540', 'G0101_12800', 'G0403_13100', 'H1B104J_13010', 
       'H1B105A_12300', 'H1B302M_12880', 'H3B101A_13240', 'H3B104J_12990', 
       'H3B104X_13200', 'H3B202X_12850', 'H3B203J_13370', 'H3B203M_12480', 
       'H3B206M_12840', 'H3B209X_13340', 'H4B111J_13570', 'H4B202J_12880', 
       'H4B204J_13330', 'H4B205J_12990', 'H4B211M_13000', 'H4B402J_12600', 
       'H4B405J_13360', 'H4B406M_13460', 'H4B412M_13240', 'H4B501J_12890', 
       'H4B503X_12670', 'H4B504J_13460', 'H4B505J_12880', 'APS55_RS03850', 
       'LDX55_06325', 'K2W83_RS06180']
GS2 = ['FHON2_13560', 'FHON2_13570', 'G0403_13120', 'H1B302M_12900', 
       'H3B101A_13260', 'H3B104X_13220', 'H3B202X_12860', 'H3B203J_13390',
       'H3B209X_13370', 'H4B504J_13480', 'H4B505J_12900', 'APS55_RS03845']
GS3 = ['H3B206M_12830', 'H4B111J_13560', 'H4B405J_13350', 'H4B406M_13450']
GS4 = ['A1003_12540']
GS2BRS = ['H3B104J_13020', 'H4B202J_12890', 'H4B204J_13350', 'LDX55_06335', 
          'K2W83_RS06185']
BRS = ['A1401_12760', 'FHON2_13550', 'G0403_13110', 'H1B302M_12890', 
       'H3B101A_13250', 'H3B104J_13000', 'H3B104J_13010', 'H3B104X_13210',
       'H3B203J_13380', 'H3B209X_13350', 'H4B204J_13340', 'H4B504J_13470',
       'H4B505J_12890', 'LDX55_06330']
S2a = ['A0901_13270', 'A1001_12300', 'A1805_12810', 'H1B104J_13000', 
       'H3B203M_12470', 'H4B206J_13400', 'H4B402J_12590', 'H4B412M_13230',
       'H4B503X_12660']
S2b = ['A1805_12800', 'H1B104J_12990', 'H4B412M_13220']
S1 = ['G0101_12790', 'H4B205J_12980', 'H4B211M_12990', 'H4B501J_12880', 
      'K2W83_RS06175']
S3 = ['A0901_13360', 'A1001_12360', 'A1404_13450', 'H3B203M_12520', 
      'H4B206J_13490', 'H4B402J_12660', 'H4B412M_13310', 'H4B503X_12760']

# =============================================================================
# In this section, we run Blast between pairs of strains following the 
# phylogeny.
# =============================================================================

# phylo_order = {1:'G0101', 2:'H4B2-05J', 3:'H4B2-11M', 4:'H1B1-04J', 
#               5:'H4B4-12M', 6:'H4B5-01J', 7:'H3B1-01A', 8:'H3B2-03J', 
#               9:'H3B1-04X', 10:'H3B2-02X', 11:'A1003', 12:'H4B2-04J', 
#               13:'H4B5-03X', 14:'H4B4-02J', 15:'H1B3-02M', 16:'Fhon2', 
#               17:'DSMZ12361', 18:'H4B2-02J', 19:'H1B1-05A', 20:'H3B1-04J', 
#               21:'A1401', 22:'A1202', 23:'H4B5-05J', 24:'H3B2-09X', 
#               25:'H4B5-04J', 26:'A1805', 27:'MP2', 28:'G0403', 29:'H3B2-06M', 
#               30:'H4B4-05J', 31:'H4B1-11J', 32:'H4B4-06M', 33:'A0901', 
#               34:'IBH001', 35:'H4B2-06J', 36:'H3B2-03M', 37:'A1404', 38:'A1001'}

# phylo_order = {1:'G0101', 2:'H4B2-05J', 3:'H4B2-11M', 4:'H4B5-01J', 
#               5:'H4B4-12M', 6:'H1B1-04J', 7:'H3B1-01A', 8:'H3B2-03J', 
#               9:'H4B2-04J', 10:'H3B1-04X', 11:'H3B2-02X', 12:'A1003',
#               13:'H1B3-02M', 14:'H4B4-02J', 15:'H4B5-03X', 16:'Fhon2', 
#               17:'H1B1-05A', 18:'H4B2-02J', 19:'DSMZ12361', 20:'H3B1-04J', 
#               21:'A1202', 22:'A1401', 23:'H4B5-05J', 24:'H3B2-09X', 
#               25:'H4B5-04J', 26:'A1805', 27:'MP2', 28:'G0403', 29:'H4B4-05J', 
#               30:'H3B2-06M', 31:'H4B1-11J', 32:'H4B4-06M', 33:'IBH001', 
#               34:'H4B2-06J', 35:'A0901', 36:'H3B2-03M', 37:'A1404', 38:'A1001'}

phylo_order = {1:'H4B4-12M', 2:'H4B5-01J', 3:'G0101', 4:'H4B2-05J', 
              5:'H4B2-11M', 6:'H1B1-04J', 7:'H4B5-03X', 8:'H4B4-02J', 
              9:'H1B3-02M', 10:'H3B2-02X', 11:'H3B1-04X', 12:'A1003',
              13:'H3B1-01A', 14:'H3B2-03J', 15:'H4B2-04J', 16:'Fhon2', 
              17:'H4B2-02J', 18:'H1B1-05A', 19:'DSMZ12361', 20:'H3B1-04J', 
              21:'H3B2-06M', 22:'H4B4-05J', 23:'H4B1-11J', 24:'H4B4-06M', 
              25:'H3B2-03M', 26:'A0901', 27:'H4B2-06J', 28:'IBH001', 
              29:'H3B2-09X', 30:'H4B5-05J', 31:'A1401', 32:'A1202', 
              33:'H4B5-04J', 34:'MP2', 35:'A1805', 36:'G0403', 
              37:'A1404', 38:'A1001'}


color_dict = {'BRS': '#84E8A7', 'GS2BRS': '#BDE384', 'GS2': '#FFE570', 
              'GS3': '#FF986F', 'GS1': '#FF707C', 'GS4': '#FFC36F',
              'S1': '#656ED4', 'S2a': '#C870FF', 'S2b': '#F87BFF', 
              'S3': '#336E74'}


projdir = os.path.expanduser('~') + '/GH_project'
outpath = f'{projdir}/blast_tabs'
folder_path = os.path.expanduser('~') + '/Akunkeei_files/fna/reverse'
folder_path2 = os.path.expanduser('~') + '/Akunkeei_files/gbff/modified_gbff'
outfig = f'{projdir}/plots/pyGenomeViz/dynamic_region.svg'

if not os.path.exists(os.path.dirname(outfig)):
    os.makedirs(os.path.dirname(outfig))


def get_tabs(phylo_dict, folder_path, outpath):
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for i in range(1, len(phylo_dict)):
        strain1_file = f'{folder_path}/{phylo_dict[i]}_genomic.fna'
        strain2_file = f'{folder_path}/{phylo_dict[i+1]}_genomic.fna'
        outfile = f'{outpath}/{i}.tab'
        cline_input = cline_blast(cmd = 'blastn', query = strain2_file, 
                                  subject = strain1_file, remote = False, 
                                  out = outfile, outfmt = 7)
        os.system(str(cline_input));
        
get_tabs(phylo_order, folder_path, outpath)

# =============================================================================
# In this section, we establish the start and end positions for every strain
# and save them to a dictionary. I also get the chromosome IDs to change them
# to strain IDs.
# =============================================================================

def get_info(phylo_dict, folder_path2, length = 0):
    '''
    Function to get the start and end coordinates to be plotted.
    
    Input parameters:
        @param phylo_dict (dict, {str: int}): Dictionary indicating the order 
        in which the strains should be plotted. This parameter should 
        correspond to the order for the Blast comparisons.
        @param folder_path2 (str): Path where the input GenBank files are.
        @param length (int): Length of the segment to be plotted.
    
    Output parameters:
        @param accession_strain (dict, {str: [str, ...]}): Dictionary with 
        strain names as keys and NCBI accessions as values.
        @param pos_dict (dict, {str: (int, int)}): Dictionary with tuples
        indicating start and end positions.
        @lengths_dict (dict, {str: int}): Dictionary with the names of the strains
        as keys and the lengths of the chromosome as values.
    '''
    #In these strains, the target gene to set the start is annotated differently,
    #so it is added manually
    locus_dict = {'IBH001': 'LDX55_06270', 'DSMZ12361': 'K2W83_RS06135', 
                  'MP2': 'APS55_RS03895'}
    accession_strain = {} #Temporary variable to store NCBI accessions
    pos_dict = {} #Dictionary to store the total genome length of each strain
    lengths_dict = {} #Diccionary of start and end positions
    for strain in phylo_dict.values(): #Loop through strain names
        genbank_file = f'{folder_path2}/{strain}_genomic.gbff' #Set name of GenBank file
        with open(genbank_file) as handle: #Read GenBank file
            for record in gbk.parse(handle): #Loop through the records (chromosomes or plasmids) in the file
                if strain not in accession_strain.keys(): #If the strain is not in the accession dictionary (we only get the first record)
                    accession_strain[strain] = f'{record.accession[0]}.1' #Associate accession to strain name
                    genome_length = len(record.sequence) #Get length of the record
                    lengths_dict[strain] = genome_length #Initialize length dictionary
                    
                if strain not in ['IBH001', 'MP2', 'DSMZ12361']: #Retrieve start position in the strains
                    for feature in record.features:  #Loop through features (mostly CDS) in the record
                        gene_name = [qual for qual in feature.qualifiers if 'gene' in qual.key] #Get the four-letter gene name
                        if len(gene_name) > 0: #If there is at least one qualifier
                            gene = gene_name[0] #Set the gene name to the first qualifier (gene name or locus tag)
                            if 'ohrR' in gene.value: #If the gene name includes ohrR
                                start = int(feature.location.split('..')[0]) + 2000 #Set start location to the beginning of ohrR + 2kb

                else: #If the strain is annotated differently
                    for feature in record.features: #Loop through genes in the strain
                        locus_tag = [qual for qual in feature.qualifiers if 'locus_tag' in qual.key] #Loop through gene information
                        if len(locus_tag) > 0 and locus_tag[0].value.replace('"', '') == locus_dict[strain]: #Get the gene with the same locus tag as in locus_dict
                            start = int(feature.location.split('..')[0].replace('complement(', '')) + 5000 #Set the start position to the start of that gene + 5kb
                            if strain == 'MP2': #If the strain is MP2
                                start -= 10000 #Set the start to minus 10kb (forward)
            if strain != 'MP2': #If the strain is not MP2 (assembly of the reverse strand)
                start = lengths_dict[strain] - (start + length) #Set the start to the opposite strand
                end = start + length #Set the end to the start + segment length
                pos_dict[strain] = (start, end) #Save the new position to a dictionary
            else: #If the strain is MP2
                end = start - length #Calculate the beginning point of the segment
                pos_dict[strain] = (end, start) #Store the results in a dictionary
    return accession_strain, pos_dict, lengths_dict
    
acc, pos, lengths = get_info(phylo_order, folder_path2, 40000)


# =============================================================================
# Now, ready for the plotting!
# =============================================================================

# Set plot style
gv = GenomeViz(
    fig_track_height = 0.42, #Height of the tracks with the representation of the CDS
    link_track_ratio = 0.8, #Size ratio between the links (Blast comparisons) and the tracks
    track_align_type = 'center', #Align tracks to the center
    )

gv.set_scale_bar(scale_size_label=(5000, '5 Kb')) #Set a legend for the scale of the plot (5kb bar)

for i in range(1, len(phylo_order.keys())+1): #Loop through strain names in the order dictionary
    strain = phylo_order[i] #Get the name of a strain
    genbank_file = f'{folder_path2}/{strain}_genomic.gbff' #Get GenBank file based on strain name
    genbk = gbk_read(genbank_file) #Read GenBank file with PyGenomeViz
    segments = dict(region1=(pos[strain][0], pos[strain][1])) #Retrieve the segment to be plotted
    track = gv.add_feature_track(name = genbk.name.replace('_genomic', '').replace('Z12361', ''), #Create track (use DSM as strain name for DSMZ12361)
                                 segments = segments, #Add segments to track
                                 label_kws = dict(color = leaf_color[strain.replace('-', '')])) #Add strain name color based on phylogroup
    for segment in track.segments: #Loop through segments in the track
        if strain != 'MP2': #If the strain is not MP2
            target_range = (lengths[strain] - segment.range[1], 
                            lengths[strain] - segment.range[0]) #Get the target region on the opposite strand
            features = genbk.extract_features(feature_type = 'CDS', 
                                              target_range = target_range) #Extract features from GenBank file

            for feature in features: #Loop through features
                new_end = lengths[strain] - feature.location.start #Reverse start
                new_start = lengths[strain] - feature.location.end #Reverse end
                print(new_start, new_end)
                new_strand = feature.location.strand*-1 #Reverse strand
                feature.location = SimpleLocation(new_start, new_end, 
                                                  new_strand) #Update feature position
        else: #If the strain is MP2
            features = genbk.extract_features(feature_type = 'CDS', 
                                              target_range = segment.range)  #Extract the features directly, based on the segment range
        #Loop through CDS
        for cds in features: #Loop through CDS
            protstart = int(cds.location.start) #Get CDS start
            end = int(cds.location.end) #Get CDS end
            strand = cds.location.strand #Get strand
            color = 'skyblue' #Set color of most CDS
            gene_name  = '' #Initialize gene name
            if 'gene' in cds.qualifiers.keys():
                gene_name = cds.qualifiers['gene'][0] #When possible, add four-letter gene name
                if 'transposase' not in cds.qualifiers['product'][0]: #If the gene is not a transposon
                    gene_name = cds.qualifiers['gene'][0]
                    if '_partial' in gene_name: #Add stars to indicate that genes are incomplete
                        gene_name = gene_name.replace('_partial', '*')
                    elif 'gene' in cds.qualifiers.keys() and '_I' in gene_name: #Retrieve the set of genes that were manually annotated in the GenBanks
                        gene_name = gene_name.replace('_I', '')
                    for key in color_dict.keys(): #Use a different color for wach type of GH gene
                        if gene_name == key or gene_name[:-1] == key:
                            color = color_dict[key]
                else: #If the gene is a transposon
                    color = 'black' #Color it in black
                    
            elif 'transposase' in cds.qualifiers['product'][0]: #Color transposases in black
                gene_name = ''
                color = 'black'
             
            if segment.start <= protstart <= end <= segment.end: #If the CDS is inside the segment to be plotted
                segment.add_feature(protstart, end, strand, label = gene_name, #Add CDS to segment (position and label)
                                  plotstyle = 'bigarrow', fc = color, lw = 1, arrow_shaft_ratio = 1, #Set arrow style, CDS color, CDS line width and arrow vs shaft ratio
                                  text_kws = dict(color = 'black', rotation = 45, #Set label properties (label and rotation)
                                                  size = 15, ymargin = 0, #Set label properties (text size and distance from the CDS)
                                                  vpos = 'top', hpos = 'left')) #Set label properties (vertical and horizontal position)
        segment.add_sublabel(f'{segment.start:,} - {segment.end:,} bp')  #Add text indicating segment range to the plot

        track.align_label = True #Align track label (strain name) to track
        track.set_segment_sep() #Set separator (//) between segments

        
# =============================================================================
# Here I modify sligthly the tab files that I created.
# =============================================================================
for i in range(1, len(phylo_order.keys())): #Loop through strains in the order in which they will be plotted     
    with open(f'{outpath}/{i}.tab') as tabfile: #Open Blast comparison file for each strain
        k = 0 #Initialize parameter to retrieve column names
        for line in tabfile: #Loop through lines in tab file
                k += 1 #Increase the value of k
                #Generate a list with the fields:
                if 'Fields' in line: #If the line contains the string Fields
                    headers = line #Assign the header string to a variable
                elif k > 10: #If we pass the line
                    break #End loop
    
    header_list = headers.split(',') #Divide comma-separated fields in the header string to create a list
    header_list[0] = header_list[0].replace('# Fields: ', '') #Remove the start of the line
    header_list = [header.strip().strip('\n') for header in header_list] #Remove spaces and line breaks from beginning
    
    #Read Blast comparisons as tab-separated file, skip 5 rows and don't set the remaining rows as headers
    file_df = pd.read_csv(f'{outpath}/{i}.tab', sep = '\t', skiprows = 5, header = None)
    file_df = file_df[:-1] #Remove the last row, that also contains a comment
    file_df.columns = header_list #Add column names
    strain2 = file_df.values[0, 0] #Retrieve the accession of the query strain
    strain2_name = list(acc.keys())[list(acc.values()).index(strain2)] #Using the accession, retrieve the key (strain name)
    file_df = file_df.replace(strain2, strain2_name) #Replace the accession with the strain name in the dataframe
    strain1 = file_df.values[0, 1] #Do the same for the subject strain
    strain1_name = list(acc.keys())[list(acc.values()).index(strain1)]
    file_df = file_df.replace(strain1, strain1_name)
    print(f'Loading comparison between {strain1_name} and {strain2_name}') #Print information about the comparison being made

    query_pos = pos[phylo_order[i + 1]] #Retrieve start and end position of the query segment
    subject_pos = pos[phylo_order[i]] #Retrieve start and end position of the subject segment
    
    tab_df = file_df.copy() #Copy file_df to a new dataframe
    subject_check = False #Check that states whether the subject has one or two segments to plot
    query_check = False #Check that states whether the query has one or two segments to plot
    
    # Figure out how to loop through segments to drop genes (idea: compare strain name with tab df and retrieve those tracks, then loop through segments)
    tab_df = tab_df.drop(tab_df.index[tab_df['q. end'] < query_pos[0]]) #Drop CDS that start after the end of the query segment
    tab_df = tab_df.drop(tab_df.index[tab_df['q. start'] > query_pos[1]]) #Drop CDS that end before the start of the query segment
    tab_df.loc[tab_df['q. start'] < query_pos[0], 'q. start'] = query_pos[0] #Trim the start of matches that overlap with the segment, but start before the segment
    tab_df.loc[tab_df['q. end'] > query_pos[1], 'q. end'] = query_pos[1] #Trim the end of matches that overlap with the segment, but end after the segment
    
    #Do the same for the subject as for the query
    tab_df = tab_df.drop(tab_df.index[tab_df['s. end'] < subject_pos[0]])
    tab_df = tab_df.drop(tab_df.index[tab_df['s. start'] > subject_pos[1]])
    tab_df.loc[tab_df['s. end'] > subject_pos[1], 's. end'] = subject_pos[1]
    tab_df.loc[tab_df['s. start'] < subject_pos[0], 's. start'] = subject_pos[0]
    
    tab_df = tab_df[~tab_df['query acc.ver'].str.startswith('#')] #Remove rows with query strain names starting with #
    tab_df = tab_df.reset_index() #Reset the index of the dataframe
    
    for j in range(len(tab_df)):
        #Get strain names from dataframe and shorten name of DSM
        strain1_name = tab_df.loc[j, 'query acc.ver'].replace('Z12361', '')
        strain2_name = tab_df.loc[j, 'subject acc.ver'].replace('Z12361', '')
        
        #Define the percentage of identity (to color the links)
        identity = tab_df.loc[j, '% identity']
        
        #Define the links (strain name, segment, start, end)
        link1 = (strain1_name, 'region1', tab_df.loc[j, 'q. start'], tab_df.loc[j, 'q. end'])
        link2 = (strain2_name, 'region1', tab_df.loc[j, 's. start'], tab_df.loc[j, 's. end'])
        
        #Add the links to the gv plot, especifying the colors of Blast matches (v indicates the variable setting color, and vmin is the minimum value)
        gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                    v = identity, vmin = 60, curve = True, alpha = 0.7) #Curve makes the matches form curves
            
#Set legend for the Blast matches (two adjacent colorbars that show the colors of forward and reverse matches with a minimum of 50% identity)        
gv.set_colorbar(['grey', 'red'], vmin = 60, bar_height = 0.05, 
                tick_labelsize = 16, alpha = 0.7)
fig = gv.plotfig(dpi = 400)

#Create legend text (label) and icons (marker sets the size, color sets marker color, ms sets marker size and ls sets border size)
#The first two variables are the data to be plotted in the x and y axis, and therefore are empty
handles = [
    Line2D([], [], marker="", color='black', label="Tracks", ms=20, ls="none"),
    Line2D([], [], marker=">", color="skyblue", label="CDS", ms=20, ls="none")
    ]

handles += [Line2D([], [], marker=">", color=color_dict[GH_type], label=f'{GH_type} glycosyl hydrolase', ms=20, ls='none') for GH_type in color_dict.keys()]

handles += [
    Line2D([], [], marker=">", color='black', label="Transposase", ms=20, ls="none"),
    Line2D([], [], marker="", color='#771853', label="", ms=20, ls="none"),
    Line2D([], [], marker="", color='black', label="Matches", ms=20, ls="none"),
    Line2D([], [], marker="s", color='grey', label="Forward match", ms=20, ls="none"),
    Line2D([], [], marker="s", color='red', label="Reverse match", ms=20, ls="none"),
    Line2D([], [], marker="", color='#771853', label="", ms=20, ls="none"),
    Line2D([], [], marker="", color='#771853', label="Strain names", lw=2, ms=20, ls="none"),
    Line2D([], [], marker="X", color='#0072B2', label="Phylogroup A", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#33B18F', label="Phylogroup B", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#D55E00', label="Phylogroup C", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#FF74D6', label="Phylogroup E", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#771853', label="Phylogroup F", ms=15, ls="none"),
    Line2D([], [], marker="X", color='black', label="Not in Dyrhage et al. (2022)", ms=15, ls="none")
]

#Assign legent to figure. bbox_to_anchor sets the position, frameon removes the 
#frame (border) of the legend box, and labelspacing increases vertical space between legends
legend = fig.legend(handles=handles, bbox_to_anchor=(1.35, 1), frameon = False,
                    fontsize = 16)

#Set legend headers to bold
legend.get_texts()[0].set_fontweight('bold')
legend.get_texts()[14].set_fontweight('bold')
legend.get_texts()[18].set_fontweight('bold')
    
fig.savefig(outfig) #Save figure to SVG
fig.savefig(outfig.replace('svg', 'png')) #Save figure to PNG
fig.savefig(outfig.replace('svg', 'pdf')) #Save figure to PDF
fig.savefig(outfig.replace('svg', 'tiff')) #Save figure to TIFF
gv.savefig_html(outfig.replace('svg', 'html')) #Save figure to HTML
