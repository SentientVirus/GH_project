#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:28:42 2023

This is a script to use PyGenomeViz to plot the region of the genome where
the NGB gene is in the representative set of strains.

@author: Marina Mota Merlo
"""
#from Bio.Blast.Applications import NcbiblastpCommandline as cline_blast
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank as gbk_read
from matplotlib.lines import Line2D
from Bio import GenBank as gbk
from Bio.SeqFeature import SimpleLocation
import os
import pandas as pd

# TO DO:
    # Add meaningful comments
    # Check if there is a phage in H4B5-05J, and double-check all phage positions
    # Try to simplify the trimming of the tab_df
    
# =============================================================================
# 1. Define inputs
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
              'H4B508X': '#0072B2', 'MP2': '#33B18F', 'IBH001': '#D55E00', 
              'DSMZ12361': '#0072B2'}

#NGB locus tags
NGB = ['A0901_05380', 'A1003_04750', 'A1202_05530', 'A1401_04720', 'A1805_04920',
       'FHON2_04830', 'G0101_04800', 'H1B105A_04750', 'H1B302M_04960', 
       'H3B101A_04720', 'H3B104X_04750', 'H3B202X_04770', 'H3B203J_04720', 
       'H3B203M_04710', 'H3B206M_05060', 'H3B209X_04660', 'H4B111J_04920', 
       'H4B202J_04590', 'H4B204J_04710', 'H4B205J_04990', 'H4B211M_04810', 
       'H4B405J_04950', 'H4B406M_05270', 'H4B503X_04680', 'H4B504J_05510', 
       'H4B505J_04980', 'K2W83_RS02365']

#Locus tag numbers of phage genes, based on phage predictions (check predictions)
phage_pos = {'MP2': (190, 500), 
             'A1202': (4780, 5460), 
             'H4B5-04J': (4840, 5460), 
             'A0901': (4670, 5350)}

#Order in which strains and corresponding comparisons should be plotted
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


vmin = 60 #Minimum sequence identity

#Paths to inputs and outputs
projdir = os.path.expanduser('~') + '/GH_project'
outpath = f'{projdir}/blast_tabs'
folder_path2 = os.path.expanduser('~') + '/Akunkeei_files/gbff'
outfig = f'{projdir}/plots/pyGenomeViz/NGB_region.svg'


# =============================================================================
# 2. Establish the start and end positions for every strain and save them to a 
# dictionary. Get the chromosome IDs to change them to strain IDs.
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
    locus_dict = {'IBH001': 'LDX55_02300', 'DSMZ12361': 'K2W83_RS02275', 
                  'MP2': 'APS55_RS00555'}
    accession_strain = {} #Temporary variable to store NCBI accessions
    lengths_dict = {} #Dictionary to store the total genome length of each strain
    pos_dict = {} #Diccionary of start and end positions
    
    for strain in phylo_dict.values(): #Loop through strain names
        genbank_file = f'{folder_path2}/{strain}_genomic.gbff' #Set name of GenBank file
        start = 0 #Initialize start position
        with open(genbank_file) as handle: #Read GenBank file
            for record in gbk.parse(handle): #Loop through the records (chromosomes or plasmids) in the file
                if strain not in accession_strain.keys(): #If the strain is not in the accession dictionary (we only get the first record)
                    accession_strain[strain] = f'{record.accession[0]}.1' #Associate accession to strain name
                    genome_length = len(record.sequence) #Get length of the record
                    lengths_dict[strain] = genome_length #Initialize length dictionary

                if strain not in ['IBH001', 'MP2', 'DSMZ12361']: #Retrieve start position in the strains
                    for feature in record.features: #Loop through features (mostly CDS) in the record
                        gene_name = [qual for qual in feature.qualifiers if 'gene' in qual.key] #Get the four-letter gene name
                        if len(gene_name) > 0: #If there is at least one qualifier
                            gene = gene_name[0] #Set the gene name to the first qualifier (gene name or locus tag)
                            if 'pspB' in gene.value: #If the gene name includes pspB
                                start = int(feature.location.split('..')[0]) + 2000 #Set start location to the beginning of pspB + 2kb
                                if start < 4*10^5: #If the start position is too large
                                    start = 0 #Reset start position
                                else: #Otherwise, stop looking for the gene
                                    break
                            elif start == 0 and 'pspA' in gene.value: #If there is a gene annotated as pspA instead
                                start = int(feature.location.split('..')[0]) + 2000 #Set the start location to the beginning of pspA + 2kb
                                if strain == 'H1B3-02M': #For this particular strain, add 10kb to the start position
                                    start += 10000
                                break

                else: #If the strain is annotated differently
                    for feature in record.features: #Loop through genes in the strain
                        locus_tag = [qual for qual in feature.qualifiers if 'locus_tag' in qual.key] #Loop through gene information
                        if len(locus_tag) > 0 and locus_tag[0].value.replace('"', '') == locus_dict[strain]: #Get the gene with the same locus tag as in locus_dict
                            start = int(feature.location.split('..')[0].replace('complement(', '')) + 2000 #Set the start position to the start of that gene + 2kb

            if strain == 'MP2': #If the strain is MP2 (assembly of the reverse strand)
                start = lengths_dict[strain] - start #Set start to length minus previous start
            #For all strains    
            end = start + length #Set the end position as the start plus the length
            pos_dict[strain] = (start, end) #Save position tuple to dictionary
            
    return accession_strain, pos_dict, lengths_dict #Give three dictionaries as output
    
acc, pos, lengths = get_info(phylo_order, folder_path2, 40000) #Run function

# =============================================================================
# 3. Create tracks with the CDS to plot
# =============================================================================

#Set plot style
gv = GenomeViz(
    fig_track_height = 0.42, #Height of the tracks with the representation of the CDS
    link_track_ratio = 0.8, #Size ratio between the links (Blast comparisons) and the tracks
    track_align_type = 'center', #Align tracks to the center
    )

gv.set_scale_bar(scale_size_label=(5000, '5 Kb')) #Set a legend for the scale of the plot (5kb bar)

for i in range(1, len(phylo_order.keys())+1): #Loop through strain names in the order dictionary
    strain = phylo_order[i] #Get the name of a strain
    genbank_file = f'{folder_path2}/{strain}_genomic.gbff' #Get GenBank file based on strain name
    genbk = gbk_read(genbank_file) #Read GenBank file with Pygenomeviz
    chromosome_length = len(genbk.records[0].seq)
    if strain in ['A1202', 'MP2', 'H4B5-04J', 'A0901']: #If the strain has a phage
        segments = dict(region1=(pos[strain][0], pos[strain][0]+16000), region2=(pos[strain][1]+12000, pos[strain][1]+36000)) #Add two segments, skipping 36kb between them
    else:
        segments = dict(region1=(pos[strain][0], pos[strain][1])) #Otherwise, get the complete 40kb segment
    track = gv.add_feature_track(name = genbk.name.replace('_genomic', '').replace('12361', ''), #Create track (use DSMZ as strain name for DSMZ12361)
                                 segments = segments, #Add segments to track
                                 label_kws = dict(color = leaf_color[strain.replace('-', '')])) #Add strain name color based on phylogroup
    
    for segment in track.segments: #Loop through segments in the track
        if strain == 'MP2': #If the strain is MP2
            target_range = (lengths['MP2'] - segment.range[1], 
                            lengths['MP2'] - segment.range[0]) #Get the target region on the opposite strand
            features = genbk.extract_features(feature_type = 'CDS', 
                                              target_range = target_range) #Extract features from GenBank file
            for feature in features: #Loop through features
                new_end = lengths['MP2'] - feature.location.start #Reverse start
                new_start = lengths['MP2'] - feature.location.end #Reverse end
                new_strand = feature.location.strand*-1 #Reverse strand
                feature.location = SimpleLocation(new_start, new_end, 
                                                  new_strand) #Update feature position
        else: #If the strain is not MP2
            features = genbk.extract_features(feature_type = 'CDS', 
                                              target_range = segment.range)  #Extract the features directly, based on the segment range
                

        for cds in features: #Loop through CDS
            protstart = int(cds.location.start) #Get CDS start
            end = int(cds.location.end) #Get CDS end
            strand = cds.location.strand #Get strand
            color = '#E3DAC9' #Set color of most CDS
            gene_name  = '' #Initialize gene name
            if 'gene' in cds.qualifiers.keys():
                gene_name = cds.qualifiers['gene'][0] #When possible, add four-letter gene name
            if 'transposase' not in cds.qualifiers['product'][0]: #If the gene is not a transposon
                if cds.qualifiers['locus_tag'][0][3:] in NGB or cds.qualifiers['locus_tag'][0] in NGB: #Check if it is in the NGB list
                    gene_name = 'NGB' #If so, add NGB as gene name
                    color = '#336E74' #Change the color of the CDS (teal)
                #If the gene is a cas gene or endonuclease
                elif 'restriction' in cds.qualifiers['product'][0].lower() or 'endonuclease' in cds.qualifiers['product'][0].lower() or 'HsdR' in cds.qualifiers['product'][0] or gene_name == 'cas2':
                    color = '#BF6D58' #Change color (brown)
                elif 'ABC' in cds.qualifiers['product'][0]: #If the gene is annotated as ABC transporter
                    color = '#8FBF58' #Change color (green)
                elif strain in phage_pos.keys(): #If the strain has a phage
                    loc_num = int(cds.qualifiers['locus_tag'][0].split('_')[1].replace('RS', '')) #Get number of the locus tag
                    if (loc_num >= phage_pos[strain][0] and loc_num <= phage_pos[strain][1]) or (cds.qualifiers['locus_tag'][0] == 'APS55_RS08015'): #If the locus tag number corresponds to a phage
                        color = '#8F58BF' #Change color (purple)

            elif 'transposase' in cds.qualifiers['product'][0]: #If the CDS is a tranposon
                color = 'black' #Change color to black
             
            if segment.start <= protstart <= end <= segment.end: #If the CDS is inside the segment to be plotted
                segment.add_feature(protstart, end, strand, label = gene_name, #Add CDS to segment (position and label)
                                  plotstyle = 'bigarrow', fc = color, lw = 1, arrow_shaft_ratio = 1, #Set arrow style, CDS color, CDS line width and arrow vs shaft ratio
                                  text_kws = dict(color = 'black', rotation = 45, #Set label properties (label and rotation)
                                                  size = 15, ymargin = 0, #Set label properties (text size and distance from the CDS)
                                                  vpos = 'top', hpos = 'left')) #Set label properties (vertical and horizontal position)

        if strain != 'MP2':
            segment.add_sublabel(f'{segment.start:,} - {segment.end:,} bp')  #Add text indicating segment range to the plot
        else:
            segment.add_sublabel(f'{chromosome_length - segment.start:,} - {chromosome_length - segment.end:,} bp')  #Add text indicating segment range to the plot

        track.align_label = True #Align track label (strain name) to track
        track.set_segment_sep() #Set separator (//) between segments

            
# =============================================================================
# 4. Plot Blast comparisons
# =============================================================================
for i in range(1, len(phylo_order.keys())): #Loop through strains in the order in which they will be plotted    
    with open(f'{outpath}/NGB_{i}.tab') as tabfile: #Open Blast comparison file for each strain
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
    file_df = pd.read_csv(f'{outpath}/NGB_{i}.tab', sep = '\t', skiprows = 5, header = None)
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
    if phylo_order[i] in ['A1202', 'MP2', 'H4B5-04J', 'A0901']: #If the query strain has a phage
        subject_check = True #Set the two-segment check to True
        subject_pos = (subject_pos[0], subject_pos[1] + 36000) #Add 36kb to the end of the segment
    if phylo_order[i+1] in ['A1202', 'MP2', 'H4B5-04J', 'A0901']: #If the subject has a phage
        query_check = True #Set the two-segment check to True
        query_pos = (query_pos[0], query_pos[1] + 36000) #Add 36kb to the end of the segment
        
    print(strain1_name, strain2_name, query_check, subject_check)
    
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
    
    if query_check: #If there are two segments to plot for the query
        seg1_q = tab_df.copy() #Copy tab_df
        seg1_q = seg1_q.drop(seg1_q.index[seg1_q['q. end'] > query_pos[0] + 16000]) #Set the first segment to start + 16kb
        seg1_q = seg1_q.reset_index() #Reset the index of the dataframe
        seg2_q = tab_df.copy() #Copy tab_df
        seg2_q = seg2_q.drop(seg2_q.index[seg2_q['q. start'] < query_pos[1] - 24000]) #Set the second segment to end - 24kb
        seg2_q = seg2_q.reset_index() #Reset the index of the dataframe
    if subject_check: #If there are two segments to plot for the subject
        seg1_s = tab_df.copy() #Copy tab_df
        seg1_s = seg1_s.drop(seg1_s.index[seg1_s['s. start'] > subject_pos[0] + 16000]) #Set the first segment to start + 16kb
        seg1_s = seg1_s.drop(seg1_s.index[seg1_s['q. end'] > query_pos[0] + 16000])
        seg1_s = seg1_s.reset_index() #Reset the index of the dataframe
        seg2_s = tab_df.copy() #Copy tab_df
        seg2_s = seg2_s.drop(seg2_s.index[seg2_s['s. end'] < subject_pos[1] - 24000]) #Set the second segment to end - 24kb
        seg2_s.loc[seg2_s['s. start'] < subject_pos[1] - 24000, 's. start'] = subject_pos[1] - 24000
        if query_check:
            seg2_s.loc[seg2_s['q. start'] < query_pos[1] - 24000, 'q. start'] = query_pos[1] - 24000
        seg2_s = seg2_s.reset_index() #Reset the index of the dataframe
    
    # if strain1_name == 'A0901':
    #     raise ValueError('Stop here!!!')
        
    if subject_check == False and query_check == False: #If none of the strains in the comparison has multiple segments
        for j in range(len(tab_df)):
            #Get strain names from dataframe and shorten name of DSMZ
            strain1_name = tab_df.loc[j, 'query acc.ver'].replace('12361', '')
            strain2_name = tab_df.loc[j, 'subject acc.ver'].replace('12361', '')
            
            #Define the percentage of identity (to color the links)
            identity = tab_df.loc[j, '% identity']
            
            #Define the links (strain name, segment, start, end)
            link1 = (strain1_name, 'region1', tab_df.loc[j, 'q. start'], tab_df.loc[j, 'q. end'])
            link2 = (strain2_name, 'region1', tab_df.loc[j, 's. start'], tab_df.loc[j, 's. end'])
            
            #Add the links to the gv plot, especifying the colors of Blast matches (v indicates the variable setting color, and vmin is the minimum value)
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = vmin, curve = True, alpha = 0.7) #Curve makes the matches form curves
            
    elif query_check == True and subject_check == False: #If the query strain has two segments
        for j in range(len(seg1_q)): #Do two comparisons, one for each segment
            link1 = (seg1_q.loc[j, 'query acc.ver'], 'region1', seg1_q.loc[j, 'q. start'], seg1_q.loc[j, 'q. end'])
            link2 = (seg1_q.loc[j, 'subject acc.ver'], 'region1', seg1_q.loc[j, 's. start'], seg1_q.loc[j, 's. end'])
            identity = seg1_q.loc[j, '% identity']
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = vmin, curve = True, alpha = 0.7)
        for j in range(len(seg2_q)):
            link1 = (seg2_q.loc[j, 'query acc.ver'], 'region2', seg2_q.loc[j, 'q. start'], seg2_q.loc[j, 'q. end'])
            link2 = (seg2_q.loc[j, 'subject acc.ver'], 'region1', seg2_q.loc[j, 's. start'], seg2_q.loc[j, 's. end'])
            identity = seg2_q.loc[j, '% identity']
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = vmin, curve = True, alpha = 0.7)
    
    elif subject_check == True and query_check == False: #If the subject strain has two segments
        for j in range(len(seg1_s)): #Do two comparisons, one for each segment
            link1 = (seg1_s.loc[j, 'query acc.ver'], 'region1', seg1_s.loc[j, 'q. start'], seg1_s.loc[j, 'q. end'])
            link2 = (seg1_s.loc[j, 'subject acc.ver'], 'region1', seg1_s.loc[j, 's. start'], seg1_s.loc[j, 's. end'])
            identity = seg1_s.loc[j, '% identity']
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = vmin, curve = True, alpha = 0.7)
        for j in range(len(seg2_s)):
            link1 = (seg2_s.loc[j, 'query acc.ver'], 'region1', seg2_s.loc[j, 'q. start'], seg2_s.loc[j, 'q. end'])
            link2 = (seg2_s.loc[j, 'subject acc.ver'], 'region2', seg2_s.loc[j, 's. start'], seg2_s.loc[j, 's. end'])
            identity = seg2_s.loc[j, '% identity']
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = vmin, curve = True, alpha = 0.7)
            
    elif subject_check == True and query_check == True:
        for j in range(len(seg1_s)): #Do two comparisons, one for each segment
            link1 = (seg1_s.loc[j, 'query acc.ver'], 'region1', seg1_s.loc[j, 'q. start'], seg1_s.loc[j, 'q. end'])
            link2 = (seg1_s.loc[j, 'subject acc.ver'], 'region1', seg1_s.loc[j, 's. start'], seg1_s.loc[j, 's. end'])
            identity = seg1_s.loc[j, '% identity']
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = vmin, curve = True, alpha = 0.7)
        for j in range(len(seg2_s)):
            link1 = (seg2_s.loc[j, 'query acc.ver'], 'region2', seg2_s.loc[j, 'q. start'], seg2_s.loc[j, 'q. end'])
            link2 = (seg2_s.loc[j, 'subject acc.ver'], 'region2', seg2_s.loc[j, 's. start'], seg2_s.loc[j, 's. end'])
            identity = seg2_s.loc[j, '% identity']
            gv.add_link(link1, link2, color = 'grey', inverted_color = 'red', 
                        v = identity, vmin = vmin, curve = True, alpha = 0.7)

#Set legend for the Blast matches (two adjacent colorbars that show the colors of forward and reverse matches with a minimum of 50% identity)
gv.set_colorbar(['grey', 'red'], vmin = vmin, bar_height = 0.05, 
                tick_labelsize = 16, alpha = 0.7) #Bar height sets the height of the bars and tick_labelsize sets the font size of tick labels
 
fig = gv.plotfig(dpi = 400) #Plot the figure with a resolution of 400 dpi

#Create legend text (label) and icons (marker sets the size, color sets marker color, ms sets marker size and ls sets border size)
#The first two variables are the data to be plotted in the x and y axis, and therefore are empty
handles = [
    Line2D([], [], marker='', color='black', label='Tracks', ms=20, ls='none'),
    Line2D([], [], marker='>', color='#E3DAC9', label='CDS', ms=20, ls='none'),
    Line2D([], [], marker='>', color='#336E74', label='NGB glycosyl hydrolase', ms=20, ls='none'),
    Line2D([], [], marker='>', color='#8F58BF', label='Phage genes', ms=20, ls='none'),
    Line2D([], [], marker='>', color='#8FBF58', label='ABC transporters', ms=20, ls='none'),
    Line2D([], [], marker='>', color='#BF6D58', label='CRISPR endonucleases', ms=20, ls='none'),
    Line2D([], [], marker='>', color='black', label='Transposase', ms=20, ls='none'),
    Line2D([], [], marker='', color='black', label='', ms=20, ls='none'),
    Line2D([], [], marker='', color='black', label='Matches', ms=20, ls='none'),
    Line2D([], [], marker='s', color='grey', label='Forward match', ms=20, ls='none'),
    Line2D([], [], marker='s', color='red', label='Reverse match', ms=20, ls='none'),
    Line2D([], [], marker='', color='#771853', label='', ms=20, ls='none'),
    Line2D([], [], marker='', color='#771853', label='Strain names', lw=2, ms=20, ls='none'),
    Line2D([], [], marker='X', color='#0072B2', label='Phylogroup A', ms=15, ls='none'),
    Line2D([], [], marker='X', color='#33B18F', label='Phylogroup B', ms=15, ls='none'),
    Line2D([], [], marker='X', color='#D55E00', label='Phylogroup C', ms=15, ls='none'),
    Line2D([], [], marker='X', color='#FF74D6', label='Phylogroup E', ms=15, ls='none'),
    Line2D([], [], marker='X', color='#771853', label='Phylogroup F', ms=15, ls='none'),
    Line2D([], [], marker='X', color='black', label='Not in Dyrhage et al. (2022)', ms=15, ls='none')
]

#Assign legent to figure. bbox_to_anchor sets the position, frameon removes the 
#frame (border) of the legend box, and labelspacing increases vertical space between legends
legend = fig.legend(handles=handles, bbox_to_anchor=(1.28, 1), frameon = False,
                    fontsize = 16)

#Set legend headers to bold
legend.get_texts()[0].set_fontweight('bold')
legend.get_texts()[8].set_fontweight('bold')
legend.get_texts()[12].set_fontweight('bold')
            
fig.savefig(outfig) #Save figure to SVG
fig.savefig(outfig.replace('svg', 'tiff')) #Save figure to TIFF
fig.savefig(outfig.replace('svg', 'png')) #Save figure to PNG
fig.savefig(outfig.replace('svg', 'pdf')) #Save figure to PDF
gv.savefig_html(outfig.replace('svg', 'html')) #Save figure to HTML
    
