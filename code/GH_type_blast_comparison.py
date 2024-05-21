#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:28:42 2023

This is the code to generate a Blast comparison plot in PygenomeViz focusing
on the region where the GH70 genes GS1, GS2 and BRS are, together with the
GH32 S1-3 genes. However, this script compares only strains that are in the
GS1 repset tree to each other, strains in the GS2 tree to each other and
strains in the BRS tree to each other.

@author: Marina Mota Merlo
"""
from Bio.Blast.Applications import NcbiblastpCommandline as cline_blast
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.transforms as transforms
from pygenomeviz import Genbank as gbk_read, GenomeViz
from Bio import GenBank as gbk
import os
import pandas as pd
from ete3 import Tree
from itertools import compress

# =============================================================================
# Color dictionary for strain names
# =============================================================================
leaf_color = {'A0901': '#A027FF', 'A1001': '#E5BA60', 'A1002': '#A027FF',
              'A1003': '#1E55F6', 'A1201': '#A027FF', 'A1202': '#00AEFF',
              'A1401': '#00AEFF', 'A1404': '#FF74D6', 'A1802': '#A027FF',
              'A1803': '#A027FF', 'A1805': '#00AEFF', 'A2001': '#A027FF',
              'A2002': '#E5BA60', 'A2003': '#E5BA60', 'A2101': '#00AEFF',
              'A2102': '#A027FF', 'A2103': '#1E55F6', 'G0101': '#1E55F6', 
              'G0102': '#1E55F6', 'G0103': '#1E55F6', 'G0401': '#1E55F6',
              'G0402': '#1E55F6', 'G0403': '#00AEFF', 'G0404': '#1E55F6',
              'G0405': '#1E55F6', 'G0406': '#00AEFF', 'G0407': '#1E55F6',
              'G0408': '#1E55F6', 'G0410': '#1E55F6', 'G0412': '#1E55F6',
              'G0414': '#1E55F6', 'G0415': '#1E55F6', 'G0417': '#1E55F6',
              'G0420': '#00AEFF', 'G0601': '#1E55F6', 'G0602': '#1E55F6',
              'G0702': '#1E55F6', 'G0801': '#1E55F6', 'G0802': '#1E55F6',
              'G0803': '#1E55F6', 'G0804': '#1E55F6', 'Fhon2': '#1E55F6',
              'H1B104J': '#1E55F6', 'H1B105A': '#1E55F6', 'H1B302M': '#1E55F6',
              'H2B105J': '#A027FF', 'H3B101A': '#1E55F6', 'H3B101J': '#1E55F6',
              'H3B101X': '#A027FF', 'H3B102A': '#1E55F6', 'H3B102X': '#A027FF',
              'H3B103J': '#1E55F6', 'H3B103X': '#A027FF','H3B103M': '#1E55F6', 
              'H3B104J': '#1E55F6', 'H3B104X': '#1E55F6', 'H3B107A': '#1E55F6', 
              'H3B109M': '#1E55F6', 'H3B110M': '#1E55F6', 'H3B111A': '#1E55F6', 
              'H3B111M': '#1E55F6', 'H3B202M': '#A027FF', 'H3B202X': '#1E55F6', 
              'H3B203J': '#1E55F6', 'H3B203M': '#A027FF', 'H3B204J': '#1E55F6', 
              'H3B204M': '#A027FF', 'H3B205J': '#1E55F6', 'H3B206M': '#A027FF', 
              'H3B207X': '#1E55F6', 'H3B208X': '#1E55F6', 'H3B209X': '#00AEFF', 
              'H4B101A': '#A027FF', 'H4B102A': '#A027FF', 'H4B103J': '#A027FF', 
              'H4B104A': '#A027FF', 'H4B111J': '#A027FF', 'H4B114J': '#A027FF', 
              'H4B116J': '#A027FF', 'H4B202J': '#1E55F6', 'H4B203M': '#A027FF', 
              'H4B204J': '#1E55F6', 'H4B205J': '#1E55F6', 'H4B206J': '#A027FF', 
              'H4B210M': '#A027FF', 'H4B211M': '#1E55F6', 'H4B303J': '#A027FF', 
              'H4B402J': '#1E55F6', 'H4B403J': '#A027FF', 'H4B404J': '#A027FF', 
              'H4B405J': '#A027FF', 'H4B406M': '#A027FF', 'H4B410M': '#A027FF', 
              'H4B411M': '#A027FF', 'H4B412M': '#1E55F6', 'H4B501J': '#1E55F6', 
              'H4B502X': '#1E55F6', 'H4B503X': '#1E55F6', 'H4B504J': '#00AEFF',
              'H4B505J': '#00AEFF', 'H4B507J': '#1E55F6', 'H4B507X': '#1E55F6', 
              'H4B508X': '#1E55F6', 'MP2': '#00AEFF', 'IBH001': 'black', 
              'DSMZ12361': 'black'}

# =============================================================================
# Locus tags
# =============================================================================
GS1 = ['A1001_12310', 'A1003_12540', 'A1202_13520', 'A1401_12750', 'A1805_12820',
       'FHON2_13540', 'G0101_12800', 'G0403_13100', 'H1B104J_13010', 
       'H1B105A_12300', 'H1B302M_12880', 'H3B101A_13240', 'H3B104J_12990', 
       'H3B104X_13200', 'H3B202X_12850', 'H3B203J_13370', 'H3B203M_12480', 
       'H3B206M_12830', 'H3B206M_12840', 'H3B209X_13340', 'H4B111J_13560', 'H4B111J_13570',
       'H4B202J_12880', 'H4B204J_13330', 'H4B205J_12990', 'H4B211M_13000',
       'H4B402J_12600', 'H4B406M_13450', 'H4B406M_13460', 'H4B412M_13240',
       'H4B501J_12890', 'H4B503X_12670', 'H4B504J_13460', 'H4B505J_12880',
       'APS55_RS03850', 'LDX55_06325', 'K2W83_RS06180']
GS2 = ['FHON2_13560', 'FHON2_13570', 'G0403_13120', 'H1B302M_12900', 
       'H3B101A_13260', 'H3B104X_13220', 'H3B202X_12860', 'H3B203J_13390',
       'H3B209X_13370', 'H4B504J_13480', 'H4B505J_12900', 'APS55_RS03845']
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
# Function definitions
# =============================================================================
def get_strain_name(locus_tag, get_strain = False):
    """Function to change locus tags that do not correspond with strain names
    to strain names. The function also removes RS from RefSeq locus tags to
    make the final tags that are plotted more readable. It also removes AKU
    from the non-RefSeq locus tags, as it indicates the species, which is
    the same for all the strains, hence redundant.
    
    Parameters
    ----------
    locus_tag : str
        The locus tag to be overwritten."""
        
    locus_tag = locus_tag.replace('LDX55', 'IBH001').replace('APS55', 'MP2').replace('K2W83', 'DSM').replace('RS', '').replace('FHON', 'Fhon').replace('AKU', '')
    if get_strain:
        locus_tag = locus_tag.split('_')[0]
        if locus_tag[0] == 'H':
            locus_tag = locus_tag[:4] + '-' + locus_tag[4:]
        elif 'DSM' in locus_tag:
            locus_tag = locus_tag.replace('DSM', 'DSMZ12361')
    return locus_tag

def get_tabs(phylo_dict, folder_path, outpath):
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for i in range(1, len(phylo_dict)):
        strain1_file = f'{folder_path}/{phylo_dict[i]}_genomic.fna'
        strain2_file = f'{folder_path}/{phylo_dict[i+1]}_genomic.fna'
        outfile = f'{outpath}/{GH_type}_{i}.tab'
        cline_input = cline_blast(cmd = 'blastn', query = strain2_file, subject = strain1_file, remote = False, out = outfile, outfmt = 7)
        os.system(str(cline_input));

def get_info(phylo_dict, tag_dict, folder_path2, length = 0): #OBS! Genome lengths are wrong because plasmid lengths are summed!
    locus_dict = {'IBH001': 'LDX55_06270', 'DSMZ12361': 'K2W83_RS06135', 'MP2': 'APS55_RS03895'}
    accession_strain = {}
    pos_dict = {}
    lengths_dict = {}
    for j in range(1, len(phylo_dict) + 1):
        strain = phylo_dict[j]
        loctag = tag_dict[j]
        genbank_file = f'{folder_path2}/{phylo_dict[j]}_genomic.gbff'
        with open(genbank_file) as handle:
            check = 0
            for record in gbk.parse(handle):
                if loctag not in accession_strain.keys():
                    accession_strain[loctag] = f'{record.accession[0]}.1'
                    genome_length = len(record.sequence)
                    lengths_dict[loctag] = [genome_length, check]
                else:
                    lengths_dict[loctag][0] = lengths_dict[loctag][0] + len(record.sequence)
                    check += len(record.sequence)
                    lengths_dict[loctag][1] = check
                if strain not in ['IBH001', 'MP2', 'DSMZ12361']:
                    for feature in record.features:
                        gene_name = [qual for qual in feature.qualifiers if 'gene' in qual.key]
                        if len(gene_name) > 0:
                            gene = gene_name[0]
                            if 'ohrR' in gene.value:
                                start = int(feature.location.split('..')[0]) + 2000

                else:
                    for feature in record.features:
                        locus_tag = [qual for qual in feature.qualifiers if 'locus_tag' in qual.key]
                        if len(locus_tag) > 0 and locus_tag[0].value.replace('"', '') == locus_dict[strain]:
                            start = int(feature.location.split('..')[0].replace('complement(', '')) + 5000
                            if strain == 'MP2':
                                start -= 10000
            if strain != 'MP2':
                start = lengths_dict[loctag][0] - (start + length)
                end = start + length
                pos_dict[loctag] = (start, end)
            else:
                end = start - length
                pos_dict[loctag] = (end, start)
    return accession_strain, pos_dict, lengths_dict

# =============================================================================
# Define inputs
# =============================================================================
outpath = os.path.expanduser('~') + '/GH_project/blast_tabs/subplots'
folder_path = os.path.expanduser('~') + '/Akunkeei_files/fna/reverse'
folder_path2 = os.path.expanduser('~') + '/Akunkeei_files/gbff/modified_gbff'
outfig = os.path.expanduser('~') + '/GH_project/plots/trees'
GH_types= ['GS1', 'GS2', 'BRS']

# =============================================================================
# In this section, we loop through GH types and run Blast between pairs of 
# strains following the phylogeny.
# =============================================================================
for GH_type in GH_types:
    treefile = os.path.expanduser('~') + f'/GH_project/data/fasta/GH70/trees/{GH_type}_repset.mafft.faa.treefile'
    t = Tree(treefile, format = 0)
    if GH_type == 'GS1':
        outnode = 'A1003_12540' # GS1
    elif GH_type == 'GS2':
        outnode = t.get_common_ancestor('H3B104X_13220', 'H4B505J_12900')
    elif GH_type == 'BRS':
        outnode = t.get_common_ancestor('H3B101A_13250', 'LDX55_06330')
    else: raise Exception(f'{GH_type} is not a valid GH70 gene type') 
    t.set_outgroup(outnode)
    generat_0r = list(t.get_leaf_names())
    phylo_order = {no + 1: generat_0r[-(no + 1)] for no in reversed(range(len(generat_0r)))}
    
    # OBS! Maybe I should re-do the Blast comparisons reversing the genome.
    
    phylo_list = [get_strain_name(leaf, True) for leaf in generat_0r]
    phylo_strains = {key: get_strain_name(phylo_order[key], True) for key in phylo_order.keys()}
    get_tabs(phylo_strains, folder_path, outpath)
    
    
    
    # =============================================================================
    # In this section, we establish the start and end positions for every strain
    # and save them to a dictionary. I also get the chromosome IDs to change them
    # to strain IDs.
    # =============================================================================
    acc, pos, lengths = get_info(phylo_strains, phylo_order, folder_path2, 40000)
    
    
    # =============================================================================
    # Now, ready for the plotting!
    # =============================================================================
    # Set plot style
    gv = GenomeViz(
        fig_track_height = 0.4,
        link_track_ratio = 3.0,
        align_type = 'center',
        tick_style = 'bar',
        )
    for i in range(1, len(phylo_strains.keys())+1):
        strain = phylo_strains[i]
        locus = phylo_order[i]
    
        genbank_file = f'{folder_path2}/{strain}_genomic.gbff'
        if strain != 'MP2':
            rev = True
        else: rev = False
        genbk = gbk_read(genbank_file, reverse = rev, min_range = pos[locus][0], max_range = pos[locus][1])
        features = genbk.extract_features('CDS')
        track_name = genbk.name.replace('_genomic', '').replace(strain, locus)
        locname = get_strain_name(locus)
            
        track = gv.add_feature_track(genbk.name.replace('_genomic', '').replace(strain, locname), size = genbk.range_size, start_pos = genbk.min_range, labelcolor = leaf_color[strain.replace('-', '')])
        #Loop through CDS
        for cds in features:
            protstart = int(cds.location.start) #Get CDS start
            end = int(cds.location.end) #Get CDS end
            strand = cds.strand #Get strand
            color = 'skyblue'
            if 'transposase' not in cds.qualifiers['product'][0] and 'gene' in cds.qualifiers.keys():
                gene_name = cds.qualifiers['gene'][0]
                if '_partial' in gene_name:
                    gene_name = gene_name.replace('_partial', '*')
                elif 'gene' in cds.qualifiers.keys() and '_I' in gene_name:
                    gene_name = gene_name.replace('_I', '')
                else:
                    gene_name = ''
                if gene_name == 'GS1' or gene_name[:-1] == 'GS1':
                    color = '#ff7575'
                elif gene_name == 'GS2' or gene_name[:-1] == 'GS2':
                    color = '#ff75b6'
                elif (gene_name == 'BRS' or gene_name[:-1] == 'BRS'):
                    color = '#e875ff'
                elif gene_name == 'GS2BRS':
                    color = '#ff75ee'
                elif gene_name  == 'S1':
                    color = '#ffb875'
                elif gene_name == 'S2a' or gene_name == 'S2b':
                    color = '#ffd775'
                elif gene_name == 'S3':
                    color = '#fff575'
            elif 'transposase' in cds.qualifiers['product'][0]:
                gene_name = ''
                color = 'black'
            else:
                gene_name = ''
            gene_n = gene_name.replace('*', '')
            if gene_n in ['branch', 'branch2', 'GTB', 'wzyC', 'GH39', 'rfbX', 'ydiL']:
                gene_name = ''
            elif gene_name == 'yifK2':
                gene_name = 'yifK'
            #This sets the CDS arrows and adds gene names
            track.add_feature(protstart, end, strand, label = gene_name, 
                              labelcolor = 'black', labelsize = 12, facecolor = color, 
                              linewidth = 1, labelrotation = 45, labelvpos = 'top', 
                              labelhpos = 'center', labelha = 'left', arrow_shaft_ratio = 1.0)
            track.set_sublabel(position = 'bottom-left')
            gene_name = ''
            
    # =============================================================================
    # Here I modify sligthly the tab files that I created.
    # =============================================================================
    # There is something wrong with the pairs in the tabs, they are not the same
    # as in the loop and that's why this doesn't work. Theory: the indexing of the
    # dictionary starts from 1, but normal indexing doesn't.
    for i in range(1, len(phylo_strains.keys())):
        tab_name = f'{outpath}/{GH_type}_{i}.tab'
        strain = phylo_strains[i]
        locus = phylo_order[i]
        next_locus = phylo_order[i + 1]
        print(locus, next_locus)
        with open(tab_name) as tabfile:
            k = 0
            for line in tabfile:
                    #print(line)
                    k += 1
                    #We can also use this loop to generate a list with the fields:
                    if 'Fields' in line:
                        headers = line
                    elif k > 10:
                        break
        header_list = headers.split(',') #Divide comma-separated fields
        header_list[0] = header_list[0].replace('# Fields: ', '') #Remove start of line
        header_list = [header.strip().strip('\n') for header in header_list] #Remove spaces and line breaks from beginning
        
        file_df = pd.read_csv(tab_name, sep = '\t', skiprows = 5, header = None) #Read as tab-separated file, skip 5 rows and don't set the remaining rows as headers
        file_df = file_df[:-1] #Remove the last row, that also contains a comment
        file_df.columns = header_list #Add column names
        strain2 = file_df.values[0, 0]
        strain2_name = list(acc.keys())[i] #list(acc.keys())[list(acc.values()).index(strain2)]
        file_df = file_df.replace(strain2, strain2_name)
        strain1 = file_df.values[0, 1]
        strain1_name = list(acc.keys())[i-1] #list(acc.keys())[list(acc.values()).index(strain1)]
        file_df = file_df.replace(strain1, strain1_name)
        # print(f'{i}:\t{strain1_name}\t{strain2_name}')
        query_pos =  (int(pos[phylo_order[i + 1]][0] - lengths[phylo_order[i + 1]][1]), int(pos[phylo_order[i + 1]][1] - lengths[phylo_order[i + 1]][1]))
        subject_pos = (int(pos[phylo_order[i]][0] - lengths[phylo_order[i]][1]), int(pos[phylo_order[i]][1] - lengths[phylo_order[i]][1]))
    
        tab_df = file_df.copy()
        # print(f'{i}:\t{tab_df.iloc[0]["query acc.ver"]}\t{tab_df.iloc[0]["subject acc.ver"]}')
        
        tab_df['query acc.ver'] =  tab_df['query acc.ver'].apply(lambda x: get_strain_name(str(x))) #tab_df['query acc.ver'].str.replace('LDX55', 'IBH001').replace('K2W83_RS', 'DSM_').replace('APS55_RS', 'MP2_').replace('FHON', 'Fhon')
        tab_df['subject acc.ver'] = tab_df['subject acc.ver'].apply(lambda x: get_strain_name(str(x))) #tab_df['subject acc.ver'].str.replace('LDX55', 'IBH001').replace('K2W83_RS', 'DSM_').replace('APS55_RS', 'MP2_').replace('FHON', 'Fhon')
        
        tab_df = tab_df.drop(tab_df[tab_df['q. end'] < query_pos[0]].index)
        tab_df = tab_df.drop(tab_df[tab_df['s. end'] < subject_pos[0]].index)
        tab_df['q. start'][tab_df['q. start'] < query_pos[0]] = query_pos[0]
        tab_df['s. start'][tab_df['s. start'] < subject_pos[0]] = subject_pos[0]
        tab_df = tab_df.drop(tab_df[tab_df['q. start'] > query_pos[1]].index)
        tab_df = tab_df.drop(tab_df[tab_df['s. start'] > subject_pos[1]].index)
        tab_df['q. end'][tab_df['q. end'] > query_pos[1]] = query_pos[1]
        tab_df['s. end'][tab_df['s. end'] > subject_pos[1]] = subject_pos[1]
        tab_df = tab_df[~tab_df['query acc.ver'].str.startswith('#')]
        tab_df = tab_df.reset_index()
        
        for j in range(len(tab_df)):
            #We define the links
            link1 = (tab_df.loc[j, 'query acc.ver'], tab_df.loc[j, 'q. start'] + lengths[next_locus][1], tab_df.loc[j, 'q. end'] + lengths[next_locus][1])
            link2 = (tab_df.loc[j, 'subject acc.ver'], tab_df.loc[j, 's. start'] + lengths[locus][1], tab_df.loc[j, 's. end'] + lengths[locus][1])
            #We define the percentage of identity (important to color the links)
            identity = tab_df.loc[j, '% identity']
            #Here we add the links to the gv plot
            gv.add_link(link1, link2, v = identity, vmin = 50, curve = True)
            gv.tick_style = 'bar' #This adds the scale of the plot (20 Kb)
            
    fig = gv.plotfig(400) #We plot the figure
    gv.set_colorbar(fig, vmin = 50, bar_height = 0.05) #We add a color bar to interpret the colors
    handles = [
        Line2D([], [], marker='', color='black', label='Tracks', ms=20, ls='none'),
        Line2D([], [], marker='>', color='skyblue', label='CDS', ms=20, ls='none'),
        Line2D([], [], marker='>', color='#ff7575', label='GS1 glycosyl hydrolase', ms=20, ls='none'),
        Line2D([], [], marker='>', color='#ff75b6', label='GS2 glycosyl hydrolase', ms=20, ls='none'),
        Line2D([], [], marker='>', color='#ff75ee', label='Double-GH70 domain', ms=20, ls='none'),
        Line2D([], [], marker='>', color='#e875ff', label='BRS glycosyl hydrolase', ms=20, ls='none'),
        Line2D([], [], marker='>', color='#ffb875', label='S1 glycosyl hydrolase', ms=20, ls='none'),
        Line2D([], [], marker='>', color='#ffd775', label='S2 glycosyl hydrolase', ms=20, ls='none'),
        Line2D([], [], marker='>', color='#fff575', label='S3 glycosyl hydrolase', ms=20, ls='none'),
        Line2D([], [], marker='>', color='black', label='Transposase', ms=20, ls='none'),
        Line2D([], [], marker='', color='#E5BA60', label='', ms=20, ls='none'),
        Line2D([], [], marker='', color='black', label='Matches', ms=20, ls='none'),
        Line2D([], [], marker='s', color='grey', label='Forward match', ms=20, ls='none'),
        Line2D([], [], marker='s', color='red', label='Reverse match', ms=20, ls='none'),
        Line2D([], [], marker='', color='#E5BA60', label='', ms=20, ls='none'),
        Line2D([], [], marker='', color='#E5BA60', label='Strain names', lw=2, ms=20, ls='none'),
        Line2D([], [], marker='X', color='#1E55F6', label='Phylogroup A', ms=15, ls='none'),
        Line2D([], [], marker='X', color='#00AEFF', label='Phylogroup B', ms=15, ls='none'),
        Line2D([], [], marker='X', color='#A027FF', label='Phylogroup C', ms=15, ls='none'),
        Line2D([], [], marker='X', color='#FF74D6', label='Phylogroup E', ms=15, ls='none'),
        Line2D([], [], marker='X', color='#E5BA60', label='Phylogroup F', ms=15, ls='none'),
        Line2D([], [], marker='X', color='black', label='Not in Dyrhage et al. (2022)', ms=15, ls='none')
    ]
    
    legend = fig.legend(handles=handles, bbox_to_anchor=(1.05, 1))
    legend.get_texts()[0].set_fontweight('bold')
    legend.get_texts()[11].set_fontweight('bold')
    legend.get_texts()[15].set_fontweight('bold')
    
    # Set the fontsize for legend labels
    legend_fontsize = 20
    for text in legend.get_texts():
        text.set_fontsize(legend_fontsize)
        
    fig.savefig(f'{GH_type}_blast.svg')
    fig.savefig(f'{GH_type}_blast.png')
