#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:28:42 2023

@author: marina
"""
from Bio.Blast.Applications import NcbiblastpCommandline as cline_blast
from pygenomeviz import Genbank as gbk_read, GenomeViz
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from Bio import GenBank as gbk
import os
import pandas as pd

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
NGB = ['A0901_05380', 'A1003_04750', 'A1202_05530', 'A1401_04720', 'A1805_04920',
       'FHON2_04830', 'G0101_04800', 'H1B105A_04750', 'H1B302M_04960', 
       'H3B101A_04720', 'H3B104X_04750', 'H3B202X_04770', 'H3B203J_04720', 
       'H3B203M_04710', 'H3B206M_05060', 'H3B209X_04660', 'H4B111J_04920', 
       'H4B202J_04590', 'H4B204J_04710', 'H4B205J_04990', 'H4B211M_04810', 
       'H4B405J_04950', 'H4B406M_05270', 'H4B503X_04680', 'H4B504J_05510', 
       'H4B505J_04980', 'K2W83_RS02365']

phage_pos = {'MP2': (190, 500), 
             'A1202': (4780, 5300+160), 
             'H4B5-04J': (4840, 5330+130), 
             'A0901': (4670, 5190+160)}

# =============================================================================
# In this section, we run Blast between pairs of strains following the 
# phylogeny.
# =============================================================================

phylo_order = {1:'G0101', 2:'H4B2-05J', 3:'H4B2-11M', 4:'H1B1-04J', 
              5:'H4B4-12M', 6:'H4B5-01J', 7:'H3B1-01A', 8:'H3B2-03J', 
              9:'H3B1-04X', 10:'H3B2-02X', 11:'A1003', 12:'H4B2-04J', 
              13:'H4B5-03X', 14:'H4B4-02J', 15:'H1B3-02M', 16:'Fhon2', 
              17:'DSMZ12361', 18:'H4B2-02J', 19:'H1B1-05A', 20:'H3B1-04J', 
              21:'A1401', 22:'A1202', 23:'H4B5-05J', 24:'H3B2-09X', 
              25:'H4B5-04J', 26:'A1805', 27:'MP2', 28:'G0403', 29:'H3B2-06M', 
              30:'H4B4-05J', 31:'H4B1-11J', 32:'H4B4-06M', 33:'A0901', 
              34:'IBH001', 35:'H4B2-06J', 36:'H3B2-03M', 37:'A1404', 38:'A1001'}


outpath = 'blast_tabs'
folder_path = '~/Akunkeei_files/fna'
folder_path2 = os.path.expanduser('~') + '/Akunkeei_files/gbff'

# OBS! Maybe I should re-do the Blast comparisons reversing the genome.

def get_tabs(phylo_dict, folder_path, outpath):
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for i in range(1, len(phylo_dict)):
        strain1_file = f'{folder_path}/{phylo_dict[i]}_genomic.fna'
        strain2_file = f'{folder_path}/{phylo_dict[i+1]}_genomic.fna'
        outfile = f'{outpath}/NGB_{i}.tab'
        cline_input = cline_blast(cmd = 'blastn', query = strain2_file, subject = strain1_file, remote = False, out = outfile, outfmt = 7)
        os.system(str(cline_input));

get_tabs(phylo_order, folder_path, outpath)



# =============================================================================
# In this section, we establish the start and end positions for every strain
# and save them to a dictionary. I also get the chromosome IDs to change them
# to strain IDs.
# =============================================================================

def get_info(phylo_dict, folder_path2, length = 0): #OBS! Genome lengths are wrong because plasmid lengths are summed!
    locus_dict = {'IBH001': 'LDX55_02300', 'DSMZ12361': 'K2W83_RS02275', 'MP2': 'APS55_RS00555'}
    accession_strain = {}
    lengths_dict = {}
    pos_dict = {}
    for strain in phylo_dict.values():
        genbank_file = f'{folder_path2}/{strain}_genomic.gbff'
        start = 0
        with open(genbank_file) as handle:
            # check = 0
            for record in gbk.parse(handle):
                if strain not in accession_strain.keys():
                    accession_strain[strain] = f'{record.accession[0]}.1'
                    genome_length = len(record.sequence)
                    lengths_dict[strain] = genome_length #[genome_length, check]
                # else:
                #     lengths_dict[strain][0] = lengths_dict[strain][0] + len(record.sequence)
                #     check += len(record.sequence)
                #     lengths_dict[strain][1] = check
                if strain not in ['IBH001', 'MP2', 'DSMZ12361']:
                    for feature in record.features:
                        gene_name = [qual for qual in feature.qualifiers if 'gene' in qual.key]
                        if len(gene_name) > 0:
                            gene = gene_name[0]
                            if 'pspB' in gene.value:
                                start = int(feature.location.split('..')[0]) + 2000
                                if start < 4*10^5:
                                    start = 0
                                else:
                                    break
                            elif start == 0 and 'pspA' in gene.value:
                                start = int(feature.location.split('..')[0]) + 2000
                                if strain == 'H1B3-02M':
                                    start += 10000
                                break



                else:
                    for feature in record.features:
                        locus_tag = [qual for qual in feature.qualifiers if 'locus_tag' in qual.key]
                        if len(locus_tag) > 0 and locus_tag[0].value.replace('"', '') == locus_dict[strain]:
                            start = int(feature.location.split('..')[0].replace('complement(', '')) + 2000
                            # if strain == 'MP2':
                            #     start -= 10000

            if strain != 'MP2':
                end = start + length
                pos_dict[strain] = (start, end)
            else:
                start = lengths_dict[strain] - start
                end = start + length
                pos_dict[strain] = (start, end)
    return accession_strain, pos_dict, lengths_dict
    
acc, pos, lengths = get_info(phylo_order, folder_path2, 40000)
                           
# OBS! There are some iisues, probably because I changed the order of the genomes. The Blast matches don't fully align with the CDS plots.


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
for i in range(1, len(phylo_order.keys())+1):
    strain = phylo_order[i]
    genbank_file = f'{folder_path2}/{strain}_genomic.gbff'
    if strain != 'MP2':
        rev = False
    else: rev = True
    if strain in ['A1202', 'MP2', 'H4B5-04J', 'A0901']:
        genbk2 = gbk_read(genbank_file, reverse = rev, min_range = pos[strain][0], max_range = pos[strain][1]+40000)
    genbk = gbk_read(genbank_file, reverse = rev, min_range = pos[strain][0], max_range = pos[strain][1])
    features = genbk.extract_features('CDS')
    track = gv.add_feature_track(genbk.name.replace('_genomic', ''), size = genbk.range_size, start_pos = genbk.min_range, labelcolor = leaf_color[strain.replace('-', '')])
    if strain in ['A1202', 'MP2', 'H4B5-04J', 'A0901']:
        features = genbk2.extract_features('CDS')
    #Loop through CDS
    for cds in features:
        protstart = int(cds.location.start) #Get CDS start
        end = int(cds.location.end) #Get CDS end
        strand = cds.strand #Get strand
        color = 'skyblue'
        gene_name  = ''
        if 'transposase' not in cds.qualifiers['product'][0]:
            if cds.qualifiers['locus_tag'][0][3:] in NGB or cds.qualifiers['locus_tag'][0] in NGB:
                gene_name = 'NGB'
                color = '#c475ff' # Add something similar for cas1, cas9 and other genes of interest
            elif 'restriction' in cds.qualifiers['product'][0].lower() or 'endonuclease' in cds.qualifiers['product'][0].lower() or 'HsdR' in cds.qualifiers['product'][0]:
                color = '#ffa375'
                if 'gene' in cds.qualifiers.keys():
                    gene_name = cds.qualifiers['gene'][0]
            elif 'ABC' in cds.qualifiers['product'][0]:
                color = '#b8ff75'
            elif strain in phage_pos.keys():
                loc_num = int(cds.qualifiers['locus_tag'][0].split('_')[1].replace('RS', ''))
                if (loc_num >= phage_pos[strain][0] and loc_num <= phage_pos[strain][1]) or (cds.qualifiers['locus_tag'][0] == 'APS55_RS08015'):
            #elif ('phage' in cds.qualifiers['product'][0]) or ('tape measure' in cds.qualifiers['product'][0]) or ('capsid' in cds.qualifiers['product'][0]) or ('terminase' in cds.qualifiers['product'][0]) or ('integrase' in cds.qualifiers['product'][0]):
                    color = '#7581ff'
            elif 'gene' in cds.qualifiers.keys():
                gene_name = cds.qualifiers['gene'][0]
                if 'cas2' in gene_name:
                    color = '#ffa375'
        elif 'transposase' in cds.qualifiers['product'][0]:
            gene_name = ''
            color = 'black'
            
        # track.add_feature(protstart, end, strand, label = gene_name, 
        #                       labelcolor = 'black', labelsize = 12, facecolor = color, 
        #                       linewidth = 1, labelrotation = 45, labelvpos = 'top', 
        #                       labelhpos = 'center', labelha = 'left', arrow_shaft_ratio = 1.0)
        # track.set_sublabel(position = 'bottom-left')
         
        
        if strain in ['A1202', 'MP2', 'H4B5-04J', 'A0901'] and (pos[strain][0]+16000 < end < pos[strain][1]+12000 or end > pos[strain][1]+36000):
            print('Exclude')
        elif strain in ['A1202', 'MP2', 'H4B5-04J', 'A0901'] and pos[strain][1]+12000 <= end < pos[strain][1]+36000:
            protstart -= 36000
            end -= 36000
            track.add_feature(protstart, end, strand, label = gene_name, 
                              labelcolor = 'black', labelsize = 12, facecolor = color, 
                              linewidth = 1, labelrotation = 45, labelvpos = 'top', 
                              labelhpos = 'center', labelha = 'left', arrow_shaft_ratio = 1.0)
            track.set_sublabel(position = 'bottom-left') 
        else:
            track.add_feature(protstart, end, strand, label = gene_name, 
                              labelcolor = 'black', labelsize = 12, facecolor = color, 
                              linewidth = 1, labelrotation = 45, labelvpos = 'top', 
                              labelhpos = 'center', labelha = 'left', arrow_shaft_ratio = 1.0)
            track.set_sublabel(position = 'bottom-left')

            
# =============================================================================
# Here I modify sligthly the tab files that I created.
# =============================================================================
for i in range(1, len(phylo_order.keys())):     
    with open(f'{outpath}/NGB_{i}.tab') as tabfile:
        k = 0
        for line in tabfile:
                #print(line)
                k += 1
                #We can also use this loop to generate a list with the fields:
                if "Fields" in line:
                    headers = line
                elif k > 10:
                    break
    header_list = headers.split(",") #Divide comma-separated fields
    header_list[0] = header_list[0].replace("# Fields: ", "") #Remove start of line
    header_list = [header.strip().strip("\n") for header in header_list] #Remove spaces and line breaks from beginning
    
    file_df = pd.read_csv(f'{outpath}/NGB_{i}.tab', sep = "\t", skiprows = 5, header = None) #Read as tab-separated file, skip 5 rows and don't set the remaining rows as headers
    file_df = file_df[:-1] #Remove the last row, that also contains a comment
    file_df.columns = header_list #Add column names
    strain2 = file_df.values[0, 0]
    strain2_name = list(acc.keys())[list(acc.values()).index(strain2)]
    file_df = file_df.replace(strain2, strain2_name)
    strain1 = file_df.values[0, 1]
    strain1_name = list(acc.keys())[list(acc.values()).index(strain1)]
    file_df = file_df.replace(strain1, strain1_name)
    print(strain1_name, strain2_name)
    
    query_pos = pos[phylo_order[i + 1]]
    subject_pos = pos[phylo_order[i]]
    tab_df = file_df.copy()
    if phylo_order[i] in ['A1202', 'MP2', 'H4B5-04J', 'A0901']:
        subject_pos = (subject_pos[0], subject_pos[1] + 36000)
    elif phylo_order[i+1] in ['A1202', 'MP2', 'H4B5-04J', 'A0901']:
        query_pos = (query_pos[0], query_pos[1] + 36000)
    
    tab_df = tab_df.drop(tab_df[tab_df["q. end"] < query_pos[0]].index)
    tab_df = tab_df.drop(tab_df[tab_df["q. start"] > query_pos[1]].index)
    tab_df["q. start"][tab_df["q. start"] < query_pos[0]] = query_pos[0]
    tab_df["q. end"][tab_df["q. end"] > query_pos[1]] = query_pos[1]

    tab_df = tab_df.drop(tab_df[tab_df["s. end"] < subject_pos[0]].index)
    tab_df = tab_df.drop(tab_df[tab_df["s. start"] > subject_pos[1]].index)
    tab_df["s. end"][tab_df["s. end"] > subject_pos[1]] = subject_pos[1]
    tab_df["s. start"][tab_df["s. start"] < subject_pos[0]] = subject_pos[0]
    
    if phylo_order[i+1] in ['A1202', 'MP2', 'H4B5-04J', 'A0901']:
        tab_df = tab_df.drop(tab_df[(tab_df["q. start"] <= query_pos[1]-24000) & (tab_df["q. start"] >= query_pos[0]+16000)].index)
        tab_df = tab_df.drop(tab_df[(tab_df["q. end"] <= query_pos[1]-24000) & (tab_df["q. end"] >= query_pos[0]+16000)].index)
        tab_df.loc[tab_df["q. start"] > query_pos[1]-24000, "q. end"] -= 36000
        tab_df.loc[tab_df["q. start"] > query_pos[1]-24000, "q. start"] -= 36000
        # Line to change links  over +20000
    elif phylo_order[i] in ['A1202', 'MP2', 'H4B5-04J', 'A0901']:
        tab_df = tab_df.drop(tab_df[(tab_df["s. start"] <= subject_pos[1]-24000) & (tab_df["s. start"] >= subject_pos[0]+16000)].index)
        tab_df = tab_df.drop(tab_df[(tab_df["s. end"] <= subject_pos[1]-24000) & (tab_df["s. end"] >= subject_pos[0]+16000)].index)
        tab_df.loc[tab_df["s. start"] > subject_pos[1]-24000, "s. end"] -= 36000
        tab_df.loc[tab_df["s. start"] > subject_pos[1]-24000, "s. start"] -= 36000
        
    
    tab_df = tab_df[~tab_df["query acc.ver"].str.startswith('#')]
    tab_df = tab_df.reset_index()
    
    for j in range(len(tab_df)):
        #We define the links
        link1 = (tab_df.loc[j, "query acc.ver"], tab_df.loc[j, "q. start"], tab_df.loc[j, "q. end"])
        link2 = (tab_df.loc[j, "subject acc.ver"], tab_df.loc[j, "s. start"], tab_df.loc[j, "s. end"])
        #We define the percentage of identity (important to color the links)
        identity = tab_df.loc[j, "% identity"]
        #Here we add the links to the gv plot
        gv.add_link(link1, link2, v = identity, vmin = 50, curve = True)
        gv.tick_style = "bar" #This adds the scale of the plot (20 Kb)
        
fig = gv.plotfig(400) #We plot the figure
gv.set_colorbar(fig, vmin = 50, bar_height = 0.05) #We add a color bar to interpret the colors

handles = [
    Line2D([], [], marker="", color='black', label="Tracks", ms=20, ls="none"),
    Line2D([], [], marker=">", color="skyblue", label="CDS", ms=20, ls="none"),
    Line2D([], [], marker=">", color='#7581ff', label="Phage genes", ms=20, ls="none"),
    Line2D([], [], marker=">", color="#ffa375", label="CRISPR endonucleases", ms=20, ls="none"),
    Line2D([], [], marker=">", color='#c475ff', label="NGB glycosyl hydrolase", ms=20, ls="none"),
    Line2D([], [], marker=">", color='#b8ff75', label="ABC transporters", ms=20, ls="none"),
    Line2D([], [], marker=">", color='black', label="Transposase", ms=20, ls="none"),
    Line2D([], [], marker="", color='black', label="", ms=20, ls="none"),
    Line2D([], [], marker="", color='black', label="Matches", ms=20, ls="none"),
    Line2D([], [], marker="s", color='grey', label="Forward match", ms=20, ls="none"),
    Line2D([], [], marker="s", color='red', label="Reverse match", ms=20, ls="none"),
    Line2D([], [], marker="", color='#E5BA60', label="", ms=20, ls="none"),
    Line2D([], [], marker="", color='#E5BA60', label="Strain names", lw=2, ms=20, ls="none"),
    Line2D([], [], marker="X", color='#1E55F6', label="Phylogroup A", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#00AEFF', label="Phylogroup B", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#A027FF', label="Phylogroup C", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#FF74D6', label="Phylogroup E", ms=15, ls="none"),
    Line2D([], [], marker="X", color='#E5BA60', label="Phylogroup F", ms=15, ls="none"),
    Line2D([], [], marker="X", color='black', label="Not in Dyrhage et al. (2022)", ms=15, ls="none")
]
legend = fig.legend(handles=handles, bbox_to_anchor=(1.05, 1))
legend.get_texts()[0].set_fontweight('bold')
legend.get_texts()[8].set_fontweight('bold')
legend.get_texts()[12].set_fontweight('bold')
# Set the fontsize for legend labels
legend_fontsize = 20
for text in legend.get_texts():
    text.set_fontsize(legend_fontsize)
            
fig.savefig('figure5.svg')

gv.savefig_html('figure5.html')
    
