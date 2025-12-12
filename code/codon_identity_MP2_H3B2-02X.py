#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:07:19 2023

This script uses a codon-by-codon alignment to check if each codon position
is identical, synonymous or non-synonymous between genes, and then plots
each codon position, which results in a big matrix. The types of genes in
each pairwise comparison are also specified.

@author: Marina Mota-Merlo
"""
# =============================================================================
# 0. Import required packages
# =============================================================================

import os, subprocess
import logging, traceback, sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import matplotlib.font_manager as font_manager

# =============================================================================
# 0. Logging
# =============================================================================

log = os.path.expanduser('~') + '/GH_project/logs/dS_matrix.log' # Log file

if os.path.exists(log): # Remove previous log file if it exists
    os.remove(log)

if not os.path.exists(os.path.dirname(log)): # Create log directory if it does not exist
    os.makedirs(os.path.dirname(log))

# Create a logger
logger = logging.getLogger()

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

# Redirect errors and stout to log file
sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 0. Input and function definition
# =============================================================================

# Paths
aa_path = os.path.expanduser('~') + '/GH_project/data/fasta'
gene_types = ['GH70', 'GH32', 'other_genes']
neighbor = 'other_genes'
comparison = 'MP2_vs_H3B2-02X'
seq_outdir = f'{aa_path}/{neighbor}/subset/{comparison}'
aln_path = f'{seq_outdir}/alignments'
codon_path = os.path.expanduser('~') + '/GH_project/data/codons/other_genes/subset'
outdir = os.path.expanduser('~') + '/GH_project/plots/substitution_subsets'
pal2nal_path = os.path.expanduser('~') + '/GH_project/pal2nal.v14/pal2nal.pl'
outplot = f'{outdir}/{comparison}_matrix_subs_plot.png'

# Threads for mafft (set to automatic)
thr = -1

# Dictionaries with info about the data to be used
strains =  ['MP2', 'H3B2-02X']

gene_order = {29: 'ohrR', 28: 'yifK', 27: 'yifK2', 26: 'yhdG', 25: 'T5', 
              24: 'CDS8', 23: 'nox', 22: 'T4', 21: 'CDS7',
              20: 'GS1', 19: 'GS2', 18: 'T3', 17: 'mhpD', 
              16: 'oppA', 15: 'CDS6', 14: 'T2', 13: 'CDS5', 
              12: 'sbnD', 11: 'CDS4', 10: 'CDS3', 9: 'CDS2', 8: 'epsE', 
              7: 'CDS1', 6: 'epsL', 5: 'ywqE', 4: 'ywqD', 3: 'ywqC', 2: 'T1', 
              1: 'tagU'}

# Translation table for bacterial codons, assumes standard codon translations
transtable = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L',
              'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I',
              'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
              'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T',
              'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
              'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Stop',
              'TAG': 'Stop', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 
              'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C',
              'TGA': 'Stop', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
              'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

# List of glycosyl hydrolase genes, using the names that I gave them
GH_genes = ['GS1', 'GS2']

# Paths to be created if they don't exist
paths2make = [aln_path, codon_path, outdir]

#First: after tagU
#Second: after CDS5 or S3
#Third: after mhpD
#Fourth: After CDS7
#Last: after CDS8
transposons = {'H3B2-02X': [False, False, False, False, False],
               'MP2': [False, False, False, False, True]}

# Positions and size of transposon representations
y_list = [0.25, 0.75]
k = 200

# Set font size and style
fontsize_max = 84
fontsize_min = 77
font = font_manager.FontProperties(family='Arial', size = fontsize_min)

# =============================================================================
# 2. Function definitions
# =============================================================================

def make_paths(my_path):
    '''This function creates a path if it doesn't exist already
    Input: Path to be created'''
    if not os.path.exists(my_path):
        os.makedirs(my_path)
        print(f'Created directory {my_path}')
        
[make_paths(path) for path in paths2make]

def get_codons(sequence):  
    '''This function takes a sequence in a codon alignment and
    splits it into codons
    Input: Nucleotide codon sequence (with or without gaps)
    Output: List of codons or positions with gaps'''           
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    return codons

def retrieve_gene_seqs(input_dir):
    '''This function retrieves the GH genes and neighboring genes for the
    strains of interest (in strain_groups).
    Input: Path were the files are
    Output: Files with the desired sequences
    
    Note that the files have to follow a specific naming scheme:
        neighboring genes: other_genes/a_kunkeei_<gene>.<extension>, where the 
        extension can be faa or fna
        GH genes: <GH family>/<GH type>_repset.<extension>, where GH type has 
        to be one of the types in GH genes (or S1) and the extension has to be
        faa or fna'''
    for nbr in os.listdir(input_dir):
        check = False
        if 'other_genes' in input_dir and nbr.startswith('a_kunkeei') and 'mafft' not in nbr:
            gene_name = nbr.split('_')[2].split('.')[0]
            check = True
        elif ('GH' in input_dir and any(gh_gene in nbr for gh_gene in GH_genes) and 'all' in nbr and 'mafft' not in nbr) and (('GH70' in input_dir and 'complete' in nbr) or 'GH32' in input_dir):
            if 'GH70' not in input_dir:
                gene_name = nbr.split('_')[0]
            else:
                gene_name = nbr.split('_')[1]
            if 'BRS2' not in nbr and 'BRS3' not in nbr and 'BRS_clade' not in nbr:
                check = True
        if check:
            suffix = nbr.split('.')[-1]
            print(f'Suffix is {suffix}')
            print(f'Retrieving gene {gene_name}')
            with open(f'{input_dir}/{nbr}') as nbr_file:
                seqs = SeqIO.parse(nbr_file, 'fasta')
                for seq in seqs:
                    seq.id = seq.id.replace('APS', '')
                    seq.id = seq.id.replace('55_RS', 'MP2_')
                    print(seq.id)
                    seq.description = seq.description.replace('55_RS', 'MP2_')
                    seq.description = seq.description.replace('APS', '')
                    strain = seq.id.split('_')[0]
                    if strain[0] == 'H':
                        strain = strain[:4] + '-' + strain[4:]
                    if strain in strains:
                        with open(f'{seq_outdir}/{comparison}_{gene_name}.{suffix}', 'a+') as seqfile:
                            print(f'Writing gene {seq.id} to {seq_outdir}/{comparison}_{gene_name}.{suffix}')
                            SeqIO.write(seq, seqfile, 'fasta')
                            
def get_unique_values(tuple_list):
    '''This function loops through a tuple with pairs of locus tags and returns 
    a string with the names of the two strains in a pair'''
    strains = []
    for n1, n2 in tuple_list:
        strain = n1.split('_')[0]
        strain2 = n2.split('_')[0]
        strains.append(f'{strain}_{strain2}')
    return strains

def format_x(x, pos):
    '''This function specifies how to display the tick labels of the x axis'''
    return f'{x+1:.0f}'

# =============================================================================
# 0. Remove outputs of prior iterations
# =============================================================================

for file in os.listdir(seq_outdir):
    if file.endswith('.faa') or file.endswith('fna'):
        os.remove(f'{seq_outdir}/{file}')

# =============================================================================
# 1. From amino acid sequences of genes in the 30kb region, including glycosyl
# hydrolases, create faa and fna files
# =============================================================================

for folder in gene_types:
    retrieve_gene_seqs(f'{aa_path}/{folder}')
    
# =============================================================================
# 2. Align sequences in the faa files and create Pal2nal alignments
# =============================================================================

for file in os.listdir(seq_outdir):
    if file.endswith('faa'):
        outfile_base = file.replace('faa', 'mafft.faa')
        print(f'Get alignments for {file.replace(".faa", "")}')
        subprocess.run(f'mafft-linsi --thread {thr} {seq_outdir}/{file} > {aln_path}/{outfile_base} 2>> {log}',
                       shell = True)
        subprocess.run(f'{pal2nal_path} {aln_path}/{outfile_base} {seq_outdir}/{file.replace("faa", "fna")} -output fasta > {aln_path}/{outfile_base.replace("mafft.faa", "pal2nal.fna")} 2>> {log};',
                       shell = True)

# =============================================================================
# 3. Start the plotting!
# =============================================================================

# Get the length of each alignment to scale gene length in the plots
length_list = []
for gene_no in range(1, len(gene_order)+1):
    gene = gene_order[gene_no]
    file = f'{comparison}_{gene}.pal2nal.fna' # Input alignment
    
    if file in os.listdir(aln_path):
        print(f'Retrieving lengths from {aln_path}/{file}')
    
        with open(f'{aln_path}/{file}') as handle:
            f = SeqIO.parse(handle, 'fasta')
            for record in f:
                aln_length = len(record.seq)//3
                length_list.append(aln_length)
                break
    else: length_list.append(0)
            
# Assign the lengths to transposons
for t_no in range(1, 6):
    gene = f'T{t_no}'
    gene_number = [key for key in gene_order if gene_order[key] == gene][0]
    if length_list[gene_number-1] == 0:
        length_list[gene_number-1] = k
        

# Create a dataframe assigning different values to identical codons (0),
# synonymous substitutions (1), non-synonymous substitutions (2) and gaps (3)
df_aln_dict = {}
fstrains = []

# Create the figure to save the plots
fig = plt.figure(constrained_layout=False, figsize=(160, 10))
spec = gridspec.GridSpec(nrows=1, ncols=29, figure=fig, width_ratios=length_list)
count = 0 # Number of the row

plt.tight_layout(h_pad = 2, w_pad = 2) # Adjust space between genes
for gene_no in range(1, len(gene_order)+1):
    gene = gene_order[gene_no]
    
    if gene.startswith('T'):
            
        ax = fig.add_subplot(spec[count, gene_no-1]) # Set the right row and column
        ax.set_ylim(0, 1)
        ax.set_xlim(0, k)
        T_no = int(gene[-1])
            
        column = df_aln.columns[0]
        if type(column) == tuple:
            strain1 = column[0].split('_')[0]
            strain2 = column[1].split('_')[0]
        else:
            strain1 = column.split('_')[0]
            strain2 = column.split('_')[1]
        if len(strain1) > 6 and '-' not in strain1 and strain1.startswith('H'):
            strain1 = strain1[:4] + '-' + strain1[4:]
        if len(strain2) > 6 and '-' not in strain2 and strain2.startswith('H'):
            strain2 = strain2[:4] + '-' + strain2[4:]
        T_presence1 = transposons[strain1][T_no-1] 
        T_presence2 = transposons[strain2][T_no-1] 
        
        print(f'Strain names: {strain1}/{strain2}')
        
        
        # Add arrows indicating gene direction
        if count == 0: # Add this only on top of the first row
            print(f'Comparison is ... {comparison}')
            ax.set_title(gene, fontsize = fontsize_min, style = 'italic', 
                         fontname = 'Arial') # Add gene names at the top of the plot
            color = 'white' # Add background color to arrows
                
            if 4 > T_no >= 2: # Invert arrows for genes in the forward strand
                startx = k
                dx = -k
            else: # Otherwise, plot the arrows facing toward the right
                startx = 0
                dx = k
            ax.arrow(startx, 1.304, dx, 0, width = 0.068, head_width = 0.068, 
                     head_length = 30, length_includes_head = True, 
                     clip_on = False, linewidth = 1, linestyle = '--',
                     edgecolor = 'black', facecolor = color)
            
        # Colors of transposon representations
        if T_presence1 == 1:
            circle_col = 'black'
        else: circle_col = '#F8F9F9'
        if T_presence2 == 1:
            circle2_col = 'black'
        else: circle2_col = '#F8F9F9'
        
        # Plot transposon presence/absence
        ax.scatter([k//2], y_list[0], s = 6000, edgecolor = 'black',  
                   facecolor = circle_col)
        ax.scatter([k//2], y_list[1], s = 6000, edgecolor = 'black',  
                   facecolor = circle2_col)
        ax.set_axis_off()
            
        print(f'Added {gene} presence/absence plot to {comparison}')
        
    else:
        
        file = f'{comparison}_{gene}.pal2nal.fna' # Input alignment
        
        if file in os.listdir(aln_path): # If the gene is present in at least two strains of a comparison
            print(f'Retrieving codons from {aln_path}/{file}')
        
            pos_dict = {}
            with open(f'{aln_path}/{file}') as handle:
                f = SeqIO.parse(handle, 'fasta')
                for record in f:
                    pos_dict[record.id] = str(record.seq) # Get every sequence in the alignment and save it to a dictionary
            
            pair_dict = {}
            locus_tags = list(pos_dict.keys())
            aln_length = len(pos_dict[record.id])//3
            for i in range(len(locus_tags) - 1):
                for j in range(i + 1, len(locus_tags)): # Loop through every pair of sequences
                    pair_dict[(locus_tags[i], locus_tags[j])] = [] # Create empty list for the pair
                    for n in range(0, len(pos_dict[locus_tags[i]]), 3): # Loop through the codons in the alignment
                        if '-' in pos_dict[locus_tags[i]][n:n+3] or '-' in pos_dict[locus_tags[j]][n:n+3]: # Append different numbers to the empty list based on the codon comparison
                            pair_dict[(locus_tags[i], locus_tags[j])].append(3) # 3 for gaps
                        elif pos_dict[locus_tags[i]][n:n+3] == pos_dict[locus_tags[j]][n:n+3]:
                            pair_dict[(locus_tags[i], locus_tags[j])].append(0) # 0 for identical codons
                        elif transtable[pos_dict[locus_tags[i]][n:n+3]] == transtable[pos_dict[locus_tags[j]][n:n+3]]:
                            pair_dict[(locus_tags[i], locus_tags[j])].append(1) # 1 for synonymous codons
                        else: pair_dict[(locus_tags[i], locus_tags[j])].append(2) # 2 for non-synonymous codons
                        
            df_aln = pd.DataFrame.from_dict(pair_dict) # Convert dictionary to dataframe
            
            fstrains = get_unique_values(list(df_aln.columns))


        df_aln.index += 1 # Make the index start with 1, not 0
        df_aln_dict[gene] = df_aln # Save dataframe to dictionary
        # Assign colors to the values in the dataframe
        cmap = []
        if 0 in df_aln_dict[gene].values:
            cmap.append('#FEFDED')
        if 1 in df_aln_dict[gene].values:
            cmap.append('#80C4E9')
        if 2 in df_aln_dict[gene].values:
            cmap.append('#FF7F3E')
        if 3 in df_aln_dict[gene].values:
            cmap.append('#F8F9F9')
        
        # Create a plot for the gene and comparison
        ax = fig.add_subplot(spec[count, gene_no-1]) # Set the right row and column
        h3 = sns.heatmap(df_aln_dict[gene].T, cmap = cmap, cbar = False, ax = ax) # Plot
        ax.patch.set(lw = 1, ec = 'black') # Border of the subplots
        gn = gene # Change some of the gene annotations
        
        # Add arrows indicating gene direction
        ax.set_title(gn, fontsize = fontsize_min, style = 'italic',
                     fontname = 'Arial') # Add gene names at the top of the plot
        color = 'white' # Add background color to arrows
            
        startx = 0
        dx = max(df_aln.index)
        
        ax.arrow(startx, -2.165*0.14, dx, 0, width = 0.068, head_width = 0.068, 
                 head_length = 30, length_includes_head = True, 
                 clip_on = False, linewidth = 1,
                 edgecolor = 'black', facecolor = color)
            
        # Add ticks to the y axis if the gene is present in the comparison
        if count == 0:
            plt.yticks(rotation=90)
            
        plt.tick_params(axis='x', which='major', labelsize = fontsize_min,
                        length = 20) # Increase font size of ax ticks
        
        if gene_no >= 27 or gene_no == 15: #For the genes that are in the forward strand
            ax.invert_xaxis() #Invert the axis
        ax.xaxis.set_major_locator(MultipleLocator(150))
        ax.xaxis.set_major_formatter(format_x)
        for label in ax.get_xticklabels():
            label.set_fontproperties(font)
        ax.set(xlabel=None) # Remove x ax label (label on the x axis of the subplot)
        ax.set(ylabel=None) # Remove y ax label

        # If the gene is not the first one, then remove the y ax ticks and labels
        if gene != 'tagU':
            plt.tick_params(left = False, labelleft = False)
            
        # Otherwise, plot the comparisons as left labels (ytick labels)
        else:
            y_ticks = get_unique_values(list(df_aln.columns))
            y_ticks = [ytick.replace('_', ' vs ') for ytick in y_ticks]
            ax.set_yticklabels(y_ticks, fontsize = fontsize_min, 
                               rotation = 'vertical', family = 'Arial')
            ax.set_ylabel('Comparison', fontsize = fontsize_max, family = 'Arial')
        
        print(f'Added {gene} plot to comparison {comparison}!')
    
count += 1
print(f'Plot for {comparison} saved to {outplot}!')
    
# Add global legend
gap_patch = mpatches.Patch(color = '#F8F9F9', label = 'Gaps', ec = 'black', 
                           lw = 0.5)

identity_patch = mpatches.Patch(color = '#FEFDED', label = 'Identical sites', 
                                ec = 'black', lw = 0.5)

S_patch = mpatches.Patch(color = '#80C4E9', label = 'Synonymous sites', 
                         ec = 'black', lw = 0.5)

NS_patch = mpatches.Patch(color = '#FF7F3E', label = 'Non-synonymous sites', 
                          ec = 'black', lw = 0.5)

t2_patch = Line2D([], [], markeredgecolor = 'black', marker = 'o', markerfacecolor = 'black',
                  color = 'white', ms = 75, label = 'Transposon in strain')

t0_patch = Line2D([], [], markeredgecolor = 'black', marker = 'o', markerfacecolor = '#F8F9F9', 
           color = 'white', ms = 75, label = 'Transposon not in strain')

patches = [gap_patch, identity_patch, S_patch, NS_patch, t2_patch, t0_patch]


legend = plt.legend(handles = patches,  handlelength = 0.35, handleheight = 1.2, 
                    title = 'Site color', loc = 'upper right', framealpha = 0, 
                    frameon = False, fontsize = fontsize_min, ncol = 6,
                    bbox_to_anchor = (1, 2), prop = font)

plt.setp(legend.get_title(),fontsize = fontsize_max, family = 'Arial')
legend._legend_box.align = 'left'

# Save plot
fig.text(0.5, -0.2, 'Codon position', ha='center', fontsize = fontsize_max,
         family = 'Arial')
plt.savefig(outplot, bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.pdf'), bbox_inches='tight')
