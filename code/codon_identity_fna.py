#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:07:19 2023

This script uses a codon-by-codon alignment to check if each codon position
is identical, synonymous or non-synonymous between genes, and then plots
each codon position, which results in a big matrix. The types of genes in
each pairwise comparison are also specified.

@author: Marina Mota Merlo
"""
# Make the y axis ticks appear only on the ohrR plot, with higher fontsize
import os, subprocess
import logging, traceback, sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.gridspec as gridspec
# =============================================================================
# 0. Logging
# =============================================================================
log = os.path.expanduser('~') + '/GH_project/logs/dS_matrix.log' # Log file

if os.path.exists(log): # Remove årevious log file if it exists
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
seq_outdir = f'{aa_path}/{neighbor}/subset'
aln_path = f'{seq_outdir}/alignments'
codon_path = os.path.expanduser('~') + '/GH_project/data/codons/other_genes/subset'
outdir = os.path.expanduser('~') + '/GH_project/plots/substitution_subsets'
pal2nal_path = os.path.expanduser('~') + '/GH_project/pal2nal.v14/pal2nal.pl'
outplot = f'{outdir}/matrix_subs_plot.png'

# Threads for mafft (set to automatic)
thr = -1

# Dictionaries with info about the data to be used
strain_groups = {'H3B2-03M': 0, 'H4B4-02J': 0, 'H4B5-03X': 0, 'H1B3-02M': 1,
                 'H3B2-09X': 1, 'H4B5-05J': 1, 'MP2': 2, 'H3B2-02X': 2,
                 'H1B1-05A': 2}

group_names = {0: 'GS1_S2-3_subset', 1: 'GS1-2_BRS', 2: 'only_GS1+GS2'}

gene_order = {27: 'ohrR', 26: 'yifK', 25: 'yifK2', 24: 'yhdG', 23: 'CDS8', 
              22: 'nox', 21: 'CDS7', 20: 'S2a', 19: 'GS1', 18: 'BRS', 
              17: 'GS2', 16: 'mhpD', 15: 'oppA', 14: 'S3', 13: 'CDS6', 12: 'CDS5', 
              11: 'sbnD', 10: 'CDS4', 9: 'CDS3', 8: 'CDS2', 7: 'epsE', 
              6: 'CDS1', 5: 'epsL', 4: 'ywqE', 3: 'ywqD', 2: 'ywqC', 1: 'tagU'}

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
GH_genes = ['GS1', 'GS2', 'BRS', 'S2a', 'S3']

# Paths to be created if they don't exist
paths2make = [aln_path, codon_path, outdir]

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
        # check2 = False
        check = False
        add_gene = True
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
                    if strain in strain_groups.keys():
                        prefix = group_names[strain_groups[strain]]
                        with open(f'{seq_outdir}/{prefix}_{gene_name}.{suffix}', 'a+') as seqfile:
                            print(f'Writing gene {seq.id} to {seq_outdir}/{prefix}_{gene_name}.{suffix}')
                            # Here (commented section) I was trying to add a gene that was split in two, but wasn't sucessful:
                                
                            # if strain != 'H1B3-02M' and prefix == group_names[1] and gene_name == GH_genes[1] and check2 == False:
                            #     with open(f'{aa_path}/{gene_types[0]}/trunc_{GH_genes[1]}.{suffix}') as gs2_trunc:
                            #         record = SeqIO.read(gs2_trunc, 'fasta')
                            #         SeqIO.write(record, seqfile, 'fasta')
                            #     check2 = True
                            if seq.id == 'H3B209X_13360':
                            #     seq_add = str(seq.seq)
                                add_gene = False
                            # elif seq.id == 'H3B209X_13370':
                            #     seq.id = seq.id.replace('70', '60-70')
                            #     seq.seq = Seq(str(seq.seq) + seq_add)
                            #     seq.description = seq.description.replace('70', '60-70')
                            #     add_gene = True
                            else: add_gene = True
                            if add_gene:
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
length_dict = {}
for comparison in group_names.values():
    length_dict[comparison] = []
    for gene_no in range(1, len(gene_order)+1):
        gene = gene_order[gene_no]
        file = f'{comparison}_{gene}.pal2nal.fna' # Input alignment
        
        if file in os.listdir(aln_path) and not (('only' in file and gene == 'BRS') or ('BRS2' in file or 'BRS3' in file or 'clade' in file)):
            print(f'Retrieving lengths from {aln_path}/{file}')
        
            with open(f'{aln_path}/{file}') as handle:
                f = SeqIO.parse(handle, 'fasta')
                for record in f:
                    aln_length = len(record.seq)//3
                    length_dict[comparison].append(aln_length)
                    break
        else: length_dict[comparison].append(0)
        
# Get the lengths of the genes that are absent in some of the comparisons
# To do this, check other comparisons
absent_genes = {}
for comparison in group_names.values():
    for gene in ['S2a', 'GS2', 'BRS', 'S3']:
        gene_number = list(filter(lambda x: gene_order[x] == gene, gene_order))[0]
        if length_dict[comparison][gene_number-1] != 0:
            absent_genes[gene] = (gene_number, length_dict[comparison][gene_number-1])
            
# Assign the lengths retrieved from other comparisons
for comparison in group_names.values():
    for gene in ['S2a', 'GS2', 'BRS', 'S3']:
        gene_number = absent_genes[gene][0]
        if length_dict[comparison][gene_number-1] == 0:
            length_dict[comparison][gene_number-1] = absent_genes[gene][1]
        

# Create a dataframe assigning different values to identical codons (0),
# synonymous substitutions (1), non-synonymous substitutions (2) and gaps (3)
df_aln_dict = {}
# pos_dict_2 = {}
fstrains = []

# Create the figure to save the plots
fig = plt.figure(constrained_layout=False, figsize=(160, 30))
spec = gridspec.GridSpec(nrows=3, ncols=27, figure=fig, width_ratios=length_dict[comparison])
count = 0 # Number of the row
for comparison in group_names.values():
    plt.tight_layout(h_pad = 2, w_pad = 2) # Adjust space between genes
    for gene_no in range(1, len(gene_order)+1):
        gene = gene_order[gene_no]
            
        file = f'{comparison}_{gene}.pal2nal.fna' # Input alignment
        
        if file in os.listdir(aln_path) and not (('only' in file and gene == 'BRS') or ('BRS2' in file or 'BRS3' in file or 'clade' in file)): # If the gene is present in at least two strains of a comparison
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
            
            # Fill in with gaps individual comparisons where the gene is not present in both strains
            if gene == 'mhpD' and group_names[0] == comparison:
                tuple1 = tuple(fstrains[1].split('_'))
                tuple2 = tuple(fstrains[2].split('_'))
                df_add = pd.DataFrame({tuple1:[3]*aln_length, tuple2:[3]*aln_length})
                df_aln = pd.concat([df_aln, df_add], axis=1)
                print(df_aln.columns)
                
            elif gene == 'GS2' and group_names[2] == comparison:
                tuple1 = ('H1B1-05A', fstrains[0].split('_')[0])
                tuple2 = ('H1B1-05A', fstrains[0].split('_')[1])
                df_add = pd.DataFrame({tuple1:[3]*aln_length, tuple2:[3]*aln_length})
                df_aln = pd.concat([df_add, df_aln], axis=1)
                print(df_aln.columns)
            
            fstrains = get_unique_values(list(df_aln.columns))

            # pos_dict_2[gene] = pos_dict
            
        else: # Create a dataframe only with gaps for the genes that are absent
            if gene == 'S2a' or gene == 'S3':
                gene_group = 0
            elif gene == 'BRS':
                gene_group = 1
                if 'only' in file:
                    strain1 = fstrains[0].split('_')[0]
                    strain2 = fstrains[0].split('_')[1]
                    fstrains = [f'H1B105A_{strain1}', f'H1B105A_{strain2}'] + fstrains
            elif gene == 'GS2' and 'only' in file:
                gene_group = 2
            elif gene == 'GS2' and 'only' not in file:
                gene_group = 1 
                
            g_no = list(filter(lambda x: gene_order[x] == gene, gene_order))[0]
            gene_len = length_dict[group_names[gene_group]][g_no-1]
            df_aln = pd.DataFrame({fstrains[0]:[3]*gene_len, fstrains[1]:[3]*gene_len, fstrains[2]:[3]*gene_len})
            
        
        # Reverse index of genes in the forward strand
        if gene_no >= 24 or gene_no == 13:
            df_aln = df_aln.iloc[::-1] #.reset_index(drop = True)
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
        ax.patch.set(lw = 5, ec = 'black') # Border of the subplots
        gn = gene #.replace('rfbX', '?wzx').replace('wzyC', '?waaL').replace('GH39', '?GH39') # Change some of the gene annotations
        
        # Add arrows indicating gene direction
        if comparison == group_names[min(group_names.keys())]: # Add this only on top of the first row
            ax.set_title(gn, fontsize = 48, style = 'italic') # Add gene names at the top of the plot
            color = 'white' #'#e5e5e5' # Add background color to arrows
            # if gene == 'GS1':
            #     color = '#ffd9d9'
            # elif gene == 'GS2':
            #     color = '#ffd6e9'
            # elif gene == 'BRS':
            #     color = '#f8d8ff'
            # elif gene == 'S2a':
            #     color = '#fff3d7'
            # elif gene == 'S3':
            #     color = '#fffcdc'
                
            if gene_no >= 24 or gene_no == 13: # Invert arrows for genes in the forward strand
                startx = max(df_aln.index)
                dx = -max(df_aln.index)
            else: # Otherwise, plot the arrows facing toward the right
                startx = 0
                dx = max(df_aln.index)
            ax.arrow(startx, -0.5, dx, 0, width = 0.2, head_width = 0.2, 
                     head_length = 30, length_includes_head = True, 
                     clip_on = False, linewidth = 5,
                     edgecolor = 'black', facecolor = color)
            
        # Add ticks to the y axis if the gene is present in the comparison
        if count == 0:
            plt.yticks(rotation=90)
            
        nbins = ax.get_xlim()[1]//50 # Number of ticks for the x axis
        plt.tick_params(axis='x', which='major', labelsize=20) # Increase font size of ax ticks
        plt.locator_params(axis='x', nbins=nbins) # Apply change in number of ticks
        ax.set(xlabel=None) # Remove x ax label (label on the x axis of the subplot)
        ax.set(ylabel=None) # Remove y ax label

        # If the gene is not the first one, then remove the y ax ticks and labels
        if gene != 'tagU':
            plt.tick_params(left = False, labelleft = False)
            
        # Otherwise, plot the comparisons as left labels (ytick labels)
        else:
            y_ticks = get_unique_values(list(df_aln.columns))
            y_ticks = [ytick.replace('_', ' vs ') for ytick in y_ticks]
            ax.set_yticklabels(y_ticks, fontsize = 36, rotation = 'horizontal')
            if comparison == group_names[max(group_names.keys())//2]:
                ax.set_ylabel('Comparison', fontsize = 48)
        
        print(f'Added {gene} plot to comparison {comparison}!')
    
    count += 1
    print(f'Plot for {comparison} saved to {outplot}!')
    
# Add global legend
gap_patch = mpatches.Patch(color = '#F8F9F9', label = 'Gaps', ec = 'black', 
                           lw = 2)

identity_patch = mpatches.Patch(color = '#FEFDED', label = 'Identical sites', 
                                ec = 'black', lw = 2)

S_patch = mpatches.Patch(color = '#80C4E9', label = 'Synonymous sites', 
                         ec = 'black', lw = 2)

NS_patch = mpatches.Patch(color = '#FF7F3E', label = 'Non-synonymous sites', 
                          ec = 'black', lw = 2)

patches = [gap_patch, identity_patch, S_patch, NS_patch]

legend = plt.legend(handles = patches,  handlelength = 1, handleheight = 3, 
                    title = 'Site color', loc = 'upper right', framealpha = 0, 
                    frameon = False, fontsize = 36, bbox_to_anchor = (6.8, 3.38))
           
plt.setp(legend.get_title(),fontsize = 48)
legend._legend_box.align = 'left'

# Save plot
fig.text(0.5, -0.05, 'Codon position', ha='center', fontsize = 36)
plt.savefig(outplot, bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.tiff'), bbox_inches='tight')
plt.savefig(outplot.replace('.png', '.svg'), bbox_inches='tight')
