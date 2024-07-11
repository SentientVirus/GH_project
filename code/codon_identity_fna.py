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

# OBS! Add the grey to the absent comparisons of mhpD too!
import os, subprocess
import logging, traceback, sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from Bio import SeqIO
import matplotlib.gridspec as gridspec
# =============================================================================
# 0. Logging
# =============================================================================
log = os.path.expanduser('~') + '/GH_project/logs/dS_matrix.log'

if os.path.exists(log):
    os.remove(log)

if not os.path.exists(os.path.dirname(log)):
    os.makedirs(os.path.dirname(log))

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

sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 0. Input and function definition
# =============================================================================
aa_path = os.path.expanduser('~') + '/GH_project/data/fasta'
gene_types = ['GH70', 'GH32', 'other_genes']
neighbor = 'other_genes'
seq_outdir = f'{aa_path}/{neighbor}/subset'
aln_path = f'{seq_outdir}/alignments'
codon_path = os.path.expanduser('~') + '/GH_project/data/codons/other_genes/subset'
outdir = os.path.expanduser('~') + '/GH_project/plots/substitution_subsets'
thr = -1
pal2nal_path = os.path.expanduser('~') + '/GH_project/pal2nal.v14/pal2nal.pl'

strain_groups = {'H3B2-03M': 0, 'H4B4-02J': 0, 'H4B5-03X': 0, 'H1B3-02M': 1,
                 'H3B2-09X': 1, 'H4B5-05J': 1, 'MP2': 2, 'H3B2-02X': 2,
                 'H3B2-03J': 2}

group_names = {0: 'GS1_S2-3_subset', 1: 'GS1-2_BRS', 2: 'only_GS1+GS2'}

gene_order = {1: 'ohrR', 2: 'yifK', 3: 'yifK2', 4: 'yhdG', 5: 'GH39', 6: 'nox',
              7: 'ydiL', 8: 'S2a', 9: 'GS1', 10: 'BRS', 11: 'GS2', 12: 'mhpD',
              13: 'oppA', 14: 'S3', 15: 'rfbX', 16: 'sbnD', 17: 'wzyC', 
              18: 'branch', 19: 'branch2', 20: 'epsE', 21: 'GTB', 22: 'epsL',
              23: 'ywqE', 24: 'ywqD', 25: 'ywqC', 26: 'tagU'}

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

GH_genes = ['GS1', 'GS2', 'BRS', 'S2a', 'S3']

paths2make = [aln_path, codon_path, outdir]


# =============================================================================
# 2. Function definitions
# =============================================================================
def make_paths(my_path):
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
    for nbr in os.listdir(input_dir):
        check = False
        if 'other_genes' in input_dir and nbr.startswith('a_kunkeei') and 'mafft' not in nbr:
            gene_name = nbr.split('_')[2].split('.')[0]
            check = True
        elif 'GH' in input_dir and any(gh_gene in nbr for gh_gene in GH_genes) and 'repset' in nbr and 'mafft' not in nbr:
            gene_name = nbr.split('_')[0]
            check = True
        if check:
            if nbr.endswith('.faa'):
                suffix = 'faa'
            elif nbr.endswith('.fna'):
                suffix = 'fna'
            print(f'Retrieving gene {gene_name}')
            with open(f'{input_dir}/{nbr}') as nbr_file:
                seqs = SeqIO.parse(nbr_file, 'fasta')
                for seq in seqs:
                    seq.id = seq.id.replace('APS', '')
                    seq.id = seq.id.replace('55_RS', 'MP2_')
                    print(seq.id)
                    seq.description = seq.description.replace('55_RS', 'MP2_')
                    strain = seq.id.split('_')[0]
                    if strain[0] == 'H':
                        strain = strain[:4] + '-' + strain[4:]
                    if strain in strain_groups.keys():
                        prefix = group_names[strain_groups[strain]]
                        with open(f'{seq_outdir}/{prefix}_{gene_name}.{suffix}', 'a+') as seqfile:
                            print(f'Writing gene {seq.id} to {seq_outdir}/{prefix}_{gene_name}.{suffix}')
                            SeqIO.write(seq, seqfile, 'fasta')
                            
def get_unique_values(tuple_list):
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
        
        if file in os.listdir(aln_path) and not ('only' in file and 'BRS' in file):
            print(f'Retrieving lengths from {aln_path}/{file}')
        
            with open(f'{aln_path}/{file}') as handle:
                f = SeqIO.parse(handle, 'fasta')
                for record in f:
                    aln_length = len(record.seq)//3
                    length_dict[comparison].append(aln_length)
                    break
        else: length_dict[comparison].append(0)
        
absent_genes = {}
for comparison in group_names.values():
    for gene in ['S2a', 'GS2', 'BRS', 'S3']:
        gene_number = list(filter(lambda x: gene_order[x] == gene, gene_order))[0]
        if length_dict[comparison][gene_number-1] != 0:
            absent_genes[gene] = (gene_number, length_dict[comparison][gene_number-1])
            
for comparison in group_names.values():
    for gene in ['S2a', 'GS2', 'BRS', 'S3']:
        gene_number = absent_genes[gene][0]
        if length_dict[comparison][gene_number-1] == 0:
            length_dict[comparison][gene_number-1] = absent_genes[gene][1]
        

# Create a dataframe assigning different values to identical codons (0),
# synonymous substitutions (1), non-synonymous substitutions (2) and gaps (3)
df_aln_dict = {}
pos_dict_2 = {}
fig = plt.figure(constrained_layout=False, figsize=(160, 30))
spec = gridspec.GridSpec(nrows=3, ncols=26, figure=fig, width_ratios=length_dict[comparison])
count = 0
for comparison in group_names.values():
    plt.tight_layout(h_pad = 2, w_pad = 2) # Adjust space between genes
    # fig = plt.figure(constrained_layout=False, figsize=(160, 10))
    # spec = gridspec.GridSpec(nrows=1, ncols=26, figure=fig, width_ratios=length_dict[comparison])
    for gene_no in range(1, len(gene_order)+1):
        gene = gene_order[gene_no]
        file = f'{comparison}_{gene}.pal2nal.fna' # Input alignment
        
        if file in os.listdir(aln_path) and not ('only' in file and 'BRS' in file):
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

            pos_dict_2[gene] = pos_dict
            
        else:
            if gene == 'S2a' or gene == 'S3':
                gene_group = 0
            elif gene == 'BRS':
                gene_group = 1
            g_no = list(filter(lambda x: gene_order[x] == gene, gene_order))[0]
            gene_len = length_dict[group_names[gene_group]][g_no-1]
            df_aln = pd.DataFrame({fstrains[0]:[3]*gene_len, fstrains[1]:[3]*gene_len, fstrains[2]:[3]*gene_len})
            
        
        df_aln.index += 1 # Make the index start with 1, not 0
        df_aln_dict[gene] = df_aln # Save dataframe to dictionary
        # Assign colors to the values in the dataframe
        cmap = []
        if 0 in df_aln_dict[gene].values:
            cmap.append('white')
        if 1 in df_aln_dict[gene].values:
            cmap.append('#67e3ff')
        if 2 in df_aln_dict[gene].values:
            cmap.append('#0079d8')
        if 3 in df_aln_dict[gene].values:
            cmap.append('#9a9a9a')
        
        # Create a plot for the gene and comparison
        ax = fig.add_subplot(spec[count, gene_no-1])
        h3 = sns.heatmap(df_aln_dict[gene].T, cmap = cmap, cbar=False, ax=ax)
        ax.patch.set(lw=5, ec='black')
        ax.set_ylabel('Comparison', fontsize = 24)
        ax.set_xlabel('Codon position', fontsize = 24)
        ax.set_title(gene, fontsize = 30)
        nbins = ax.get_xlim()[1]//25
        if count == 0:
            plt.yticks(rotation=90)
        plt.tick_params(axis='x', which='major', labelsize=20) # Increase font size of ax ticks
        plt.tick_params(axis='y', which='major', labelsize=11)
        plt.locator_params(axis='x', nbins=nbins) 
        print(f'Added {gene} plot to comparison {comparison}!')
    
    count += 1
    print(f'Plot for {comparison} saved to {outdir}/matrix_dNdSplot.png!')
    # Save plot
plt.savefig(f'{outdir}/matrix_dNdSplot.png', bbox_inches='tight')
