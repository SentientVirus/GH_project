#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:49:25 2025

Python script to make a strain phylogeny based on core genes for any set of
strains.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================
from Bio import SeqIO
import os, subprocess
# from Bio.Nexus import Nexus

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================
inpath = os.path.expanduser('~') + '/GH_project/all_core/38_strains/fasta'
# outpath = inpath.replace('fasta', 'nexus')
# nexus_out = f'{outpath}/concatenated.nex'
fasta_out = f'{inpath}/concatenated.faa'
log = os.path.expanduser('~') + '/GH_project/logs/iqtree/38_strains_core_iqtree.log'
outdir = inpath.replace('fasta', 'trees')
aln_dir = inpath.replace('fasta', 'alignment')
threads = 48

if not os.path.exists(aln_dir):
    os.makedirs(aln_dir)

infiles = [f'{inpath}/{file}' for file in os.listdir(inpath) if file.endswith('.mafft.faa')]

# =============================================================================
# 2. Loop through input files to generate the right inputs for IQtree
# =============================================================================
for file in infiles:
    outfile = file.replace(inpath, aln_dir)
    records = []
    with open(file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            record.id = record.id.split('_')[0].replace('K2W83', 'DSM').replace('LDX55', 'IBH001').replace('FHON', 'fhon')
            if record.id.startswith('H'):
                record.id = record.id[:4] + '-' + record.id[4:]
            record.description = record.id
            records.append(record)
            print(record.id)
    with open(outfile, 'w') as out_file:
        SeqIO.write(records, out_file, 'fasta')
            
# =============================================================================
# 3. Run IQtree
# =============================================================================
#Run IQtree2 on the folder with input alignments, restricting the possible models
#to a subset and generating 1000 UFBootstrap replicates
subprocess.run(f'iqtree2 -nt AUTO -ntmax {threads} -p {aln_dir} --prefix 38strains -st AA -mset Q.pfam,LG,WAG,JTT -bb 1000 -bnni > {log}', shell = True)
subprocess.run(f'mkdir -p {outdir}', shell = True) #Create the output directory if it doesn't exist
subprocess.run(f'mv 38strains.* {outdir}', shell = True) #Move the IQtree output to the desired directory

## Below is the other script that I tried to concatenate the alignments before making the tree
# import pandas as pd

# def fasta_reader(file):
#     fasta_df = pd.read_csv(file, sep='>', lineterminator='>', header=None)
#     fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split('\n', 1, \
#                                         expand=True)
#     new_acc = fasta_df['Accession'].to_list()
#     new_accs = [acc.split(' ')[0].split('_')[0].replace('K2W83', 'DSM').replace('LDX55', 'IBH001').replace('FHON', 'fhon') for acc in new_acc]
#     new_accs = [acc[:4] + '-' + acc[4:] if acc.startswith('H') else acc for acc in new_accs]
#     fasta_df['Accession'] = new_accs
#     fasta_df.drop(0, axis=1, inplace=True)
#     fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True)
#     return fasta_df

# df = pd.concat(fasta_reader(i) for i in infiles)
# agg_func = {'Sequence': 'sum'}
# df_new = df.groupby(df['Accession']).aggregate(agg_func)

# acc = [f'>{strain}' for strain in list(df_new.index)]
# df_new = pd.concat([pd.Series(acc, index=df_new.index, name='Accession'), df_new], axis=1)

# df_new.to_csv(fasta_out, sep='\n', index=None, header=None)

# subprocess.run(f'iqtree -nt AUTO -ntmax {threads} -s {fasta_out} -st AA -mset Q.pfam,LG,WAG,JTT -bb 1000 -bnni > {log}', shell = True)
# subprocess.run(f'mkdir -p {outdir}', shell = True)
# subprocess.run(f'mv {fasta_out}.* {outdir}', shell = True)
