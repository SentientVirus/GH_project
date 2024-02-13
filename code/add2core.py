#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 17:10:01 2023

This is a script to add the genes from strain IBH001 to the files with a 
set of core genes.

@author: Marina Mota Merlo
"""
import os
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

indir = 'core_genes'
outdirs = ['untrimmed_core', 'trimmed_core']
match_file = 'IBH001_tags.tsv'
cds_file = '../Akunkeei_files/cds/IBH001_cds_from_genomic.fna'
record_dict = {}

# =============================================================================
# Create output directories if they don't exist
# =============================================================================

for directory in outdirs:
    if not os.path.exists(directory):
            os.makedirs(directory)

# =============================================================================
# Define functions to be used in the code
# =============================================================================
def trimal_gappyout(in_file, out_file):
    # Trim alignments using TrimAl with the -gappyout option
    command = f'trimal -in {in_file} -out {out_file} -gappyout;'
    # Call TrimAl using os
    os.system(command)
    
def mafft_linsi(extra_seq, unaligned, aligned, log):
    with open(log, 'w') as logfile:
        logfile.write('')
    command = f'GREPDB="mafft-linsi --thread 64 --add {extra_seq} {unaligned} > {aligned} 2>> {log}"; /bin/bash -c "$GREPDB"'
    os.system(command)


# =============================================================================
# 1. Create a dictionary with all the SeqRecords of IBH001 matched to the
# A1401 IDs
# =============================================================================
with open(cds_file) as cds:
    cds_obj = SeqIO.parse(cds, 'fasta')
    gen_obj = sorted(list(cds_obj), key = lambda a: a.description.split('locus_tag=')[1].split(']')[0])
    
with open(match_file) as matches:
    df = pd.read_csv(match_file, sep = '\t')
    for row in df.iterrows():
        for record in gen_obj:
            if row[1][1] in record.description:
                record.id = record.description.split('locus_tag=')[1].split(']')[0]
                record_dict[row[1][0]] = record

# =============================================================================
# 2. Create files with the genes I want to add to the sequence alignment and
# align them to the original alignments
# =============================================================================
for file in os.listdir(indir):
    if file.endswith('genes.fasta'):
        ref_tag = file.split('_genes')[0]
        add_sq = record_dict[ref_tag]
        add_sq.description = ''
        fna_file = 'temp.fna'
        faa_file = 'temp.faa'
        with open(fna_file, 'w') as fna, open(faa_file, 'w') as faa:
            SeqIO.write(add_sq, fna, 'fasta')
            add_sq.seq = add_sq.seq.translate(to_stop = True)
            SeqIO.write(add_sq, faa, 'fasta')
        prot_file = file.replace('genes', 'prot')
        # with open(f'{indir}/{file}') as genes, open(f'{indir}/{prot_file}') as prots, open('temp.fna', 'a') as fna, open('temp.faa', 'a') as faa:
        #     gene_read = SeqIO.parse(genes, 'fasta')
        #     prot_read = SeqIO.parse(prots, 'fasta')
        #     SeqIO.write(gene_read, fna, 'fasta')
        #     SeqIO.write(prot_read, faa, 'fasta')
        outfile = f'{outdirs[0]}/{file}'
        prot_out = f'{outdirs[0]}/{prot_file}'
        mafft_linsi(fna_file, f'{indir}/{file}', outfile, f'{outdirs[0]}/{file.split(".")[0]}.log')
        mafft_linsi(faa_file, f'{indir}/{prot_file}', prot_out, f'{outdirs[0]}/{prot_file.split(".")[0]}.log')
        
# =============================================================================
# 3. Call trimal from Python to trim the alignment (using os.shell)
# =============================================================================

        trimal_gappyout(outfile, outfile.replace('un', ''))
        trimal_gappyout(prot_out, prot_out.replace('un', ''))
