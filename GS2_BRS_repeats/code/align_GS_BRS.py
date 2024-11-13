#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:31:16 2024

@author: Marina Mota-Merlo

This script retrieves the sequences of GS2_BRS genes and those of G0403 GS2 and
BRS and A1401 BRS and performs an alignment of all the sequences. Requires
Mafft and Biopython and PyMSAViz
"""

import logging, traceback, sys
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pymsaviz import MsaViz

# =============================================================================
# 0. Logging
# =============================================================================
log = os.path.expanduser('~') + '/GH_project/GS2_BRS_repeats/logs/align_sequences.log'

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
# 1. Defining inputs and functions
# =============================================================================
strains_of_interest = ['A1401', 'DSMZ12361', 'H3B1-01J', 'H3B1-04J', 'H4B2-04J', 
                       'G0403', 'H4B2-02J', 'IBH001']
genes_of_interest = ['A1401_12760', 'G0403_13110', 'G0403_13120', 
                     'H3B101J_13010', 'H3B104J_13020', 'H4B202J_12890', 
                     'H4B204J_13350', 'IBH001_06335', 'DSM_06185']
types = {'AKUH3B101J_13010': 'GS2_BRS, alpha-1,3', 
         'AKUA1401_12760': 'BRS, incomplete alpha-1,3',
         'AKUH3B104J_13020': 'GS2_BRS, alpha-1,3',
         'AKUG0403_13110': 'BRS, alpha-1,2', 'AKUG0403_13120': 'GS2',
         'AKUH4B204J_13350': 'GS2_BRS, alpha-1,3', 
         'AKUH4B202J_12890': 'GS2_BRS, alpha-1,2', 
         'IBH001_06335': 'GS2_BRS, alpha-1,2',
         'DSM_06185': 'GtfZ, alpha-1,3'}


faa = os.path.expanduser('~') + '/Akunkeei_files/faa'
out_faa = os.path.expanduser('~') + '/GH_project/GS2_BRS_repeats/files/input_seqs.faa'
out_mafft = out_faa.replace('.faa', '.mafft.faa')
out_aln_fig = os.path.expanduser('~') + '/GH_project/GS2_BRS_repeats/plots/GS2_BRS_aln.png'

out_faa_path = os.path.dirname(out_faa)
log_path = os.path.dirname(log)
fig_path = os.path.dirname(out_aln_fig)
paths = [out_faa_path, log_path, fig_path]

threads = 12

[os.makedirs(file_path) for file_path in paths if not os.path.exists(file_path)]
    
with open(out_faa, 'w+') as f:
    f.write('')
    
with open(log, 'w+') as flog:
    flog.write('')

# =============================================================================
# 2. Save input sequences to fasta file
# =============================================================================
for strain in strains_of_interest:
    # if any([strain in file for strain in strains_of_interest]):
    with open(f'{faa}/{strain}_protein.faa') as strain_faa:
        print(f'Reading file {faa}/{strain}_protein.faa')
        fr = SeqIO.parse(strain_faa, 'fasta')
        for seq in fr:
            locus_tag = ''
            description_list = seq.description.split(' ')
            for item in description_list:
                if ('AKU' in item) or ('WP_220382103.1' in item) or ('UZX33033.1' in item):
                    locus_tag = item
                    locus_tag = locus_tag.replace('UZX33033.1', 'IBH001_06335')
                    locus_tag = locus_tag.replace('WP_220382103.1', 'DSM_06185')
            if any([gene in locus_tag for gene in genes_of_interest]):
                print(f'Retrieving gene with locus tag {locus_tag}...')
                gene_name = types[locus_tag].split(',')[0]
                description = types[locus_tag]
                locus_tag = locus_tag.replace('AKU', '')
                record = SeqRecord(seq.seq, id = locus_tag, name = gene_name,
                                   description = description)
                print(f'Creating record with id {locus_tag} and length {len(seq.seq)}')
                with open(out_faa, 'a+') as f:
                    SeqIO.write(record, f, 'fasta')
            
            

# =============================================================================
# 3. Run Mafft-linsi and save results
# =============================================================================
with open(out_faa) as fasta_seqs:
    subprocess.run(f'mafft-linsi --thread {threads} {out_faa} > {out_mafft} 2>> {log};', shell = True)
    
# =============================================================================
# 4. Generate plot
# =============================================================================
mv = MsaViz(out_mafft, wrap_length=120, show_count=True)
mv.savefig(out_aln_fig)
mv.savefig(out_aln_fig.replace('png', 'svg')
mv.savefig(out_aln_fig.replace('png', 'tiff')
