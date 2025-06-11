# -*- coding: utf-8 -*-
"""
Created on Mon May  12 17:16 2025

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================
import os, sys
import logging, traceback
from Bio import SeqIO
import subprocess

# =============================================================================
# 0. Logging
# =============================================================================
workdir = os.path.expanduser('~') + '/GH_project/motifs'
logfile = f'{workdir}/logs/01-get_complete_codon_aln.log'

if not os.path.exists(os.path.dirname(logfile)):
    os.makedirs(os.path.dirname(logfile))
    
with open(logfile, 'w') as log_out:
    log_out.write('')

logging.basicConfig(filename = logfile, level = logging.INFO,
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

sys.stdout = open(logfile, 'a')

# =============================================================================
# 1. Define inputs
# =============================================================================
indir = os.path.expanduser('~') + '/GH_project/data/fasta/GH70'
in_file = f'{indir}/complete_GH70_repset.faa'
aln_file = in_file.replace('faa', 'mafft.faa')
fna_file = in_file.replace('faa', 'fna')
pal2nal_path = os.path.expanduser('~') + '/GH_project/pal2nal.v14'
codon_file = in_file.replace('fasta/GH70', 'codons').replace('.faa', '_codon.fna')
threads = 12

subprocess.run(f'mafft-linsi --thread {threads} {in_file} > {aln_file} 2>> {logfile}', shell = True)
subprocess.run(f'{pal2nal_path}/pal2nal.pl {aln_file} {fna_file} -output fasta > {codon_file} 2>> {logfile}', shell = True)