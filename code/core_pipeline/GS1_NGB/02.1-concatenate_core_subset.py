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

# =============================================================================
# 1. Define inputs and outputs
# =============================================================================
threads = 48
GHs = ['NGB'] #['GS1', 'NGB']

for GH in GHs:
    inpath = os.path.expanduser('~') + f'/GH_project/all_core/38_strains/fasta/{GH}'
    fasta_out = f'{inpath}/concatenated.faa'
    log = os.path.expanduser('~') + f'/GH_project/logs/iqtree/38_strains_{GH}_core_iqtree.log'
    outdir = inpath.replace('fasta', 'trees')
    aln_dir = inpath.replace('fasta', 'alignment')
    
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
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
    #All nuclear and general matrices are considered -mset DCMut,JTT,JTTDCMut,VT,LG,WAG,Q.pfam,PMB
    subprocess.run(f'iqtree2 -nt AUTO -ntmax {threads} -p {aln_dir} --prefix 38strains_{GH} -st AA -m MFP -msub nuclear -bb 1000 -bnni -redo > {log}', shell = True)
    subprocess.run(f'mkdir -p {outdir}', shell = True) #Create the output directory if it doesn't exist
    subprocess.run(f'mv 38strains_{GH}.* {outdir}', shell = True) #Move the IQtree output to the desired directory
