#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:41:49 2025

Python script to get the core genes for a subset of 38 strains and produce
an alignment of the amino acid sequences.

Requires Biopython and Mafft.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================
import os
from Bio import SeqIO
import subprocess

# =============================================================================
# 1. Define inputs and outputs, create paths
# =============================================================================
outpath = os.path.expanduser('~') + '/GH_project/all_core/38_strains/fasta' #Path to save outputs

representatives = {'GS1': ['A1001', 'A1202', 'A1401', 'A1805', 'K2W83', 
                           'fhon2', 'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A',
                           'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 
                           'H3B2-02X', 'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 
                           'H3B2-09X', 'H4B1-11J', 'H4B2-02J', 'H4B2-04J', 
                           'H4B2-05J', 'H4B2-11M', 'H4B4-02J', 'H4B4-05J', 
                           'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 'H4B5-03X', 
                           'H4B5-04J', 'H4B5-05J', 'LDX55', 'MP2'], 
                   'NGB': ['A0901', 'A1003', 'A1202', 'A1401', 'A1805', 
                           'K2W83', 'fhon2', 'G0101', 'H1B1-05A', 'H1B3-02M', 
                           'H3B1-01A', 'H3B1-04X', 'H3B2-02X', 'H3B2-03J', 
                           'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 'H4B1-11J', 
                           'H4B2-02J', 'H4B2-04J', 'H4B2-05J', 'H4B2-11M',
                           'H4B4-05J', 'H4B4-06M', 'H4B5-03X', 'H4B5-04J', 
                           'H4B5-05J']} #Set of strains of interest

GHs = ['NGB'] #, 'GS1']

inpath = os.path.expanduser('~') + '/GH_project/all_core/fasta' #Path to input files

if not os.path.exists(outpath):
    os.makedirs(outpath)
        
[os.makedirs(f'{outpath}/{GH}') for GH in GHs if not os.path.exists(f'{outpath}/{GH}')]
    
faas = [file for file in os.listdir(inpath) if file.endswith('faa') and 'mafft' not in file] #Raw fasta files with the protein sequences of core genes (including all strains)
fnas = [file.replace('faa', 'fna') for file in faas] #Same, but with nucleotide sequences

# =============================================================================
# 2. Create fasta files with the nucleotide and amino acid sequences of core
# genes
# =============================================================================
for file in faas: #Loop through protein fasta
    with open(f'{inpath}/{file}') as handle: #Open input file
        loop = list(SeqIO.parse(handle, 'fasta'))
        for GH in GHs:
            records = [] #Create an empty list to store the records
            for record in loop: #Loop through the records in the input file
                if any([repre.replace('-', '').upper() in record.id.upper() for repre in representatives[GH]]): #If the records are in representative strains
                    print(record.id) #Print the record
                    records.append(record) #Add the record to the list
            with open(f'{outpath}/{GH}/{file}', 'w') as out_file: #Open the output file
                SeqIO.write(records, out_file, 'fasta') #Write the records to the output file
                    
    with open(f'{inpath}/{file.replace("faa", "fna")}') as handle_fna: #Same as above, but for nucleotide sequences
        loop_fna = list(SeqIO.parse(handle_fna, 'fasta'))
        for GH in GHs:
            records_fna = []
            for record in loop_fna:
                if any([repre.replace('-', '').upper() in record.id.upper() for repre in representatives[GH]]):
                    print(record.id)
                    records_fna.append(record)
                    with open(f'{outpath}/{GH}/{file.replace("faa", "fna")}', 'w') as out_file:
                        SeqIO.write(records_fna, out_file, 'fasta')
                
# =============================================================================
# 3. Create alignments of the amino acid fasta files
# =============================================================================
for GH in GHs:
    for file in os.listdir(f'{outpath}/{GH}'): #Loop through the files generated in the previous step
        if file.endswith('.faa') and 'mafft' not in file: #Check that the input file is not an alignment
            out_aln = file.replace('.faa', '.mafft.faa') #Set the name of the output file
            subprocess.run(f'mafft-linsi --thread 12 {outpath}/{GH}/{file} > {outpath}/{GH}/{out_aln}', 
                           shell = True) #Create the alignments using 12 threads
