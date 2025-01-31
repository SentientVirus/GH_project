#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:41:49 2025

Python script to get the core genes for a subset the subset of strains which
contains each type of GH gene and produce an alignment of the amino acid 
sequences.

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
workdir = os.path.expanduser('~') + '/GH_project'
inpath = f'{workdir}/all_core/fasta' #Path to input files

# gene_types = ['GS1', 'GS2', 'BRS', 'NGB', 'S2a', 'S3']

gene_dict = {'GS1': ['A1001', 'A1202', 'A1401', 'A1805', 'K2W83', 'fhon2', 
                     'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A', 'H1B3-02M', 
                     'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 'H3B2-02X', 
                     'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 
                     'H4B1-11J', 'H4B2-02J', 'H4B2-04J', 'H4B2-05J', 
                     'H4B2-11M', 'H4B4-02J', 'H4B4-05J', 'H4B4-06M', 
                     'H4B4-12M', 'H4B5-01J', 'H4B5-03X', 'H4B5-04J', 
                     'H4B5-05J', 'LDX55', 'MP2'],
             
             'GS2': ['K2W83', 'G0403', 'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 
                     'H3B1-04X', 'H3B2-02X', 'H3B2-03J', 'H4B2-02J', 
                     'H4B2-04J', 'H4B5-04J', 'H4B5-05J', 'LDX55', 'MP2'],
             
             'BRS': ['A1401', 'K2W83', 'fhon2', 'G0403', 'H1B3-02M', 
                     'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 'H3B2-02X', 
                     'H3B2-03J', 'H4B2-02J', 'H4B2-04J', 'H4B5-04J', 
                     'H4B5-05J', 'LDX55'],
             
             'NGB': ['A0901', 'A1003', 'A1202', 'A1401', 'A1805', 'K2W83', 'fhon2', 
                         'G0101', 'H1B1-05A', 'H1B3-02M', 'H3B1-01A', 'H3B1-04X', 
                         'H3B2-02X', 'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 
                         'H4B1-11J', 'H4B2-04J', 'H4B2-05J', 'H4B2-11M', 'H4B4-05J', 
                         'H4B4-06M', 'H4B5-03X', 'H4B5-04J', 'H4B5-05J'],
             
             'S2a': ['A0901', 'A1001', 'A1805', 'H1B1-04J', 'H3B2-03M',
                     'H4B2-06J', 'H4B4-02J', 'H4B4-12M', 'H4B5-03X'],
             
             'S3': ['A0901', 'A1001', 'A1404', 'H3B2-03M', 'H4B2-06J', 
                    'H4B4-02J', 'H4B4-12M', 'H4B5-03X']
             }
    
faas = [file for file in os.listdir(inpath) if file.endswith('faa') and 'mafft' not in file] #Raw fasta files with the protein sequences of core genes (including all strains)
fnas = [file.replace('faa', 'fna') for file in faas] #Same, but with nucleotide sequences

# =============================================================================
# 2. Create fasta files with the nucleotide and amino acid sequences of core
# genes
# =============================================================================
for gtype in gene_dict.keys():
    outpath = os.path.expanduser('~') + f'/GH_project/all_core/{gtype}/fasta' #Path to save outputs
    
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    for file in faas: #Loop through protein fasta
        with open(f'{inpath}/{file}') as handle: #Open input file
            records = [] #Create an empty list to store the records
            for record in SeqIO.parse(handle, 'fasta'): #Loop through the records in the input file
                if any([repre.replace('-', '').upper() in record.id.upper() for repre in gene_dict[gtype]]): #If the records are in representative strains
                    print(record.id) #Print the record
                    records.append(record) #Add the record to the list
            with open(f'{outpath}/{file}', 'w') as out_file: #Open the output file
                SeqIO.write(records, out_file, 'fasta') #Write the records to the output file
                        
        with open(f'{inpath}/{file.replace("faa", "fna")}') as handle_fna: #Same as above, but for nucleotide sequences
            records_fna = []
            for record in SeqIO.parse(handle_fna, 'fasta'):
                if any([repre.replace('-', '').upper() in record.id.upper() for repre in gene_dict[gtype]]):
                    print(record.id)
                    records_fna.append(record)
                    with open(f'{outpath}/{file.replace("faa", "fna")}', 'w') as out_file:
                        SeqIO.write(records_fna, out_file, 'fasta')
                
# =============================================================================
# 3. Create alignments of the amino acid fasta files
# =============================================================================
    for file in os.listdir(outpath): #Loop through the files generated in the previous step
        if file.endswith('.faa') and 'mafft' not in file: #Check that the input file is not an alignment
            out_aln = file.replace('.faa', '.mafft.faa') #Set the name of the output file
            in_fna = file.replace('.faa', '.fna')
            out_pal2nal = file.replace('.faa', '.pal2nal')
            outpath2 = outpath.replace('fasta', 'codons')
            
            if not os.path.exists(outpath2):
                os.makedirs(outpath2)
                
            subprocess.run(f'mafft-linsi --thread 12 {outpath}/{file} > {outpath}/{out_aln}', 
                           shell = True) #Create the alignments using 12 threads
            subprocess.run(f'{workdir}/pal2nal.v14/pal2nal.pl {outpath}/{out_aln} {outpath}/{in_fna} -output paml -nogap > {outpath2}/{out_pal2nal};',
                           shell = True)
