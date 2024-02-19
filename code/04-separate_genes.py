#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 16:43:13 2022

Script to divide GHs by type and subset of strains of interest.

@author: Marina Mota Merlo

"""

# =============================================================================
# Import modules
# =============================================================================

import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import as_fasta
import logging, traceback

# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# Define Snakemake input
# =============================================================================

GS1 = snakemake.params.GS1
GS2 = snakemake.params.GS2
BRS = snakemake.params.BRS
short = snakemake.params.short
NGB = snakemake.params.NGB
S1 = snakemake.params.S1
S2a = snakemake.params.S2a
S2b = snakemake.params.S2b
S3 = snakemake.params.S3

input_files = snakemake.input
GH70s = snakemake.params.GH70s
GH32s = snakemake.params.GH32s
representatives = snakemake.params.repr
subset = snakemake.params.subset

# =============================================================================
# Define path to outputs and dictionary of variables
# =============================================================================

path = 'data/fasta'

myVars = locals()

# =============================================================================
# Function to create file if it does not exist, otherwise append to file
# =============================================================================

def create_write(path, prefix, filename):
    with open(f'{path}/{prefix}/{filename}', 'a+') as outfile:
        outfile.write(as_fasta(record))
        
# =============================================================================
# Loop through inputs and save the information to outputs
# =============================================================================

for file in input_files: #Loop through Snakemake inputs
    if 'GH70' in file: #If GH70 is in the file
        gene_list = GH70s #Get all GH70 types
        prefix = 'GH70' #Set GH70 as prefix for outputs
    elif 'GH32' in file:
        gene_list = GH32s
        prefix = 'GH32'
    if 'complete' in file:
        add = 'complete_'
    else: add = ''
    with open(file) as seq_file: #Open input file
        for record in SeqIO.parse(seq_file, 'fasta'): #Read it as fasta
            check = False #To check if the strain is in the subset
            check2 = True #To check that GH70 domain is complete
            sname = record.id.split('_')[0] #Get strain name from locus tag
            if sname[0] == 'H':
                sname = sname[:4] + '-' + sname[4:]
            elif sname == 'APS55':
                sname = 'MP2'
            elif sname == 'LDX55':
                sname = 'IBH001'
            elif sname == 'K2W83':
                sname = 'DSMZ12361'
            elif sname == 'FHON2':
                sname = 'Fhon2'
                
            filename = f'{add}{prefix}_all.{file[-3:]}' #Define output file name
            with open(f'{path}/{prefix}/{filename}', 'a+') as outfile: #Create the file if it does not exist, otherwise append to file
                outfile.write(as_fasta(record))
            print(f'{record.id} added to {path}/{prefix}/{filename}')
            
            for gene_type in gene_list: #Loop through gene types
                if record.id in myVars[gene_type]: #Retrieve locus tags per type
                    if gene_type not in ['GS1', 'GS2', 'BRS', 'NGB'] or (file.endswith('.faa') and len(record.seq) > 700) or (file.endswith('.fna') and len(record.seq) > 2100):
                        filename = f'{add}{gene_type}_all.{file[-3:]}' #File for a particular GH type
                        create_write(path, prefix, filename)
                        print(f'{record.id} added to {path}/{prefix}/{filename}')
                    
                        if check:
                            filename = f'{add}{gene_type}_subset.{file[-3:]}' #File for a particular GH type in the subset
                            create_write(path, prefix, filename)
                            print(f'{record.id} added to {path}/{prefix}/{filename}')
            
            if sname in representatives: #Get strain name and check if it is among the representatives
                if sname in subset:
                    check = True
                filename = f'{add}{prefix}_repset.{file[-3:]}' #Define output file name
                with open(f'{path}/{prefix}/{filename}', 'a+') as outfile: #Create the file if it does not exist, otherwise append to file
                    outfile.write(as_fasta(record))
                print(f'{record.id} added to {path}/{prefix}/{filename}')
                
                if check: 
                    filename = f'{add}{prefix}_subset.{file[-3:]}'
                    create_write(path, prefix, filename)
                    print(f'{record.id} added to {path}/{prefix}/{filename}')
                    
                for gene_type in gene_list: #Loop through gene types
                    if record.id in myVars[gene_type]: #Retrieve locus tags per type
                        if gene_type not in ['GS1', 'GS2', 'BRS', 'NGB'] or (file.endswith('.faa') and len(record.seq) > 700) or (file.endswith('.fna') and len(record.seq) > 2100):
                            filename = f'{add}{gene_type}_repset.{file[-3:]}' #File for a particular GH type
                            create_write(path, prefix, filename)
                            print(f'{record.id} added to {path}/{prefix}/{filename}')
                        
                            if check:
                                filename = f'{add}{gene_type}_subset.{file[-3:]}' #File for a particular GH type in the subset
                                create_write(path, prefix, filename)
                                print(f'{record.id} added to {path}/{prefix}/{filename}')
                    
                
                
                    
                
        
