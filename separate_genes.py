#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 16:43:13 2022

@author: marina
"""
import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import as_fasta

GS1 = snakemake.params.GS1
GS2 = snakemake.params.GS2
BRS = snakemake.params.BRS
short = snakemake.params.short
NCB = snakemake.params.NCB
S1 = snakemake.params.S1
S2a = snakemake.params.S2a
S2b = snakemake.params.S2b
S3 = snakemake.params.S3

input_files = snakemake.input
GH70s = snakemake.params.GH70s
GH32s = snakemake.params.GH32s
representatives = snakemake.params.repr
subset = snakemake.params.subset

path = 'data/fasta'

myVars = locals()
# with open('try.txt', 'w') as file:
#     file.write('')

# for GH70 in GH70s:
#     with open('try.txt', 'a') as file:
#         file.write(GH70)

#record_dict = {}

for file in input_files:
    if 'GH70' in file:
        gene_list = GH70s
        prefix = 'GH70'
    elif 'GH32' in file:
        gene_list = GH32s
        prefix = 'GH32'
    with open(file) as seq_file:
        for record in SeqIO.parse(seq_file, 'fasta'):
            for gene_type in gene_list:
                if record.id.split('_')[0] in representatives:
                    filename = f'{prefix}_repset.{file[-3:]}'
                    if filename not in os.listdir(f'{path}/{prefix}'):
                        mode = 'w'
                    else:
                        mode = 'a'
                    with open(f'{path}/{prefix}/{filename}', mode) as outfile:
                        outfile.write(as_fasta(record))
                    print(f'{record.id} added to {path}/{prefix}/{filename}')
                    if record.id.split('_')[0] in subset:
                        filename = f'{prefix}_subset.{file[-3:]}'
                        if filename not in os.listdir(f'{path}/{prefix}'):
                            mode = 'w'
                        else:
                            mode = 'a'
                        with open(f'{path}/{prefix}/{filename}', mode) as outfile:
                            outfile.write(as_fasta(record))
                        print(f'{record.id} added to {path}/{prefix}/{filename}')
                    if record.id in myVars[gene_type]:
                        filename = f'{gene_type}_repset.{file[-3:]}'
                        if filename not in os.listdir(f'{path}/{prefix}'):
                            mode = 'w'
                        else:
                            mode = 'a'
                        with open(f'{path}/{prefix}/{filename}', mode) as outfile:
                            outfile.write(as_fasta(record))
                        print(f'{record.id} added to {path}/{prefix}/{filename}')
                        if record.id.split('_')[0] in subset:
                            filename = f'{gene_type}_subset.{file[-3:]}'
                            if filename not in os.listdir(f'{path}/{prefix}'):
                                mode = 'w'
                            else:
                                mode = 'a'
                            with open(f'{path}/{prefix}/{filename}', mode) as outfile:
                                outfile.write(as_fasta(record))
                            print(f'{record.id} added to {path}/{prefix}/{filename}')
                    
                
                
                    
                
        