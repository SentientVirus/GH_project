#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 17:00:37 2023

@author: marina
"""
import os
from Bio import SeqIO

infolder = os.path.expanduser('~') + '/Akunkeei_files/fna'
outfolder = infolder + '/reverse'

infiles = os.listdir(infolder)
for file in infiles:
    if os.path.isfile(f'{infolder}/{file}') and 'MP2' not in file:
        print(file)
        seqs = SeqIO.parse(f'{infolder}/{file}', 'fasta')
        new_seqs = []
        for fna in seqs:
            fna.seq = fna.seq.reverse_complement()
            new_seqs.append(fna)
            
            with open(f'{outfolder}/{file}', 'w') as outfile:
                SeqIO.write(new_seqs, outfile, 'fasta')
    elif os.path.isfile(f'{infolder}/{file}') and 'MP2' in file:
        print(file)
        seqs = SeqIO.parse(f'{infolder}/{file}', 'fasta')
        new_seqs = []
        for fna in seqs:
            fna.seq = fna.seq.reverse_complement()
            new_seqs.append(fna)
            
            with open(f'{infolder}/{file}', 'w') as outfile:
                SeqIO.write(new_seqs, outfile, 'fasta')
