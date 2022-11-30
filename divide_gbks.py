#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 10:01:24 2022

@author: marina
"""
# =============================================================================
# Script to divide multi-record GenBank files into single-record files.
# =============================================================================

import os

indir = 'data/genbanks'
outdir = 'gbks'
newfile = {}

#Function to parse GenBank files, it outputs record names and record info
def parse_GenBank(file):
    if file.endswith('.gbk'):
        with open(file, 'r') as handle:
            f = (line for line in handle)
            for line in f:
                if 'LOCUS' in line:
                    locus_name = line.split()[1]
                    newfile = ''
                    print(locus_name)
                newfile += line
    return locus_name, newfile

#Results are saved into a dictionary
for file in sorted(os.listdir(indir)):
    locus_name, new_file = parse_GenBank(f'{indir}/{file}')
    newfile[locus_name] = new_file

#Individual records are saved in separate files                
for locus_name in newfile.keys():
    with open(f'{outdir}/{locus_name}.gbk', 'w') as outgbk:
        outgbk.write(newfile[locus_name])
