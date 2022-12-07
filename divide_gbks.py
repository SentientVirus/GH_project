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
new_dict = {}
newfile = {}

#Function to parse GenBank files, it outputs record names and record info
def parse_GenBank(file, locus_dict):
    with open(file, 'r') as handle:
        f = (line for line in handle)
        check = False
        for line in f:
            if 'LOCUS' in line:
                locus_name = line.split()[1]
                locus_dict[locus_name] = ''
                print(locus_name)
            elif '_1' not in locus_name and 'ORIGIN' in line:
                check = True
            elif '//' in line:
                check = False
            if check:
                line = '  ' + line
            locus_dict[locus_name] += line
    return locus_dict

#Results are saved into a dictionary
for file in sorted(os.listdir(indir)):
    if file.endswith('.gbk'):
        new_dict = parse_GenBank(f'{indir}/{file}', newfile)
        newfile.update(new_dict)

#Individual records are saved in separate files                
for locus_name in newfile.keys():
    with open(f'{outdir}/{locus_name}.gbk', 'w') as outgbk:
        outgbk.write(newfile[locus_name])
