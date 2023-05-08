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

indir = '/'.join(snakemake.input[0].split('/')[:-1]) #Directory of input files
outdir = '/'.join(snakemake.output[0].split('/')[:-1]) #Directory of outputs
in_suffix = snakemake.input[0].split('.')[-1] #Duffix of inputs (gbff)
out_suffix = snakemake.output[0].split('.')[-1] #Suffix of outputs (gbk)

new_dict = {}
newfile = {}

#Function to parse GenBank files, it outputs record names and record info
def parse_GenBank(file, locus_dict):
    locus_dict = {}
    locus_2_strain = {}
    with open(file, 'r') as handle:
        n = 1
        f = (line for line in handle)
        # check = False
        for line in f:
            if 'LOCUS' in line:
                locus_name = line.split()[1]
                locus_dict[locus_name] = ''
            elif 'DEFINITION' in line:
                locus_2_strain[locus_name] = f'{line.split()[4]}_{n}'
                if '12361' in line:
                    locus_2_strain[locus_name] = f'DSMZ12361_{n}'
                n += 1
                print(locus_2_strain[locus_name])
            # elif '_1' not in locus_name and 'ORIGIN' in line:
            #     check = True
            # elif '//' in line:
            #     check = False
            # if check and 'ORIGIN' not in line:
            #     line = '  ' + line
            locus_dict[locus_name] += line
        final_locus_dict = {locus_2_strain[locname]: locus_dict[locname] for locname in locus_dict.keys()}
    return final_locus_dict

#Results are saved into a dictionary
for file in sorted(os.listdir(indir)):
    if file.endswith(in_suffix):
        new_dict = parse_GenBank(f'{indir}/{file}', newfile)
        newfile.update(new_dict)

#Individual records are saved in separate files                
for locus_name in newfile.keys():
    with open(f'{outdir}/{locus_name}.{out_suffix}', 'w') as outgbk:
        outgbk.write(newfile[locus_name])
