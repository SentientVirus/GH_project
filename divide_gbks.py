#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 10:01:24 2022

@author: marina
"""
import os

indir = 'data/genbanks'
outdir = 'gbks'
newfile = {}

for file in sorted(os.listdir(indir)):
    if file.endswith('.gbk'):
        with open(f'{indir}/{file}', 'r') as handle:
            f = (line for line in handle)
            # for record in GenBank.parse(handle):
            #     print(record)
            for line in f:
                if 'LOCUS' in line:
                    locus_name = line.split()[1]
                    newfile[locus_name] = ''
                    print(locus_name)
                newfile[locus_name] += line
                
for locus_name in newfile.keys():
    with open(f'{outdir}/{locus_name}', 'w') as outgbk:
        outgbk.write(newfile[locus_name])
