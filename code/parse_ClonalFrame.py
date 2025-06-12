#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 15:46:33 2025

Script to parse ClonalFrameML results and count the total lengths of the 
genomes in the dataset.

@author: Marina Mota-Merlo
"""

import time
start_time = time.time()

import os

workdir = os.path.expanduser('~') + '/GH_project'
phylogroups = ['A', 'B-C']

def get_lengths(infile):
    length_dict = {}
    with open(infile) as handle:
        file_reader = [line for line in handle if not line.startswith('#') and not line.startswith('=')]
        for l in file_reader:
            if l[0] == '>':
                strain = l.split()[-1].replace('.fna', '')
            if strain not in length_dict:
                length_dict[strain] = 0
            elif l[0] != '>':
                length_dict[strain] += len(l.replace('-', '').strip())
    return length_dict

for i in phylogroups:
    infile = f'{workdir}/xmfa_sslcb500_genome_len/phylogroup{i}/AK38gh_phylogr{i}_sslcb500_mod.fna.xmfa'
    outfile = infile.replace('_sslcb500_mod.fna.xmfa', '_lengths.tab')

    length_dict = get_lengths(infile)
    with open(outfile, 'w') as handle:
        handle.write('Strain\tGenome_length(bp)\n')
        [handle.write(f'{strain.replace("DSM", "DSMZ")}\t{length_dict[strain]}\n') for strain in sorted(length_dict.keys())]
    
print("--- %s seconds ---" % (time.time() - start_time))