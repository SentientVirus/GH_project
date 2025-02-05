#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 15:43:05 2024

@author: marina
"""
import logging, traceback, sys
from Bio import SeqIO
import os
import subprocess
import pandas as pd

prefixes = ['E', 'F']

for prefix in prefixes:
    pos_dir = os.path.expanduser('~') + f'/GH_project/RDP5_analysis/files/tab/{prefix}'
    out_dir = os.path.expanduser('~') + f'/GH_project/RDP5_analysis/files/fna/{prefix}'
    intab = f'{pos_dir}/{prefix}.tab'
    
    
    with open(f'{pos_dir}/{prefix}_aln.tab', 'w+') as gpos:
        gpos.write('locus_tag\tstrain\tgene_name\tstart\tend\tstrand\n')
        
    
    # =============================================================================
    # Section to make files with the gene positions in the alignment
    # =============================================================================
    group = prefix
    outtab = f'{pos_dir}/{group}_aln.tab'
    with open(f'{out_dir}/{group}.mafft.fasta') as in_aln, open(intab) as in_tab:
        tab_df = pd.read_csv(in_tab, sep = '\t')
        aln_reader = SeqIO.parse(in_aln, 'fasta')
        new_df = tab_df.copy()
        for record in aln_reader:
            print(record.id)
            # nogap = 0
            # count = 0
            for index, row in tab_df.iterrows():
                nogap = 0
                count = 0
                if row['strain'] in record.id:
                    print(row['locus_tag'], record.id, group)
                    print(row['gene_name'], row['start'], row['end'])
                    for nt in record.seq:
                        count += 1
                        if nt != '-':
                            nogap += 1
                        if nogap == row['start']:
                            new_df.loc[index, 'start'] = count
                        elif nogap == row['end']:
                            new_df.loc[index, 'end'] = count
                            print(row['gene_name'], new_df.loc[index, 'start'], new_df.loc[index, 'end'])
                            print(index, count, nogap)
                            # index += 1
                            # break
                            
    new_df.to_csv(outtab, sep = '\t', index = False)