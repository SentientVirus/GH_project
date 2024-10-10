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


pos_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab/all_subsets_nogap'
out_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/fna/all_subsets_nogap'


with open(f'{pos_dir}/all_subsets_positions_aln.tab', 'w+') as gpos:
    gpos.write('locus_tag\tstrain\tgene_name\tstart\tend\tstrand\n')
    

# =============================================================================
# Section to make files with the gene positions in the alignment
# =============================================================================
group = 'all_subsets'
outtab = f'{pos_dir}/{group}_positions_aln.tab'
with open(f'{out_dir}/{group}_seqs.mafft.fasta') as in_aln, open(f'{pos_dir}/{group}_positions.tab') as in_tab:
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
                            
                        
                    
        
            
        


# with open(f'{out_dir}/{group}_seqs.mafft.fasta', 'a+') as last_outfile:
#     SeqIO.write(new_record, last_outfile, 'fasta')