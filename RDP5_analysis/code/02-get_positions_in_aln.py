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

# =============================================================================
# 0. Logging
# =============================================================================
log = os.path.expanduser('~') + '/GH_project/RDP5_analysis/logs/get_A-D_alnpos.log'

if os.path.exists(log):
    os.remove(log)

if not os.path.exists(os.path.dirname(log)):
    os.makedirs(os.path.dirname(log))

logger = logging.getLogger()

logging.basicConfig(filename = log, level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(log, 'a')

# =============================================================================
# 1. Define inputs and paths
# =============================================================================

pos_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab'
out_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/fna'

for prefix in ['A', 'B', 'C', 'D']:
    
    outtab = f'{pos_dir}/{prefix}/{prefix}_aln.tab'
    intab = outtab.replace('_aln', '')
    inaln = f'{out_dir}/{prefix}/{prefix}.mafft.fasta'
    
    
    with open(outtab, 'w+') as gpos:
        gpos.write('locus_tag\tstrain\tgene_name\tstart\tend\tstrand\n')
    
    # =============================================================================
    # 2. Section to make files with the gene positions in the alignment
    # =============================================================================
    
    with open(inaln) as in_aln, open(intab) as in_tab:
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
                    print(row['locus_tag'], record.id)
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