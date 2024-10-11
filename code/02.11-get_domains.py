#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 19:34:41 2024

@author: marina
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord as seqr
import os

# =============================================================================
# 0. Read inputs and output from Snakemake and set domain to retrieve
# =============================================================================
in_file = snakemake.input['annot'] 
GH70_text = 'Glycosyl hydrolase family 70'
fasta_file = snakemake.input['seq'] 
outfile = snakemake.output[0]

# =============================================================================
# 1. Get the positions of the domains in the proteins and save to dictionaru
# =============================================================================
domain_dict = {}
with open(in_file) as interpro:
    domains = pd.read_csv(interpro, sep = '\t', header = None)
    for index, line in domains.iterrows():
        if line[5] == GH70_text:
            if line[0] not in domain_dict.keys():
                domain_dict[line[0]] = (line[6], line[7])
            else:
                domain_dict[line[0]+'_2'] = (line[6], line[7])
                
                
# =============================================================================
# 2. Retrieve the domains from the complete protein sequences
# =============================================================================
my_prots = []                
with open(fasta_file) as prots:
    for prot in SeqIO.parse(prots, 'fasta'):
        new_prot = seqr(prot.seq, id = prot.id, name = prot.name, description = prot.description)
        pos = domain_dict[prot.id]
        new_prot.seq = prot.seq[pos[0]:pos[1]+1]
        my_prots.append(new_prot)

        print(f'{new_prot.id} will be added to {outfile}')

# =============================================================================
# 3. Save the domains to a file
# =============================================================================
with open(outfile, 'w') as handle:
    SeqIO.write(my_prots, handle, 'fasta')
