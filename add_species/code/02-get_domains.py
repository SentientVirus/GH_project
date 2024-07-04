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

in_file = 'results/interproscan/GH70_species.tsv'
GH70_text = 'Glycosyl hydrolase family 70'
fasta_file = 'files/outgroups_ingroups.fasta'
outfile = 'results/domains/outgroup_domains.fasta'
outdir = os.path.dirname(outfile)

if not os.path.exists(outdir):
    os.makedirs(outdir)

domain_dict = {}
with open(in_file) as interpro:
    domains = pd.read_csv(interpro, sep = '\t', header = None)
    for index, line in domains.iterrows():
        if line[5] == GH70_text:
            if line[0] not in domain_dict.keys():
                domain_dict[line[0]] = (line[6], line[7])
            else:
                domain_dict[line[0]+'_2'] = (line[6], line[7])
                
my_prots = []                
with open(fasta_file) as prots:
    for prot in SeqIO.parse(prots, 'fasta'):
        new_prot = seqr(prot.seq, id = prot.id, name = prot.name, description = prot.description)
        pos = domain_dict[prot.id]
        new_prot.seq = prot.seq[pos[0]:pos[1]]
        my_prots.append(new_prot)
        if prot.id + '_2' in domain_dict.keys():
            new_prot2 = seqr(prot.seq, id = prot.id + '_2', name = prot.name, description = prot.description + '_2')
            pos2 =  domain_dict[new_prot2.id]
            new_prot2.seq = prot.seq[pos2[0]:pos2[1]]
            my_prots.append(new_prot2)
            
SeqIO.write(my_prots, outfile, 'fasta')
