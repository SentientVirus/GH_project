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

# OBS! Modify this to retrieve only one domain from the GTFC proteins, merging the two domains that are predicted in Interproscan
in_file = os.path.expanduser('~') + '/GH_project/add_species/results/interproscan/GH70_species.tsv'
GH70_text = 'Glycosyl hydrolase family 70'
fasta_file = os.path.expanduser('~') + '/GH_project/add_species/files/ingroups.faa'
outfile = os.path.expanduser('~') + '/GH_project/add_species/results/domains/ingroup_domains.faa'
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
            elif 'P_kribbensis' in line[0] or 'E_undae' in line[0] or 'Exiguobacterium' in line[0]:
                domain_dict[line[0]] = [domain_dict[line[0]], (line[6], line[7])] #(max(domain_dict[line[0]][0], line[6]), 
                                        #max(domain_dict[line[0]][1], line[7]))#[domain_dict[line[0]], (line[6], line[7])]
            else:
                domain_dict[line[0]+'_2'] = (line[6], line[7])
                
my_prots = []                
with open(fasta_file) as prots:
    for prot in SeqIO.parse(prots, 'fasta'):
        new_prot = seqr(prot.seq, id = prot.id, name = prot.name, description = prot.description)
        if type(domain_dict[prot.id]) != list:
            pos = domain_dict[prot.id]
            new_prot.seq = prot.seq[pos[0]:pos[1]+1]
        else:
            pos1 = (min(domain_dict[prot.id][0][0], domain_dict[prot.id][1][0]),
                    min(domain_dict[prot.id][0][1], domain_dict[prot.id][1][1]))
            pos2 = (max(domain_dict[prot.id][0][0], domain_dict[prot.id][1][0]),
                    max(domain_dict[prot.id][0][1], domain_dict[prot.id][1][1]))
            new_prot.seq = prot.seq[pos2[0]:pos2[1]+1]+prot.seq[pos1[0]:pos1[1]+1]
        my_prots.append(new_prot)
        print(f'{new_prot.id} will be added to {outfile}')
        if prot.id + '_2' in domain_dict.keys():
            new_prot2 = seqr(prot.seq, id = prot.id + '_2', name = prot.name, description = prot.description + '_2')
            pos2 =  domain_dict[new_prot2.id]
            new_prot2.seq = prot.seq[pos2[0]:pos2[1]+1]
            my_prots.append(new_prot2)
            print(f'{new_prot2.id} will be added to {outfile}')
            
SeqIO.write(my_prots, outfile, 'fasta')
