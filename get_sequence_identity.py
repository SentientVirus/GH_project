#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:20:09 2022

@author: marina
"""

# =============================================================================
# This script reads a fasta file with nucleotide information and calculates the
# overall percentage of identity between every pair of sequences in the file.
# =============================================================================

import os, re
from Bio import SeqIO
# from Bio.Seq import Seq
# from statistics import mean
from Bio.Emboss.Applications import NeedleCommandline

pathname = 'data/fasta'
GH_type = ['GH70', 'GH32']
out_path = 'tables'
out_file = 'percentage_identity_subset.tab'

if not os.path.exists(out_path):
    os.makedirs(out_path)

def needle_align_code(query_seq, target_seq):
    needle_cline = NeedleCommandline(asequence="asis:" + query_seq,
                                     bsequence="asis:" + target_seq,
                                     aformat="simple",
                                     gapopen=10,
                                     gapextend=0.5,
                                     outfile='stdout'
                                     )
    out_data, err = needle_cline()
    out_split = out_data.split("\n")
    p = re.compile("\((.*)\)")
    return p.search(out_split[25]).group(1).replace("%", "")


with open(f'{out_path}/{out_file}', 'w') as out:
    out.write('locus1\tlocus2\tGH_type\t%identity\n')


pos_dict = {}
for GH in GH_type:
    workdir = f'{pathname}/{GH}'
    for file in os.listdir(workdir):
        if file.endswith('subset.fna') and GH not in file and 'complete' not in file:
            GH_type = file.split('_')[0]
            if GH_type == 'NCB':
                GH_type = 'NGB'
            with open(f'{workdir}/{file}') as handle:
                record_list = []
                for record in SeqIO.parse(handle, 'fasta'):
                    record_list.append(record)
            for i in range(len(record_list)-1):
                for j in range(i+1, len(record_list)):
                    if GH_type != 'GS2' or (len(record_list[i].seq) > 2000 and len(record_list[j].seq) > 2000):
                        perc_ident = needle_align_code(record_list[i].seq, record_list[j].seq)
                        pos_dict[(record_list[i].id, record_list[j].id)] = perc_ident
                        with open(f'{out_path}/{out_file}', 'a') as out:
                            out.write(f'{record_list[i].id}\t{record_list[j].id}\t{GH_type}\t{perc_ident}\n')
        