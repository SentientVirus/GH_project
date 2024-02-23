#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:20:09 2022

This script reads a fasta file with nucleotide information and calculates the
overall percentage of identity between every pair of sequences in the file.

@author: Marina Mota Merlo

"""


import os, re
from Bio import SeqIO
# from Bio.Seq import Seq
# from statistics import mean
from Bio.Emboss.Applications import NeedleCommandline

home = os.path.expanduser('~')
pathname = f'{home}/GH_project/data/fasta'
GH_type = ['GH70', 'GH32']
out_path = f'{home}/GH_project/tables'
suffixes = ['subset', 'repset', 'all']

if not os.path.exists(out_path):
    os.makedirs(out_path)

def needle_align_code(query_seq, target_seq):
    needle_cline = NeedleCommandline(asequence="asis:" + query_seq,
                                     bsequence="asis:" + target_seq,
                                     sprotein=True,
                                     aformat="simple",
                                     gapopen=10,
                                     gapextend=0.5,
                                     outfile='stdout'
                                     )
    out_data, err = needle_cline()
    out_split = out_data.split("\n")
    p = re.compile("\((.*)\)")
    return p.search(out_split[25]).group(1).replace("%", "")

def GH_type_dict(GH_file, GH_dict):
    GH_type = os.path.basename(file).split('_')[0]
    with open(GH_file) as GH_list:
        for record in SeqIO.parse(GH_list, 'fasta'):
            GH_dict[record.id] = GH_type
            
GH_dict = {}
GH70_subtypes = ['GS1', 'GS2', 'BRS', 'NGB', 'short']
GH32_subtypes = ['S1', 'S2a', 'S2b', 'S3']
for subtype in GH70_subtypes:
    file = f'{pathname}/GH70/{subtype}_all.faa'
    GH_type_dict(file, GH_dict)
for subtype in GH32_subtypes: 
    file = f'{pathname}/GH32/{subtype}_all.faa'
    GH_type_dict(file, GH_dict)

for suffix in suffixes:
    out_file = f'percentage_identity_{suffix}.tab'
    with open(f'{out_path}/{out_file}', 'w') as out:
        out.write('locus1\tlocus2\tGH_type\t%identity\n')
    
    
    pos_dict = {}
    for GH in GH_type:
        workdir = f'{pathname}/{GH}'
        for file in os.listdir(workdir):
            if file.endswith(f'{suffix}.fna') and 'complete' not in file and ('GH32' in file or 'short' in file or 'GH70_functional' in file):
                print(file)
                # GH_name = file.split('_')[0]
                # print(GH_name)
                # if GH_name == 'NCB':
                #     GH_name = 'NGB'
                with open(f'{workdir}/{file}') as handle:
                    record_list = []
                    for record in SeqIO.parse(handle, 'fasta'):
                        record_list.append(record)
                for i in range(len(record_list)-1):
                    for j in range(i+1, len(record_list)):
                        # if GH_name != 'GS2' or (len(record_list[i].seq) > 2000 and len(record_list[j].seq) > 2000):
                        GH_name = f'{GH_dict[record_list[i].id]}/{GH_dict[record_list[j].id]}'
                        perc_ident = needle_align_code(record_list[i].seq, record_list[j].seq)
                        pos_dict[(record_list[i].id, record_list[j].id)] = perc_ident
                        with open(f'{out_path}/{out_file}', 'a') as out:
                            out.write(f'{record_list[i].id}\t{record_list[j].id}\t{GH_name}\t{perc_ident}\n')
        
