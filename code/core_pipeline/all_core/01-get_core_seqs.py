#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:09:59 2024

@author: marina
"""

import os, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

infile = os.path.expanduser('~') + '/GH_project/all_core/SCPOs.txt'
outdir = os.path.expanduser('~') + '/GH_project/all_core/fasta'
inpath = os.path.expanduser('~') + '/Akunkeei_files/cds'
inpath2 = os.path.expanduser('~') + '/Akunkeei_files/faa'
MP2_cds = os.path.expanduser('~') + '/LAB/annotations_220613/MP2.gbk'

MP2_records = {} #record for record in SeqIO.parse(MP2_faa, 'fasta')]
check = False
for record in SeqIO.parse(MP2_cds, 'genbank'):
    if check == False:
        full_seq = record.seq
        check = True
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers.keys():
                loctag = feature.qualifiers['locus_tag'][0]
                MP2_records[loctag] = feature.location
                print(feature.location.start, feature.location.strand)

MP2_records_fna = []
for locus in MP2_records.keys():
    position = MP2_records[locus]
    start = position.start
    end = position.end
    strand = position.strand
    gene_seq = full_seq[start:end]
    if strand == -1:
        gene_seq = gene_seq.reverse_complement()
    record = SeqRecord(gene_seq, id = locus, name = locus, description = locus)
    MP2_records_fna.append(record)

if not os.path.exists(outdir):
    os.makedirs(outdir)

SCPO_dict = {}
with open(infile) as SCPOs:
    for SCPO in SCPOs:
        scpo_ids = SCPO.split(' ')
        scpo_list = [scpon.split('|')[1] for scpon in scpo_ids]
        pattern = '[AG]\d{4}'
        scpo_filtered = [scpo for scpo in scpo_list if scpo.startswith('H') or re.search(pattern, scpo) != None or scpo.split('_')[0] in ['LDX55', 'K2W83', 'MP2', 'fhon2']]
        scpo_name = scpo_filtered[6]
        SCPO_dict[scpo_filtered[0]] = scpo_filtered
        outfile = f'{outdir}/{scpo_name}.fna'
        outfile2 = outfile.replace('.fna', '.faa')
        with open(outfile, 'w') as out_handle:
            out_handle.write('')
            
        record_list = []
        prot_list = []
        
        for scpo in scpo_filtered:
            strain = scpo.split('_')[0]
            # if strain.startswith('H'):
            #     strain = strain[:4] + '-' + strain[4:]
            if strain == 'K2W83':
                strain = 'DSMZ12361'
            elif strain == 'LDX55':
                strain = 'IBH001'
            elif strain == 'fhon2':
                strain = 'Fhon2'
            if strain != 'MP2':
                infile = f'{inpath}/{strain}_cds_from_genomic.fna'
                infile2 = f'{inpath2}/{strain}_protein.faa'
                    
                for record in SeqIO.parse(infile, format = 'fasta'):
                    locus_tag = record.description.split('locus_tag=')[1].split(']')[0]
                    if locus_tag.startswith('AKU'):
                        locus_tag = locus_tag[3:]
                    if locus_tag == scpo.replace('-', '').upper():
                        print(locus_tag, scpo)
                        record.id = locus_tag
                        record_list.append(record)

            else:
                for record in MP2_records_fna:
                    if record.id == scpo:
                        print(record.id, scpo)
                        record_list.append(record)

        SeqIO.write(record_list, outfile, 'fasta')
        
        for record in record_list:
            record.seq = record.seq.translate(stop_symbol = '')
        SeqIO.write(record_list, outfile2, 'fasta')
        
#Here I loop through the cds files with all faa sequences of genes and retrieve
#the SCPOs and add them to the corresponding files one by one, since they should be in alphabetical order
#I should use strain A0901 as the key in the beginning, and update the name of the files later if I want