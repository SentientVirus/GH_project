#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 17:38:21 2023

@author: marina
"""
import os
from Bio import SeqIO
from Bio import pairwise2

# =============================================================================
# Script to add the IBH001 core genes to compute the mean core dS
# =============================================================================

indir = 'core_genes'
tag_dict = {}
cds_path = '../Akunkeei_files/cds' #'../Akunkeei_files/cds/IBH001_cds_from_genomic.fna'
outfile = 'IBH001_tags.tsv'
query = 'IBH001'
reference = 'H4B2-06J'
count1 = 0
count2 = 0

ref_dict = {}
locus_list = []
H4B206J_dict = {}
percent_threshold = 0.48

# =============================================================================
# 1. Retrieve locus tags of core genes in A1401
# =============================================================================
for file in os.listdir(indir):
    if file.endswith('genes.fasta'):
        locus_tag = file.split('_genes')[0]
        #locus_list.append(locus_tag)
        with open(f'{indir}/{file}') as cds:
            fasta_seqs = SeqIO.parse(cds, format = 'fasta')
            key = [f.id for f in fasta_seqs if 'H4B2-06J' in f.id][0].replace('-', '')
            H4B206J_dict[key] = locus_tag
       
with open(f'{cds_path}/{reference}_cds_from_genomic.fna') as strain:
    fasta_seqs = SeqIO.parse(strain, format = 'fasta')
    for f in fasta_seqs:
        loctag = f.description.split('locus_tag=AKU')[1].split(']')[0].replace('UNKNOWN', '')
        if loctag in H4B206J_dict.keys():
            ref_dict[H4B206J_dict[loctag]] = f.seq
        

# =============================================================================
# 2. Get those locus tags in IBH001 by adding them to a dictionary and making
# sure that I get the most similar sequence (a couple of them seem to be split)
# =============================================================================
j = 0
for tag in sorted(ref_dict.keys()):
    tag_dict[tag] = []
    with open(f'{cds_path}/{query}_cds_from_genomic.fna') as cds_list:
        IBH_seqs = SeqIO.parse(cds_list, format = 'fasta')
        sorted_IBH = sorted(list(IBH_seqs), key = lambda a: a.description.split('locus_tag=')[1].split(']')[0])
        for i in range(j, len(sorted_IBH)):
            aln = pairwise2.align.globalxx(sorted_IBH[i].seq.translate(to_stop = True), ref_dict[tag].translate(to_stop = True))[0]
            matches = [nt_ref == nt_IBH for nt_ref, nt_IBH in zip(aln[0], aln[1])]
            similarity = sum(matches)/len(matches)
            if similarity > percent_threshold and tag_dict[tag] == []:
                loc_tag = sorted_IBH[i].description.split('locus_tag=')[1].split(']')[0]
                tag_dict[tag] = (loc_tag, similarity)
                j = i + 1
                count1 += 1
                if similarity > 0.9:
                    break
            elif similarity > percent_threshold and similarity > tag_dict[tag][1]:
                loc_tag = sorted_IBH[i].description.split('locus_tag=')[1].split(']')[0]
                tag_dict[tag] = (loc_tag, similarity)
                j = i + 1
                count2 += 1
                if similarity > 0.9:
                    break
            #     count1 += 1
            #     break
    # if len(tag_dict[tag]) == 0:
    #     count1 += 1
    # elif len(tag_dict[tag]) > 2:
    #     count2 += 1

missing_genes = count1/len(tag_dict)*100
extra_genes = count2/len(tag_dict)*100

#print(f'{missing_genes:.2f}% of the core genes are missing.')
print(f'{missing_genes:.2f}% of the core genes are missing and {extra_genes:.2f}% are duplicated.')

# Sequences are aligned AND trimmed in Karl's files, maybe I should try to retrieve the original genes in A1401 and save them to a dictionary
# Maybe try with H4B2-06J?

# =============================================================================
# 3. Create output file from dictionary
# =============================================================================

with open(outfile, 'w') as outf:
    outf.write('reference_tag\tIBH001_tag\n')
    {outf.write(f'{key}\t{tag_dict[key][0]}\n') for key in tag_dict.keys()}

