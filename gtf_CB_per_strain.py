# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import os, logging, traceback
import csv
from Bio import SeqIO
output_file = snakemake.output[0]
outdir = output_file.split('/')[:-1]
outdir = '/'.join(outdir)
if not os.path.exists(outdir):
    os.makedirs(outdir)

#Group on the upper part of the phylogeny
GS1 = snakemake.params.GS1
GS2 = snakemake.params.GS2
BRS = snakemake.params.BRS
short = snakemake.params.short
NGB = snakemake.params.NGB
S1 = snakemake.params.S1
S2a = snakemake.params.S2a
S2b = snakemake.params.S2b
S2 = S2a + S2b
S3 = snakemake.params.S3
phylo = snakemake.input.tree
GH70 = snakemake.input.GH70

value_dict = {}
with open(GH70) as lengths:
    for record in SeqIO.parse(lengths, 'fasta'):
        # if record.id[0] == 'H':
        #     record.id = record.id[:4] + '-' + record.id[4:]
        if len(record.seq) >= 750: #Add only 0.5 to the count for proteins that are not full-length
            value_dict[record.id] = 1
#            if 'G0102' in record.id:
#                print(f'{record.id}: {len(record.seq)}')
        elif record.id not in short:
            value_dict[record.id] = 0.5


all_GH70 = GS1 + GS2 + BRS + short + NGB

all_GH32 = S1 + S2a + S2b + S3

count = [0, 0, 0, 0, 0, 0, 0, 0]
with open(output_file, 'w') as tab:
    tab.write('strain\tGS1\tGS2\tBRS\tNGB\tshort_gene\ttotal_GH70\tS1\tS2\tS3\ttotal_GH32\tloci\n')

with open(phylo) as phy: #Order by phylogeny
    p = (tip for tip in phy)
    for strain in p:
        loci = ''
        strain = strain.strip()
        #strains.append(strain)
        if len(strain) > 0:
            strain = strain.replace('-', '')
            count[0] = sum(value_dict[s] for s in GS1 if strain in s)
            count[1] = sum(value_dict[s] for s in GS2 if strain in s)
            count[2] = sum(value_dict[s] for s in BRS if strain in s)
            count[3] = sum([1 for s in NGB if strain in s])
            count[4] = sum([1 for s in short if strain in s])
            count = [int(num) for num in count] #Convert floats to ints
            total = sum(count)
            for GH70 in all_GH70: #Add loci names
                if strain in GH70:
                    loci += f'{GH70}, '
            count[5] = sum([1 for s in S1 if strain in s])
            count[6] = sum([1 for s in S2 if strain in s])
            count[7] = sum([1 for s in S3 if strain in s])
            total2 = sum(count[5:8])
            #Add also cell wall-binding domains
            '''with open(f'../Alignment/interproscan/{strain}.tsv') as interpro:
                annot = csv.reader(interpro, delimiter ='\t')
                cb = 0
                cwb = 0
                for row in annot:
                    if row[3] == 'Pfam':
                        if 'Choline-binding' in row[5]:
                            cb += 1
                        elif 'cell wall binding' in row[5]:
                            cwb += 1
            count[5] = cb
            count[6] = cwb
            count[7] = cb + cwb'''
            #Write info to file
            with open(output_file, 'a') as tab:
                tab.write(f'{strain}\t{count[0]}\t{count[1]}\t{count[2]}\t'+
                          f'{count[3]}\t{count[4]}\t{total}\t{count[5]}\t'+
                          f'{count[6]}\t{count[7]}\t{total2}\t{loci[:-2]}\n')
            #strain = next(p)
        count = [0, 0, 0, 0, 0, 0, 0, 0]
