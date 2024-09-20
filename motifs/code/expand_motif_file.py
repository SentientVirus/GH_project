# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 14:56:09 2021

@author: usuario
"""
import csv

short = ['A1202_14500', 'A1404_14290', 'A1805_13670', 'H4B2-04J_14250', 
         'H2B1-05J_13740', 'H4B1-01A_13720','A0901_14250', 'A1001_13210', 
         'G0403_13980', 'fhon2_14510', 'A1003_13400', 'H1B1-04J_13940', 
         'H1B3-02M_13810', 'H4B4-02J_13530', 'H4B2-02J_13700', 'H3B1-01A_14180', 
         'H3B1-02X_14140', 'H3B1-01J_13900', 'H3B1-03M_14110', 'H3B1-03M_14110']

locus_dict = {}
presence_dict = {}
filename = 'GH70_motifs.tsv' #'gtf_motifs.tsv'
filename2 = 'GH70_motif_presence.tsv' #'motif_presence.tsv'
with open(filename, 'r') as full_motifs, open(filename2, 'r') as present_motifs:
    motifs = csv.reader(full_motifs, delimiter="\t")
    presence = csv.reader(present_motifs, delimiter="\t")
    n = 0
    for row in motifs:
        locus_dict[f'{row[0]}_{n}'] = [row[0], row[1], row[2], row[3], row[4], row[5], row[6]]
        n += 1
    for locus in presence:
        presence_dict[locus[0]] = [locus[1], locus[2], locus[3], locus[4],
                                   locus[5], locus[6], locus[7], locus[8]]
second_dict = {}
same_seqs = '../GH_sequences/same_domain_sequences.txt' #'same_protein_sequences.txt'        
with open(same_seqs) as same:
    gen = (line for line in same)
    for line in gen:
        if line[0] in 'AGHMf': #First letters of strain names
            words = line.replace(' \n', '')
            words = words.split(' ')
            domain = words[0].split('_')
            if len(domain) > 2:
                domain = ('_').join(domain[:-1])
            else:
                domain = words[0]
            if domain not in short:
                for i in range(1, len(words)):
                    keys = list(locus_dict.keys())
                    for key in keys:
                        if domain == key:
                            #print(words[0])
                            num = keys.index(key)
                            second_dict[f'{words[i]}_{num}'] = [words[i]] + locus_dict[f'{domain}_{num}'][1:]
                    presence_dict[words[i]] = presence_dict[domain]
                
locus_dict.update(second_dict)

outfile = 'all_GH70_motifs.tsv' #'all_gtf_motifs.tsv'
outfile2 = 'all_GH70_motif_presence.tsv' #'all_motif_presence.tsv'
with open(outfile, 'w') as out:
    [out.write(f'{el}\t') for el in locus_dict['locus_0']]
    out.write('\n')
    keys = list(locus_dict.keys())
    keys = keys[1:]
    for key in sorted(keys):
        [out.write(f'{el}\t') for el in locus_dict[key]]
        out.write('\n')
        
with open(outfile2, 'w') as out2:
    key = 'locus'
    out2.write(f'{key}\t')
    [out2.write(f'{el}\t') for el in presence_dict[key]]
    out2.write('\n')
    presence_dict.pop(key)
    for key in sorted(presence_dict.keys()):
        out2.write(f'{key}\t')
        [out2.write(f'{el}\t') for el in presence_dict[key]]
        out2.write('\n')
        
                
