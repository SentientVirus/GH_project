# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:39:34 2021

@author: usuario
"""
import re, os
from Bio import SeqIO

GH_types = ['S1', 'S2a', 'S2b', 'S3']
family = 'GH32'

class motif:
    def __init__(self, num, seq, start, end): #Define motif class
        self.num = num #Motif number (1-7)
        self.seq = seq #Motif sequence
        self.start = start #Start position of the motif in a protein sequence
        self.end = end #End position of the motif in the sequence
    def __str__(self): #What print(class) returns
        return f'<Motif {self.num} from position {self.start} to {self.end}>'

class prot_sequence: #All motifs associated to protein (1 sequence/position per motif)
    def __init__(self, locus, I, II, III, IV, V, VI, VII, VIII):
        self.locus = locus #Locus tag
        self.I = I #Domain 1
        self.II = II #Domain 2
        self.III = III
        self.IV = IV
        self.V = V
        self.VI = VI
        self.VII = VII #Domain 7
        self.VIII = VIII #Domain 8
        self.n_motifs = 0 #Number of motifs
        self.bool_list = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    def calculate_n_motifs(self): #Function to calculate number of motifs in the protein
        motif_list = [self.I, self.II, self.III, self.IV, self.V, self.VI, self.VII, self.VIII]
        bool_list = [0 + len(motif) for motif in motif_list]
        self.bool_list = [0] + bool_list
        self.n_motifs = sum(bool_list)
    def __str__(self): #What print(class) returns
        return f'<{self.locus} with {self.n_motifs} motifs>'
    
def search_motif(seq_str, num): #Function to search for any motif
    motif_list = ['', r'[FY]SND[IV]QS', r'WY[HG][AVL][TIK]T[KT]DF',
                  r'A[TS]G[SAT][VI][LVIF](.{1,1})N',
                  r'FR(.{0,1})[VIP][VYQ][VIQ](.{2,2})',
                  r'ETP[ILN][VI][KRQ]T', r'D[DAF](.{2,3})D(.{1,1})YA', 
                  r'[NG]W[SL]G[NDS]W[LIV]Y[STW]', r'DN[SAT][TS][LIV]E']
    
    motifs = []

    motif_info = re.finditer(motif_list[num], seq_str)
    for motif_match in motif_info:
        if motif_match is not None: #If there is a match, store new values
            motif_seq = motif_match.group(0)
            start = motif_match.start() + 1
            end = motif_match.end() + 1
        else: #Otherwise use generic values
            motif_seq = ''
            start = 0
            end = 0
        motifs.append([num, motif_seq, start, end])
    return motifs

infolder = os.path.expanduser('~') + '/GH_project/data/fasta/GH32'
workdir = os.path.expanduser('~') + '/GH_project/motifs'
output_file = f'{workdir}/GH32_motif_presence_pos.tsv'
output_file2 = f'{workdir}/GH32_motifs_pos.tsv'
prot_list = []

with open(output_file, 'w') as out_file, open(output_file2, 'w') as out_file2:
    out_file.write('locus\ttype\tmotif_I\tmotif_II\tmotif_III\tmotif_IV\tmotif_V\tmotif_VI\tmotif_VII\tmotif_VIII\n')
    out_file2.write('locus\ttype\tmotif_type\tsequence\tstart\tend\ttotal_length\n')

tag_type = {}
for GH_type in GH_types:
    with open(f'{infolder}/{GH_type}_repset.fna') as handle:
        aa_seqs = SeqIO.parse(handle, 'fasta')
        for fasta in aa_seqs:
            tag_type[fasta.id] = GH_type

with open(f'{infolder}/complete_GH32_repset.fna') as handle:
    codon_seqs = SeqIO.parse(handle, 'fasta')
    for fasta in codon_seqs:
        name, sequence = fasta.id, str(fasta.seq.translate())
        if 'K2W83' in name:
            name = name.replace('K2W83_RS', 'DSM_')
        elif 'LDX55' in name:
            name = name.replace('LDX55', 'IBH001')
        elif 'APS55' in name:
            name = name.replace('APS55_RS', 'MP2_')
        category = tag_type[fasta.id]
        motif_list = ['']
        for i in range(1, 9):
            motif_matches = search_motif(sequence, i)
            motif_list.append(motif_matches)
            for mot in motif_matches:
                mot = motif(*mot)
                if len(mot.seq) > 0:
                    with open(output_file2, 'a') as out_file2:
                        out_file2.write(f'{name}\t{category}\t{mot.num}\t{mot.seq}\t{mot.start}\t{mot.end}\t{len(sequence)}\n')
        prot = prot_sequence(name, motif_list[1], motif_list[2], motif_list[3], motif_list[4], motif_list[5], motif_list[6], motif_list[7], motif_list[8])
        prot.calculate_n_motifs()
        prot_list.append(prot)
        bool_list = list(map(int, prot.bool_list))
        with open(output_file, 'a') as out_file:
            out_file.write(f'{prot.locus}\t{category}\t{bool_list[1]}\t{bool_list[2]}\t'+
                           f'{bool_list[3]}\t{bool_list[4]}\t{bool_list[5]}\t'+
                           f'{bool_list[6]}\t{bool_list[7]}\n')
        
