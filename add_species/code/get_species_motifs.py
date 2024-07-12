# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:39:34 2021
Script to tetrieve the motifs in sequences from other species.
Modified on the 12th of July 2024 from the code to retrieve the motifs from 
GH70 genes. In the original script, an environment with Mafft and Pal2nal was 
used first to get the codon alignments of the sequences. In this case, such an
alignment is not needed as input.

@author: Marina Mota-Merlo
"""
import re
import os
from Bio import SeqIO

GH_types = ['GS1', 'GS2', 'BRS', 'NGB']

# =============================================================================
# 0. Object class definition
# =============================================================================
class motif:
    def __init__(self, num, seq, start, end): #Define motif class
        self.num = num #Motif number (1-7)
        self.seq = seq #Motif sequence
        self.start = start #Start position of the motif in a protein sequence
        self.end = end #End position of the motif in the sequence
    def __str__(self): #What print(class) returns
        return f'<Motif {self.num} from position {self.start} to {self.end}>'

class prot_sequence: #All motifs associated to protein (1 sequence/position per motif)
    def __init__(self, locus, I, II, III, IV, V, VI, VII):
        self.locus = locus #Locus tag
        self.I = I #Domain 1
        self.II = II #Domain 2
        self.III = III
        self.IV = IV
        self.V = V
        self.VI = VI
        self.VII = VII #Domain 7
        self.n_motifs = 0 #Number of motifs
        self.bool_list = [0, 0, 0, 0, 0, 0, 0, 0]
    def calculate_n_motifs(self): #Function to calculate number of motifs in the protein
        motif_list = [self.I, self.II, self.III, self.IV, self.V, self.VI, self.VII]
        bool_list = [0 + len(motif) for motif in motif_list]
        self.bool_list = [0] + bool_list
        self.n_motifs = sum(bool_list)
    def __str__(self): #What print(class) returns
        return f'<{self.locus} with {self.n_motifs} motifs>'
    
# =============================================================================
# 0. Function definition    
# =============================================================================
# Adjust this function to make sure that all the motifs are identified
def search_motif(seq_str, num):
    '''Function to search for any motifs by using string patterns
    Input: Amino acid sequence
    Output: Motifs'''
    motif_list = ['', r'[MAVT]D(.{1,1})V(.{1,1})[ND]Q', r'(.{1,1})R(.{1,1})DA(.{2,2})[NFS][MVI][DHNS]',
                  r'H[LI][SHNA][LI][LV]E(.*?)(([WGP](.{3,3})[DGTS])|(SLEAAT))',
                  r'[FIY][VIL][RHN][AV]HD(.{1,2})[AVIS]Q(.{2,2})[IV]',
                  r'[EDA][FL]LL[AG][NSD][DQ][IVE]DNSN[PV]', 
                  r'[LW]G[IVF](.{3,3})[EQ][FLM][AP][PG][HQ]Y',
                  r'[VTS][VTI]PR[VMI][YF]YGD']
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

# =============================================================================
# 0. Input definitions
# =============================================================================
# infolder = '../data/codons'
infile = os.path.expanduser('~') + '/GH_project/add_species/results/domains/ingroup_domains.fasta'
outdir = os.path.expanduser('~') + '/GH_project/add_species/results/motifs'
output_file = 'other_species_GH70_motif_presence.tsv' #'GH70_motif_presence.tsv' #'motif_presence.tsv'
output_file2 = 'other_species_GH70_motifs.tsv' #'GH70_motifs.tsv' #'gtf_motifs.tsv'
prot_list = []

# =============================================================================
# 1. Create output files and directories
# =============================================================================
if not os.path.exists(outdir):
    os.makedirs(outdir)

with open(f'{outdir}/{output_file}', 'w') as out_file, open(f'{outdir}/{output_file2}', 'w') as out_file2:
    out_file.write('locus\ttype\tmotif_I\tmotif_II\tmotif_III\tmotif_IV\tmotif_V\tmotif_VI\tmotif_VII\n')
    out_file2.write('locus\ttype\tmotif_type\tsequence\tstart\tend\ttotal_length\n')

# =============================================================================
# 2. Assign a type to the sequence if possible
# =============================================================================
# OBS! Modify tag assignment for a single input file
tag_type = {}
with open(infile) as handle:
    aa_seqs = SeqIO.parse(handle, 'fasta')
    for fasta in aa_seqs:
        if 'Dsr' in fasta.id and not fasta.id.endswith('_2'):
            tag_type[fasta.id] = 'GS2'
        elif fasta.id.endswith('_2') or 'Brs' in fasta.id:
            tag_type[fasta.id] = 'BRS'
            to_add = re.search(fasta.id, r'a1,[0-9]')
            if to_add != None:
                tag_type[fasta.id] += to_add
        elif 'Asr' in fasta.id:
            tag_type[fasta.id] = 'AS'
        elif 'Rsr' in fasta.id:
            tag_type[fasta.id] = 'RS'
        elif 'Msr' in fasta.id:
            tag_type[fasta.id] = 'MS'
        else: tag_type[fasta.id] = '???'

# for GH_type in GH_types:
#     with open(f'{infolder}/{GH_type}_codon.fna') as handle:
#         aa_seqs = SeqIO.parse(handle, 'fasta')
#         for fasta in aa_seqs:
#             tag_type[fasta.id] = GH_type
#         tag_type['A1003_12540'] = 'GS'

# Loop through the sequences and retrieve the motifs
with open(infile) as handle:
    aa_seqs = SeqIO.parse(handle, 'fasta')
    for fasta in aa_seqs:
        name, sequence = fasta.id, str(fasta.seq) #str(fasta.seq.translate())
        category = tag_type[fasta.id]
        motif_list = ['']
        for i in range(1, 8):
            motif_matches = search_motif(sequence, i)
            motif_list.append(motif_matches)
            for mot in motif_matches:
                mot = motif(*mot)
                if len(mot.seq) > 0:
                    with open(f'{outdir}/{output_file2}', 'a') as out_file2:
                        out_file2.write(f'{name}\t{category}\t{mot.num}\t{mot.seq}\t{mot.start}\t{mot.end}\t{len(sequence)}\n')
        prot = prot_sequence(name, motif_list[1], motif_list[2], motif_list[3], motif_list[4], motif_list[5], motif_list[6], motif_list[7])
        prot.calculate_n_motifs()
        prot_list.append(prot)
        bool_list = list(map(int, prot.bool_list))
        with open(f'{outdir}/{output_file}', 'a') as out_file:
            out_file.write(f'{prot.locus}\t{category}\t{bool_list[1]}\t{bool_list[2]}\t'+
                           f'{bool_list[3]}\t{bool_list[4]}\t{bool_list[5]}\t'+
                           f'{bool_list[6]}\t{bool_list[7]}\n')
        
