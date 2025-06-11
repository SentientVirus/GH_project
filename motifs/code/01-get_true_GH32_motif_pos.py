# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:39:34 2021

Script to retrieve the sequence motifs in GH70 genes and the positions of the
motifs in the complete genes. The code generates two output files, one with
motif sequence and position and another one with motif presence/absence.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================

import re, os
from Bio import SeqIO

# =============================================================================
# 1. Create object classes to store motifs and protein information
# =============================================================================

class motif:
    '''
    Class to store the motif number, sequence and position.
    '''
    
    def __init__(self, num, seq, start, end): #Define motif class
        '''
        Requires:
            - num: motif number (int, 1-8)
            - seq: motif sequence (str, amino acids)
            - start: start position of the motif in the gene (int, in amino acids)
            - end: end position of the motif in the complete gene (int, in amino acids)
        '''
        self.num = num #Motif number (1-8)
        self.seq = seq #Motif sequence
        self.start = start #Start position of the motif in a protein sequence
        self.end = end #End position of the motif in the sequence
        
    def __str__(self): #What print(class) returns
        return f'<Motif {self.num} from position {self.start} to {self.end}>'

class prot_sequence: #All motifs associated to protein (1 sequence/position per motif)
    '''
    Class to store the information of all the motifs together associated to the
    locus tag of the protein.
    '''
    
    def __init__(self, locus, I, II, III, IV, V, VI, VII, VIII):
        '''
        Requires:
            - locus: locus tag (str)
            - I-VIII: motif information (motif class)
        '''
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
        self.bool_list = [0, 0, 0, 0, 0, 0, 0, 0, 0] #Motif presence/absence
        
    def calculate_n_motifs(self):
        '''
        Function to calculate the number of motifs in the protein
        '''
        motif_list = [self.I, self.II, self.III, self.IV, self.V, self.VI, self.VII, self.VIII] #Save the motifs to a list
        bool_list = [0 + len(motif) for motif in motif_list] #Update the count of each motif
        self.bool_list = [0] + bool_list #Add 0 at the beginning so that index corresponds with motif number
        self.n_motifs = sum(bool_list) #Get the total number of motifs
        
    def __str__(self): #What print(class) returns
        return f'<{self.locus} with {self.n_motifs} motifs>'
    
# =============================================================================
# 2. Define the function to retrieve the motifs
# =============================================================================

def search_motif(seq_str, num): #Function to search for any motif
    '''
    Function that takes a protein sequence and a motif number and returns
    the motif sequence and position.
    Inputs:
        seq_str: Protein sequence (str) in amino acids
        num: Motif number (int)
    Output:
        motifs: List with motif number (int), sequence (str), start (int) and 
        end (int) positions.
    '''
    #Regular expressions for each motif
    motif_list = ['', r'[FY]SND[IV]QS', r'WY[HG][AVL][TIK]T[KT]DF',
                  r'A[TS]G[SAT][VI][LVIF](.{1,1})N',
                  r'FR(.{0,1})[VIP][VYQ][VIQ](.{2,2})',
                  r'ETP[ILN][VI][KRQ]T', r'D[DAF](.{2,3})D(.{1,1})YA', 
                  r'[NG]W[SL]G[NDS]W[LIV]Y[STW]', r'DN[SAT][TS][LIV]E']
    
    motifs = [] #Initialize list with motif information

    motif_info = re.finditer(motif_list[num], seq_str) #Search for motif sequence
    for motif_match in motif_info: #Loop through the matches
        if motif_match is not None: #If there is a match, store new values
            motif_seq = motif_match.group(0) #Retrieve motif sequence
            start = motif_match.start() + 1 #Add one to the start (index starts at 0, but positions start at 1)
            end = motif_match.end() + 1 #Add one to the end
        else: #Otherwise use generic values
            motif_seq = ''
            start = 0
            end = 0
        motifs.append([num, motif_seq, start, end]) #Save the information to a list
    return motifs #Return list with motif information

# =============================================================================
# 3. Define input variables
# =============================================================================

family = 'GH32' #Gene family
infolder = os.path.expanduser('~') + f'/GH_project/data/fasta/{family}' #Path to inputs
workdir = os.path.expanduser('~') + '/GH_project/motifs' #Working directory
output_file = f'{workdir}/{family}_motif_presence_pos.tsv' #Path to output file with motif presence/absence information
output_file2 = f'{workdir}/{family}_motifs_pos.tsv' #Path to output file with motif sequence and position
infile = f'{infolder}/complete_{family}_repset.fna' #Input file to search for motifs
GH_types = ['S1', 'S2a', 'S2b', 'S3'] #GH32 types to be assigned to the genes
tag_type = {} #Initialize dictionary to assign GH type to each locus tag
prot_list = [] #Initialize list to store protein information

# =============================================================================
# 4. Create/empty output files and write headers
# =============================================================================

with open(output_file, 'w') as out_file, open(output_file2, 'w') as out_file2: #Open output files in write mode (overwrites files)
    out_file.write('locus\ttype\tmotif_I\tmotif_II\tmotif_III\tmotif_IV\tmotif_V\tmotif_VI\tmotif_VII\tmotif_VIII\n') #Add headers
    out_file2.write('locus\ttype\tmotif_type\tsequence\tstart\tend\ttotal_length\n')

# =============================================================================
# 5. Assign a gene type to each locus tag
# =============================================================================

for GH_type in GH_types: #Loop through GH types
    with open(f'{infolder}/{GH_type}_repset.fna') as handle: #Open .fna files with fasta sequences and locus tags
        aa_seqs = SeqIO.parse(handle, 'fasta') #Parse file
        for fasta in aa_seqs: #Loop through fasta records
            tag_type[fasta.id] = GH_type #Assign the right GH type to the locus tag (ID of the record)

# =============================================================================
# 6. Search for motifs in GH32 sequences and save results to the output files
# =============================================================================

with open(infile) as handle: #Open input file with all GH32 nucleotide sequences
    codon_seqs = SeqIO.parse(handle, 'fasta') #Parse the file contents
    for fasta in codon_seqs: #Loop through fasta records in the file
        name, sequence = fasta.id, str(fasta.seq.translate()) #Retrieve the gene ID and sequence
        if 'K2W83' in name: #If the locus tag comes from DSMZ 12361
            name = name.replace('K2W83_RS', 'DSMZ_') #Modify it to show the strain name
        elif 'LDX55' in name: #Same for IBH001
            name = name.replace('LDX55', 'IBH001')
        elif 'APS55' in name: #Same for MP2
            name = name.replace('APS55_RS', 'MP2_')
        category = tag_type[fasta.id] #Retrieve gene type
        motif_list = [''] #Create list with empty string to store motif information
        for i in range(1, 9): #Loop from 1 to 8 to retrieve the 8 motifs in GH32
            motif_matches = search_motif(sequence, i) #Retrieve motif information
            motif_list.append(motif_matches) #Add motif information to another list
            for mot in motif_matches: #Loop through the motif matches that were found
                mot = motif(*mot) #Convert them to motif objects
                if len(mot.seq) > 0: #If the motif length is bigger than 0
                    with open(output_file2, 'a') as out_file2: #Open the output file to store motif sequence and position
                        out_file2.write(f'{name}\t{category}\t{mot.num}\t{mot.seq}\t{mot.start}\t{mot.end}\t{len(sequence)}\n') #Write motif information to file
        prot = prot_sequence(name, motif_list[1], motif_list[2], motif_list[3], motif_list[4], motif_list[5], motif_list[6], motif_list[7], motif_list[8]) #Create a protein sequence object with all the motifs
        prot.calculate_n_motifs() #Get number of motifs
        prot_list.append(prot) #Add protein to a list
        bool_list = list(map(int, prot.bool_list)) #Create a presence/absence list by converting boolean values to 0 or 1 integers
        with open(output_file, 'a') as out_file: #Open output file
            #Write presence/absence information to the file
            out_file.write(f'{prot.locus}\t{category}\t{bool_list[1]}\t{bool_list[2]}\t'+
                           f'{bool_list[3]}\t{bool_list[4]}\t{bool_list[5]}\t'+
                           f'{bool_list[6]}\t{bool_list[7]}\n')
        
