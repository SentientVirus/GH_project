# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:39:34 2021

Script to retrieve the sequence motifs in GH70 genes and the positions of the
motifs in the complete genes. The code generates two output files, one with
motif sequence and position and another one with motif presence/absence.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required Python modules
# =============================================================================

import re, os
from Bio import SeqIO

# =============================================================================
# 1. Define classes to store motif information
# =============================================================================

class motif:
    '''
    Class to store the motif number, sequence and position.
    '''
    
    def __init__(self, num, seq, start, end): #Define motif class
        '''Requires:
            - num: motif number (int, 1-7)
            - seq: motif sequence (str, amino acids)
            - start: start position of the motif in the gene (int, in amino acids)
            - end: end position of the motif in the complete gene (int, in amino acids)'''
        self.num = num #Motif number (1-7)
        self.seq = seq #Motif sequence
        self.start = start #Start position of the motif in a protein sequence
        self.end = end #End position of the motif in the sequence
        
    def __str__(self): #What print(class) returns
        return f'<Motif {self.num} from position {self.start} to {self.end}>'

class prot_sequence:
    ''''
    Class to store the information of all the motifs together associated to the
    locus tag of the protein.
    '''
    def __init__(self, locus, I, II, III, IV, V, VI, VII):
        '''
        Requires:
            - locus: locus tag (str)
            - I-VII: motif information (motif class)
        '''
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
        
    def calculate_n_motifs(self):
        '''
        Function to calculate the number of motifs in the protein
        '''
        motif_list = [self.I, self.II, self.III, self.IV, self.V, self.VI, self.VII] #Save the motifs to a list
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
    motif_list = ['', r'[MAVT]D(.{1,1})V(.{1,1})[ND]Q', r'(.{1,1})R(.{1,1})DA(.{2,2})[NFS][MVI][DHNS]',
                  r'H[LI][SHNA][LI][LV]E(.*?)(([WGP](.{3,3})[DGTS])|(SLEAAT))',
                  r'[FIY][VIL][RHN][AV]HD(.{1,2})[AVIS]Q(.{2,2})[IV]',
                  r'[EDA][FL]LL[AG][NSD][DQ][IVE]DNSN[PV]', 
                  r'[LW]G[IVF](.{3,3})[EQ][FLM][AP][PG][HQ]Y',
                  r'[VTS][VTI]PR[VMI][YF]YGD']
    
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

family = 'GH70'
infolder = os.path.expanduser('~') + f'/GH_project/data/fasta/{family}' #Path to inputs
workdir = os.path.expanduser('~') + '/GH_project/motifs' #Working directory
output_file = f'{workdir}/{family}_motif_presence_pos.tsv' #Path to output file with motif presence/absence information
output_file2 = f'{workdir}/{family}_motifs_pos.tsv' #Path to output file with motif sequence and position
infile = f'{infolder}/complete_{family}_repset.fna' #Input file to search for motifs
GH_types = ['GS1', 'GS2', 'GS3', 'GS4', 'BRS', 'NGB'] #GH70 types to be assigned to the genes
tag_type = {} #Initialize dictionary to assign GH type to each locus tag
prot_list = [] #Initialize list to store protein information

# =============================================================================
# 4. Create/empty output files and write headers
# =============================================================================

with open(output_file, 'w') as out_file, open(output_file2, 'w') as out_file2: #Open output files in write mode (overwrites files)
    out_file.write('locus\ttype\tmotif_I\tmotif_II\tmotif_III\tmotif_IV\tmotif_V\tmotif_VI\tmotif_VII\n') #Write headers
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
    codon_seqs = SeqIO.parse(handle, 'fasta') #Parse file contents
    for fasta in codon_seqs: #Loop through fasta records
        if fasta.id in tag_type.keys() or f'{fasta.id}_1' in tag_type.keys(): #If the locus tag is in the GH type dictionary (or the individual domains are)
            name, sequence = fasta.id, str(fasta.seq.translate()) #Store the locus tag and sequence of the gene
            name = name.replace('K2W83_RS', 'DSMZ_').replace('LDX55', 'IBH001').replace('APS55_RS', 'MP2_') #Make locus tags show the strain name
            if f'{fasta.id}_1' in tag_type.keys(): #If the gene type dictionary contains locus tags ending in _1 for the gene
                category = 'GS2_BRS' #Assign it to the GS2_BRS category (double-domain)
            else: category = tag_type[fasta.id] #Otherwise, simply assign it to a category based on the locus tag
            motif_list = [''] #Create a list with an empty string at index 0
            for i in range(1, 8): #Loop from 1 to 7
                motif_matches = search_motif(sequence, i) #Retrieve each of the motifs
                motif_list.append(motif_matches) #Append the motif to a list
                for mot in motif_matches: #Loop through motif matches
                    mot = motif(*mot) #Convert results to a motif object
                    if len(mot.seq) > 0: #If the length of the matching sequence is > 0
                        with open(output_file2, 'a') as out_file2: #Open the output file in append mode
                            out_file2.write(f'{name}\t{category}\t{mot.num}\t{mot.seq}\t{mot.start}\t{mot.end}\t{len(sequence)}\n') #Write motif sequence and position
            prot = prot_sequence(name, motif_list[1], motif_list[2], motif_list[3], motif_list[4], motif_list[5], motif_list[6], motif_list[7]) #Create a protein object with all the motifs
            prot.calculate_n_motifs() #Calculate the number of motifs in the protein
            prot_list.append(prot) #Add the protein object to a list
            bool_list = list(map(int, prot.bool_list)) #Convert the motif presence/absence list from boolean to integers
            with open(output_file, 'a') as out_file: #Open the output file in append mode
                #Write presence/absence information
                out_file.write(f'{prot.locus}\t{category}\t{bool_list[1]}\t{bool_list[2]}\t'+
                               f'{bool_list[3]}\t{bool_list[4]}\t{bool_list[5]}\t'+
                               f'{bool_list[6]}\t{bool_list[7]}\n')
