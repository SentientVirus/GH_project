# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:39:34 2021
Script to tetrieve the motifs in sequences from other species.
Modified on the 12th of July 2024 from the code to retrieve the motifs from 
GH70 genes. In the original script, an environment with Mafft and Pal2nal was 
used first to get the codon alignments of the sequences. In this case, such an
alignment is not needed as input.
Environment: alignment_tree.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================

import re
import os
from Bio import SeqIO

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
        self.I = I #Domain 1 [[1, seq, start, end]]
        self.II = II #Domain 2 [[2, seq, start, end]]
        self.III = III
        self.IV = IV
        self.V = V
        self.VI = VI
        self.VII = VII #Domain 7 [[7, seq, start, end]]
        self.n_motifs = 0 #Number of motifs
        self.bool_list = [0, 0, 0, 0, 0, 0, 0, 0] #List to store motif presence/absence
    def calculate_n_motifs(self): #Function to calculate number of motifs in the protein
        motif_list = [self.I, self.II, self.III, self.IV, self.V, self.VI, self.VII] #List of lists of lists with motif information
        bool_list = [0 + len(motif) for motif in motif_list] #Add the length (always 1) of each list with motif information to the bool_list
        self.bool_list = [0] + bool_list #Add a 0 in the beginning so that indexes correspond to motif numbers
        self.n_motifs = sum(bool_list) #Calculate the overall number of motifs
    def __str__(self): #What print(class) returns
        return f'<{self.locus} with {self.n_motifs} motifs>'
    
# =============================================================================
# 0. Function definition    
# =============================================================================

def search_motif(seq_str, num):
    '''Function to search for any motifs by using string patterns
    Input: Amino acid sequence
    Output: Motifs'''
    
    #Regular expressions for each motif
    motif_list = ['', r'[MNAVTED]D(.{1,1})V(.{1,1})[ND][QL]', 
                  r'(.{1,1})R(.{1,1})DA(.{2,2})[DNFSHY][MVIYL][DHNSK]',
                  r'[HY][VLI][QSHNATV][LVIY][LVNI]E(.*?)(([WGPS](.{3,3})[MDGTSEV])|(SLEAAT))',
                  r'[FIYM][AMTVSIL][TRHN][SNAV]HD(.{1,2})[RAVISEPT][NKQ](.{2,2})[IVL]',
                  r'[EDA][FLYM][LV][LVI][AG][VNMSD][DQ][LIVE][AD][NL][SQ]N[PVT]', 
                  r'[LW]G[IVF](.{3,3})[GEQW][FLM][AP][PG][HQAS][YF]',
                  r'[MVTS][VTI][PT][TRQ][VMITL][YF]Y[GS]D'] 
    
    motifs = [] #Create empty list to store the motifs
    motif_info = re.finditer(motif_list[num], seq_str) #Find all the possible matches to the motif in the protein
    
    for motif_match in motif_info: #Loop through the matches
        if motif_match is not None: #If there is a match, store new values
            motif_seq = motif_match.group(0) #Convert the match object to a string with the sequence
            start = motif_match.start() + 1 #Set the motif start (index starts at 1)
            end = motif_match.end() + 1 #Set the motif end
        else: #Otherwise use generic values
            motif_seq = '' #Empty sequence
            start = 0 
            end = 0
        motifs.append([num, motif_seq, start, end]) #Add list with motif information to the motifs list
    return motifs

# =============================================================================
# 0. Input definitions
# =============================================================================

infile1 = os.path.expanduser('~') + '/GH_project/add_species/files/blastp_formatted.faa' #Fasta file with full-protein sequences
infile2 = os.path.expanduser('~') + '/GH_project/add_species/files/gtfB-like_outgroups.faa' #Second file with extra sequences
outdir = os.path.expanduser('~') + '/GH_project/add_species/results/motifs' #Output directory
output_file = 'other_species_GH70_motif_presence.tsv' #Output file with motif presence/absence
output_file2 = 'other_species_GH70_motifs.tsv' #Output file with motif sequences and positions
GH_types = ['GS1', 'GS2', 'BRS', 'NGB'] #Types of genes to be considered

# =============================================================================
# 1. Create output files and directories
# =============================================================================
if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create it

with open(f'{outdir}/{output_file}', 'w') as out_file, open(f'{outdir}/{output_file2}', 'w') as out_file2: #Create or overwrite output files
    out_file.write('locus\ttype\tmotif_I\tmotif_II\tmotif_III\tmotif_IV\tmotif_V\tmotif_VI\tmotif_VII\n') #Add headers to the presence/absence file
    out_file2.write('locus\ttype\tmotif_type\tsequence\tstart\tend\ttotal_length\n') #Add header to the file with motif sequences and positions

# =============================================================================
# 2. Assign a type to the sequence if possible
# =============================================================================

tag_type = {} #Dictionary to store gene types
prot_list = [] #List to store protein information

for infile in [infile1, infile2]: #Loop through input files
    with open(infile) as handle: #Read input files
        aa_seqs = SeqIO.parse(handle, 'fasta') #Parse file contents
        for fasta in aa_seqs: #Loop through fasta records in the file
            print(fasta.id) #Print the id of the fasta record
            if 'DSR' in fasta.id: #If the id includes certain annotations
                tag_type[fasta.id] = 'DSR' #Assign the DRS (dextransucrase) gene type
            elif 'BRS' in fasta.id: #Do the same for other gene types
                tag_type[fasta.id] = 'BRS' #Branching sucrase
                to_add = re.search(fasta.id, r'a1,[0-9]') #Search for the subtype of BRS
                if to_add != None:  #If the subtype is indicated
                    tag_type[fasta.id] += to_add #Add it to the BRS tag
            elif 'ASR' in fasta.id: 
                tag_type[fasta.id] = 'ASR' #Alternansucrase
            elif 'RSR' in fasta.id:
                tag_type[fasta.id] = 'RSR' #Reuteransucrase
            elif 'MSR' in fasta.id: 
                tag_type[fasta.id] = 'MSR' #Mutansucrase
            elif '4,6_gtf' in fasta.id:
                tag_type[fasta.id] = 'GTFB_4,6' #4,6-glucanotransferase, GtfB-like
            elif '4,3_gtf' in fasta.id:
                tag_type[fasta.id] = 'GTFB_4,3' #4,3-glucanotransferase, GtfB-like
            else: tag_type[fasta.id] = 'ND' #Not determined (unknown)

            name, sequence = fasta.id, str(fasta.seq) #Get the name and sequence of the gene
            category = tag_type[fasta.id] #Retrieve the gene category
            motif_list = [''] #Create an empty motif list
            for i in range(1, 8): #Loop through the motif numbers
                motif_matches = search_motif(sequence, i) #Search for a motif in the protein
                motif_list.append(motif_matches) #Add the match to the motif list
                for mot in motif_matches: #Loop throuh motif information
                    mot = motif(*mot) #Convert the motif into an object of class motif
                    if len(mot.seq) > 0: #If the length of the motif is > 0
                        if len(motif_list[1]) > 1: #If the length of the list with motif I is > 1
                            category = 'DSR_BRS' #Change the category to DSR_BRS
                        with open(f'{outdir}/{output_file2}', 'a') as out_file2: #Open the output file with motif sequence and position information
                            out_file2.write(f'{name}\t{category}\t{mot.num}\t{mot.seq}\t{mot.start}\t{mot.end}\t{len(sequence)}\n') #Write the motif information into the file
            
            #Create an object of class protein_sequence with all the motifs and the ID of the fasta record as name
            prot = prot_sequence(name, motif_list[1], motif_list[2], 
                                 motif_list[3], motif_list[4], motif_list[5], 
                                 motif_list[6], motif_list[7])
            prot.calculate_n_motifs() #Get the motif presence/absence information
            prot_list.append(prot) #Append the protein sequence to the list (only for debugging purposes)
            
            bool_list = list(map(int, prot.bool_list)) #Get the presence/absence list
            with open(f'{outdir}/{output_file}', 'a') as out_file: #Open the output file with presence/absence information in append-only mode
                out_file.write(f'{prot.locus}\t{category}\t{bool_list[1]}\t{bool_list[2]}\t'+
                               f'{bool_list[3]}\t{bool_list[4]}\t{bool_list[5]}\t'+
                               f'{bool_list[6]}\t{bool_list[7]}\n') #Write bool list as a row in the file
            
