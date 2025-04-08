#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 19:34:41 2024

Script to retrieve domain amino acid sequences of GH70 proteins, including 
proteins from other species that are related to A. kunkeei.
Environment: alignment_tree.yml

@author: Marina Mota-Merlo
"""
# =============================================================================
# 0. Load required libraries
# =============================================================================

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord as seqr
import os

# =============================================================================
# 1. Define paths to input files
# =============================================================================

in_file = os.path.expanduser('~') + '/GH_project/add_species/results/interproscan/blastp_species.tsv' #File with domain annotations
fasta_file = os.path.expanduser('~') + '/GH_project/add_species/files/blastp_formatted.faa' #File with full-protein sequences
outfile = os.path.expanduser('~') + '/GH_project/add_species/results/domains/blastp_domains.faa' #Output file with domain sequences
outdir = os.path.dirname(outfile) #Output directory
GH70_text = 'Glycosyl hydrolase family 70' #String that indicates the presence of a GH70 domain

if not os.path.exists(outdir): #If the output directory does not exist
    os.makedirs(outdir) #Create the output directory
    
# =============================================================================
# 2. Create dictionary with the locus tags and positions of domains
# =============================================================================

domain_dict = {} #Create an empty dictionary to store domain information
with open(in_file) as interpro: #Open input file with domain annotations
    domains = pd.read_csv(interpro, sep = '\t', header = None) #Read the input file as a dataframe
    for index, line in domains.iterrows(): #Loop through the dataframe
        if line[5] == GH70_text: #If the domain is a GH70
            if line[0] not in domain_dict.keys(): #If the locus tag is not in the dictionary keys
                domain_dict[line[0]] = (line[6], line[7]) #Add the domain position to the dictionary
            else: #If the locus tag is already in the dictionary (GS2_BRS)
                domain_dict[line[0]+'_1'] = (line[6], line[7]) #Store the current values in the BRS domain
                domain_dict[line[0]+'_2'] =  domain_dict[line[0]] #Store the past values in the GS2 domain
                del domain_dict[line[0]] #Remove the full gene from dictionary keys
                
# =============================================================================
# 3. Retrieve domain sequences               
# =============================================================================

my_prots = [] #Create an empty list to store the domain sequences    
with open(fasta_file) as prots: #Open the fasta file with input full-protein sequences
    for prot in SeqIO.parse(prots, 'fasta'): #Loop through records in the file
        if f'{prot.id}_1' in domain_dict.keys(): #If the locus tag_1 in the keys of the domain dictionary
            prot.id = prot.id + '_1' #Modify the current protein id to match it
            prot.description = prot.description + '_1' #Do the same to the description
            
        new_prot = seqr(prot.seq, id = prot.id, name = prot.name, description = prot.description) #Create a new protein object
        pos = domain_dict[prot.id] #The position is the value associated to the protein id key
        new_prot.seq = prot.seq[pos[0]:pos[1]+1] #The sequence of the new protein can be retrieved from the original protein using pos
        my_prots.append(new_prot) #Add the new protein to the list
        
        print(f'{new_prot.id} will be added to {outfile}') #Print progress
        if prot.id.endswith('_1'): #If the protein id ends with _1
            new_prot2 = seqr(prot.seq, id = prot.id[:-2] + '_2', name = prot.name, description = prot.description[:-2] + '_2') #Create a protein sequence for the other domain
            pos2 = domain_dict[new_prot2.id] #Retrieve the position of the second domain
            new_prot2.seq = prot.seq[pos2[0]:pos2[1]+1] #Change the protein record sequence into the domain sequence
            my_prots.append(new_prot2) #Add domain to the list
            print(f'{new_prot2.id} will be added to {outfile}') #Print progress
            
# =============================================================================
# 4. Save information to output file in fasta format
# =============================================================================

SeqIO.write(my_prots, outfile, 'fasta')
