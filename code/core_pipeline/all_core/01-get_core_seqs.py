#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:09:59 2024

Script to retrieve the sequences of core genes (both amino acid and nucleotide)
from a dataset of 106 Apilactobacillus kunkeei strains and save each gene in
a separate file.

@author: Marina Mota-Merlo
"""
# =============================================================================
# 0. Import required modules
# =============================================================================

import os, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# =============================================================================
# 1. Add paths to input files
# =============================================================================

infile = os.path.expanduser('~') + '/GH_project/all_core/SCPOs.txt'
outdir = os.path.expanduser('~') + '/GH_project/all_core/fasta'
inpath = os.path.expanduser('~') + '/Akunkeei_files/cds'
inpath2 = os.path.expanduser('~') + '/Akunkeei_files/faa'
MP2_cds = os.path.expanduser('~') + '/LAB/annotations_220613/MP2.gbk'

# =============================================================================
# 2. Create output directory if it does not exist
# =============================================================================

if not os.path.exists(outdir):
    os.makedirs(outdir)

# =============================================================================
# 3. Retrieve the records from MP2 (reversed compared to the other strains)
# =============================================================================

MP2_records = {} #Create an empty dictionary
check = False #Create a boolean and set it to false
for record in SeqIO.parse(MP2_cds, 'genbank'): #Loop through the records in the MP2 GenBank
    if check == False: #If the boolean is false
        full_seq = record.seq #Retrieve the sequence (chromosome)
        check = True #Set the boolean to true
        for feature in record.features: #Loop through the features in the record
            if 'locus_tag' in feature.qualifiers.keys(): #If the feature includes a locus tag
                loctag = feature.qualifiers['locus_tag'][0] #Retrieve the locus tag
                MP2_records[loctag] = feature.location #Save the position of the record
                print(feature.location.start, feature.location.strand) #Print the start and strand of the record

MP2_records_fna = [] #Create an empty list
for locus in MP2_records.keys(): #Loop through the retrieved locus tags
    position = MP2_records[locus] #Get the position of the feature
    start = position.start #retrieve the start position
    end = position.end #Retrieve the end position
    strand = position.strand #Retrieve the strand
    gene_seq = full_seq[start:end] #Retrieve the feature nucleotide sequence
    if strand == -1: #If the feature is in the reverse strand
        gene_seq = gene_seq.reverse_complement() #Get the reverse complement of the sequence
    record = SeqRecord(gene_seq, id = locus, name = locus, description = locus) #Save the feature information to a SeqRecord object
    MP2_records_fna.append(record) #Add the record to the list

# =============================================================================
# 4. Retrieve the records from all genomes
# =============================================================================

SCPO_dict = {} #Create an empty dictionary
with open(infile) as SCPOs: #Loop through the file with SCPO information
    for SCPO in SCPOs: #Loop through the lines in the file
        scpo_ids = SCPO.split(' ') #Split by spaces
        scpo_list = [scpon.split('|')[1] for scpon in scpo_ids] #Split by | and get the second values
        pattern = '[AG]\d{4}' #Define the strings to retrieve
        scpo_filtered = [scpo for scpo in scpo_list if scpo.startswith('H') or re.search(pattern, scpo) != None or scpo.split('_')[0] in ['LDX55', 'K2W83', 'MP2', 'fhon2']] #Get the SCPOs of the 106 strains only
        scpo_name = scpo_filtered[6] #Get the name of the SCPO (based on strain A1401)
        SCPO_dict[scpo_filtered[0]] = scpo_filtered #Assign all the sequences to the same dictionary key (locus tag in strain A0901)
        outfile = f'{outdir}/{scpo_name}.fna' #Define the name of the output file
        outfile2 = outfile.replace('.fna', '.faa') #Define the name of the output protein sequence file
        with open(outfile, 'w') as out_handle: #Open the output file in write mode
            out_handle.write('') #Create the file or remove previous file content 
            
        record_list = [] #Create an empty list to store the nucleotide sequences
        prot_list = [] #Create an empty list to store the protein sequences
        
        for scpo in scpo_filtered: #Loop through SCPO locus tags
            strain = scpo.split('_')[0] #Get strain name from locus tag
            if strain == 'K2W83': #If the first part of the locus tag is K2W83
                strain = 'DSMZ12361' #Set the strain name to DSMZ12361
            elif strain == 'LDX55': #Same, but for strain IBH001
                strain = 'IBH001'
            elif strain == 'fhon2': #Same, but for Fhon2
                strain = 'Fhon2'
            if strain != 'MP2': #If the strain is not MP2
                infile = f'{inpath}/{strain}_cds_from_genomic.fna' #Read the nucleotide fasta file
                    
                for record in SeqIO.parse(infile, format = 'fasta'): #Loop through nucleotide records
                    locus_tag = record.description.split('locus_tag=')[1].split(']')[0] #Get the locus tag from the record
                    if locus_tag.startswith('AKU'): #If the locus tag starts with AKU
                        locus_tag = locus_tag[3:] #Remove it from the locus tag
                    if locus_tag == scpo.replace('-', '').upper(): #If the locus tag is the same as the SCPO name without - in uppercase
                        print(locus_tag, scpo) #Print the locus tag and SCPO
                        record.id = locus_tag #Set the record id to the locus tag
                        record_list.append(record) #Add the record to a list

            else: #If the strain is MP2
                for record in MP2_records_fna: #Loop through the MP2 records
                    if record.id == scpo: #If the record id is the SCPO name
                        print(record.id, scpo) #Print both
                        record_list.append(record) #Add the record to the list

        SeqIO.write(record_list, outfile, 'fasta') #Write the records to the nucleotide output file
        
        for record in record_list: #Loop through the records
            record.seq = record.seq.translate(stop_symbol = '') #Translate the records
        SeqIO.write(record_list, outfile2, 'fasta') #Write the translated records to a file
