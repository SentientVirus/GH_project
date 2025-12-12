#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:20:09 2022

This script reads a fasta file with nucleotide information and calculates the
overall percentage of identity between every pair of sequences in the file.

@author: Marina Mota-Merlo

"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os, re
from Bio import SeqIO
import multiprocessing 
import time 
from functools import partial
from Bio.Emboss.Applications import NeedleCommandline

# =============================================================================
# 1. Define input files and variables
# =============================================================================

np = 36
start_time = time.time()  
home = os.path.expanduser('~')
pathname = f'{home}/GH_project/data/fasta'
GH_type = ['GH70', 'GH32']
out_path = f'{home}/GH_project/tables'
suffixes = ['repset', 'all']

#Create output directory if it doesn't exist
if not os.path.exists(out_path):
    os.makedirs(out_path)
    
# =============================================================================
# 2. Create the functions to generate outputs
# =============================================================================

def needle_align_code(query_seq, target_seq):
    '''Function to perform a Needle alignment between a pair of sequences
    and retrieve the percentage of identity
    
    Inputs:
        - query_seq (str): One of the two sequences to compare
        - target_seq (str): The second sequence'''
        
    #Define parameters to run Needle
    needle_cline = NeedleCommandline(asequence="asis:" + query_seq,
                                      bsequence="asis:" + target_seq,
                                      sprotein=True,
                                      aformat="simple",
                                      gapopen=10,
                                      gapextend=0.5,
                                      outfile='stdout'
                                      )
    out_data, err = needle_cline() #Run Needle
    out_split = out_data.split("\n") #Dive the output into lines
    p = re.compile("\((.*)\)") #Define a regular expression to parse the information in each line
    return p.search(out_split[25]).group(1).replace("%", "") #Retrieve the percentage of identity

def GH_type_dict(GH_file: str, GH_dict: dict):
    '''Function to parse a fasta file, retrieve the headers (locus tags) and
    assign them a gene type.
    
    Inputs:
        - GH_file (str): Fasta file with biological sequences
        - GH_dict (dict): Dictionary to which the sequences should be added'''
        
    GH_type = os.path.basename(GH_file).split('_')[0] #Retrieve gene type from filename
    with open(GH_file) as GH_list: #Open the fasta file
        for record in SeqIO.parse(GH_list, 'fasta'): #Loop through the record in the file
            GH_dict[record.id] = GH_type #Assign the gene type to each record ID (locus tag)
            
def GH_assign_subtypes(pathname: str, GH_fam: str, GH_dict: dict):
    '''Function to create a dictionary where the keys and locus tags and the
    values are gene subtypes.
    
    Inputs:
        - pathname (str): Path to files.
        - GH_fam (str): Gene family (GH70 or GH32)
        - GH_dict (dict): Dictionary to be filled'''
        
    if GH_fam == 'GH70': #If the gene family is GH70
        GH_subtypes = ['GS1', 'GS2', 'GS3', 'GS4', 'BRS', 'NGB', 'short'] #Define subtypes
    elif GH_fam == 'GH32': #If the gene family is GH32
        GH_subtypes = ['S1', 'S2a', 'S2b', 'S3'] #Define subtypes
    else: raise NameError() #Else, raise an error
    for subtype in GH_subtypes: #Loop through the subtypes
        file = f'{pathname}/{subtype}_all.faa' #Open the file with the subtype
        GH_type_dict(file, GH_dict) #Parse it to assign the subtype to locus tags
    return GH_dict
    
def create_out_tab(path: str, suffix: str):
    '''Function to create a tab file to store the percentages of identity and
    write the header of the file.
    
    Inputs:
        - path (str): Path to the output file
        - suffix (str): Suffix of the file (complete dataset vs representative strains)'''
        
    out_file = f'percentage_identity_{suffix}.tab' #Define output file name
    with open(f'{path}/{out_file}', 'w') as out: #Open output file
        out.write('locus1\tlocus2\tGH_type\t%identity\n') #Write header
                        
def process_input(infile: str, path: str, suffix: str, GH_dict: dict, pos_dict: dict):
    '''Function to calculate the percentages of identity and write them to the
    output file.
    
    Inputs:
        - in_list (list): List of input files
        - path (str): Path to output file
        - suffix (str): Suffix of output file
        - GH_dict (dict): GH type dictionary
        - pos_dict (dict): Dictionary to store the percentages of identity'''
        
    out_file = f'{path}/percentage_identity_{suffix}.tab' #Path to output file
    with open(infile) as handle: #Open output file
        record_list = [record for record in SeqIO.parse(handle, 'fasta') if record.id in GH_dict.keys()] #Get a list of records
        for i in range(len(record_list)-1): #Loop through the records
            for j in range(i+1, len(record_list)): #Loop again to retrieve all pairs
                GH_name = f'{GH_dict[record_list[i].id]}/{GH_dict[record_list[j].id]}' #Get the types of each gene
                GH_name = GH_name.replace('BRS', 'BrS').replace('short', 'PTL') #Replace certain names
                perc_ident = needle_align_code(record_list[i].seq, record_list[j].seq) #Calculate the percentage of identity
                pos_dict[(record_list[i].id, record_list[j].id)] = perc_ident #Store them in a dictionary
                with open(out_file, 'a') as out: #Open output file
                    out.write(f'{record_list[i].id}\t{record_list[j].id}\t{GH_name}\t{perc_ident}\n') #Write the results to the file
    return out_file
    
# =============================================================================
# 3. Apply the functions and generate outputs
# =============================================================================

[create_out_tab(out_path, suffix) for suffix in suffixes] #Create output files

GH_dict = {} #Create an empty dictionary
for GH in GH_type: #Loop through GH types
    workdir = f'{pathname}/{GH}' #Get the working directory
    GH_dict = GH_assign_subtypes(workdir, GH, GH_dict) #Create the dictionary to assign gene subtypes
    for suffix in suffixes: #Loop through the suffixes
        pos_dict = {} #Create an empty dictionary
        
        #Retrieve all inputs
        inputs = [f'{workdir}/{file}' for file in os.listdir(workdir) if file.endswith(f'{suffix}.fna') and 'complete' not in file and ('GH' in file) and 'GH70_functional' not in file]
        print(inputs)
        
        args = [out_path, suffix, GH_dict, pos_dict] #Create a list with function arguments
        
        if __name__ == '__main__': #Run the function to process inputs with multiprocessing
            pool = multiprocessing.Pool() 
            pool = multiprocessing.Pool(processes=np)
            outputs = pool.map(partial(process_input, path = out_path, suffix = suffix, GH_dict = GH_dict, pos_dict = pos_dict), inputs)
            print("Input: {}".format(inputs))
            print("Output: {}".format(outputs))
    
end_time = time.time() - start_time #Calculate end time
print(f'This script took {end_time:2f} with {np} processes') #Print end time and number of processes
