#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:38:36 2024

Script to parse CodeML results and retrieve the locus tags of the core genes
that have dS > 10.

@author: Marina Mota-Merlo
"""
# =============================================================================
# 0. Load required libraries
# =============================================================================
import re, os

# =============================================================================
# 1. Function to parse CodeML results
# =============================================================================
def parse_codeml_output(infile, outfile):
    
    #Strains not to be considered
    to_exclude = ['fhon13', 'A1001', 'A1404', 'A2002', 'A2003', 'EFB6', 
                   'FF306', 'LAan', 'LAce', 'LAdo', 'LAfl', 'LAko', 'LAla', 
                   'LAni', 'LAnu', 'LMbe', 'LMbo', 'YH15']
    
    check = False #Boolean to check if the gene should be appended to a list
    
    filename = os.path.basename(infile) #Get base name of the input file
    filename = filename.replace('.txt', '') #Remove file extension (only gene name remains)

    status = False #Boolean to indicate when to start parsing the file

    results = open(infile, 'r') #Open input file in read mode
    
    for i in results: #Loop through the lines in the file
        if i.startswith('pairwise comparison, codon frequencies'): #Until the right line is reached
            status = True #Boolean set to True to start parsing the file

        if status: #If status is set to True
            if i[0].isdigit(): #If the first character in the line is a digit
                line = i.rstrip() #Remove spaces at the end of the line
                line2 = re.sub('\(', '', line) #Remove opening parentheses
                line3 = re.sub('\)', '', line2) #Remove closing parentheses
                spaces = line3.split(' ') #Separate the line string by spaces to create a list
                first = spaces[1] #First locus tag (second element of the list)
                second = spaces[4] #Second locus tag (fifth element of the list)

            if i.startswith('t='): #If the line starts with t=
                line = i.rstrip() #Remove spaces at the end of the line
                line0 = re.sub('=', '', line) #Remove equal signs from the line
                line1 = re.sub('\s+', '\t', line0) #Remove any non-numerican characters
                tabs = line1.split('\t') #Create a list by splitting the line by tabulations
                dnds = float(tabs[7]) #dN/dS is the eigth element of the list
                dn = float(tabs[9]) #dN is the tenth element of the list
                ds = float(tabs[11]) #dS is the twelfth element of the list
                
                val_list = [first.split('_')[0], second.split('_')[0]] #Get the strain names from the locus tags
                
                if not any(excl in val_list for excl in to_exclude) and ds > 10: #Filter out basal and incomplete strains, and get pairs with ds > 10
                    with open(outfile, 'a') as handle: #Open the output file in append mode
                        print(f'Pair of {first}, {second} of gene {filename} has ds {ds}') #Print pair information
                        file_info = f'{first}\t{second}\t{dn}\t{ds}\t{dnds}\n' #Create string to be added to file
                        handle.write(file_info) #Write string to file
                    check = True #Boolean to add gene to a list
    
    if check: #If there is at least one pair with dS > 10
        gene = filename #Save the gene locus tag to a variable
    else: gene = '' #Otherwise, the gene variable is an empty string
                
    return gene #Return gene variable


# =============================================================================
# 2. Define inputs and outputs
# =============================================================================
infolder = os.path.expanduser('~') + '/GH_project/results/core_genes' #Folder with the input files
outfolder = f'{infolder}/high_dS' #Folder to save the output
outfile = f'{outfolder}/pairwise_high_dS.tsv' #Output file

# =============================================================================
# 3. Create outputs
# =============================================================================
if not os.path.exists(outfolder): #If the output folder does not exist
    os.makedirs(outfolder) #Create it
    
with open(outfile, 'w+') as handle: #Open output file
    handle.write('locus1\tlocus2\tdN\tdS\tw\n') #Remove all file contents if there is any, add headers

# =============================================================================
# 4. Implement parsing function
# =============================================================================
genes = [] #List of genes with at least one pair with dS > 10
for file in os.listdir(infolder): #Loop through the contents of the input folder
    if os.path.isdir(f'{infolder}/{file}') and file.startswith('A1401'): #If the element is a subfolder whose name starts with A1401
        infile = f'{infolder}/{file}/{file}.txt' #Get the path to the input file
        if f'{file}.txt' in os.listdir(f'{infolder}/{file}'): #If the input file is in the subfolder
            gene = parse_codeml_output(infile, outfile) #Run the parsing function
            if gene != '': #If the gene variable is not an empty string
                genes.append(gene) #Append the variable to the list
