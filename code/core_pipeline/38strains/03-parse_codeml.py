#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Created on Tue June 21 2022

Script from https://github.com/faylward/dnds, slightly modified by Julia 
Pedersen and Marina Mota. It parses the CodeML output and calculates some
statistics for each subset of genes.

"""

# =============================================================================
# 0. Import required packages
# =============================================================================

import re, os
from numpy import mean, median


# =============================================================================
# 1. Load inputs
# =============================================================================

indir = os.path.expanduser('~') + '/GH_project/all_core/38_strains/results'
basal_strains = ['A1001', 'A1404'] #Strains to exclude

# =============================================================================
# 2. Define function to parse CodeML output
# =============================================================================

def parse_codeml_output(infile, outfile1, outfile2):

    dn_list = [] #Empty list to store dN values
    ds_list = [] #Empty list to store dS values
    dnds_list = [] #Empty list to store dN/dS values

    status = False #Boolean to indicate when to start parsing the file

    outfile = open(outfile1, 'w') #Overwrite output file

    outfile.write('locus1\tlocus2\tdN\tdS\tw\n') #Add header to output file

    results = open(infile, 'r') #Open input file in read-only mode
    
    for i in results:
        if i.startswith('pairwise comparison, codon frequencies'):
            status = True #Start parsing the file

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

#                if ds < 50 and ds > 0.01 and dnds < 99:
#                    outfile.write(f'{first}\t{second}\t{dn}\t{ds}\t{dnds}\n')
                    
                if not any(basal in loctag.split('_')[0] for basal in basal_strains for loctag in [first, second]): #Filter out pairs that include basal strains
                    outfile.write(f'{second}\t{first}\t{dn}\t{ds}\t{dnds}\n') #Write pairwise metrics to file
                    print(f'Added dN, dS and w for pair {first} vs {second}.') #Print progress
                
                    ds_list.append(ds) #Basal values are not considered for global metrics
                    dn_list.append(dn)
                    dnds_list.append(dnds)
    
    with open(outfile2, 'w') as outno2: #Open the other output file in write mode (overwrites the file)
        filtered_ds = [ds for ds in ds_list if ds >= 0] #Create a list of dS values
        filtered_dn = [dn_list[i] for i in range(len(ds_list)) if ds_list[i] >= 0] #Create a list of dN values
        filtered_dnds = [dnds_list[j] for j in range(len(ds_list)) if ds_list[j] >= 0] #Create a list of dN/dS values
        print(f'Filtered dS:\n{filtered_ds}')
        outno2.write(f'Mean_dN\t{mean(filtered_dn):.4f}\n') #Write mean dN to file
        outno2.write(f'Median_dN\t{median(filtered_dn):.4f}\n') #Write median dN to file
        outno2.write(f'Mean_dS\t{mean(filtered_ds):.4f}\n') #Write mean dS to file
        outno2.write(f'Median_dS\t{median(filtered_ds):.4f}\n') #Write median dS to file
        outno2.write(f'Mean_w\t{mean(filtered_dnds):.4f}\n') #Write mean dN/dS to file
        outno2.write(f'Median_w\t{median(filtered_dnds):.4f}\n') #Write median dN/dS to file

# =============================================================================
# 3. Run function in list comprehension
# =============================================================================

infiles = [f'{indir}/{core_gene}/{core_gene}.txt' for core_gene in os.listdir(indir) if os.path.isdir(f'{indir}/{core_gene}')]
outfiles1 = [f'{indir}/{core_gene}/dNdS.tsv' for core_gene in os.listdir(indir) if os.path.isdir(f'{indir}/{core_gene}')]
outfiles2 = [f'{indir}/{core_gene}/stats.tsv' for core_gene in os.listdir(indir) if os.path.isdir(f'{indir}/{core_gene}')]

[parse_codeml_output(infiles[i], outfiles1[i], outfiles2[i]) for i in range(len(infiles))]

