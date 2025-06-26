#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:18:40 2021

From the CodeML output, this script creates a tab file
with the dN and dS values for all the comparisons of
core genes, and it calculates statistics for each strain
pair (not general statistics for all core genes).

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
from numpy import mean, median
import re

# =============================================================================
# 1. Set path to inputs and outputs
# =============================================================================

inpath = os.path.expanduser('~') + '/GH_project/all_core/results'
outpath = f'{inpath}/global'
finalfile = f'{outpath}/core_pairwise_metrics.tsv'

if not os.path.exists(outpath): #Create output directory if it doesn't exist
    os.makedirs(outpath)

# =============================================================================
# 2. Define function to parse the CodeML output (.txt)
# =============================================================================

def parse_codeml_output(infile, pairwise_dict):
    '''

    Parameters
    ----------
    infile : str
        Path to the input file.
    pairwise_dict : dict
        Dictionary to be filled with pairwise evolutionary metrics.

    Returns
    -------
    new_pairwise_dict : dict
        Dictionary to store the pairwise evolutionary metrics between strains.

    '''

    status = 0 #Integer used as boolean
    
    new_pairwise_dict = pairwise_dict.copy() #Copy the input dictionary
    
    results = open(infile, 'r') #Open the input file
    print('results: ', results) #Print file information
    
    for i in results: #Loop through lines in the results
        if i.startswith('pairwise comparison, codon frequencies'): #If it is the first line
            status = 1 #Change the status check

        if status == 1: #If status is equal to 1
            if i[0].isdigit(): #If the first value in the line is 0
                line = i.rstrip() #Remove trailing whitespaces
                line2 = re.sub('\(', '', line) #Remove opening parentheses
                line3 = re.sub('\)', '', line2) #Remove closing parentheses
                spaces = line3.split(' ') #Divide the line by spaces
                first = spaces[1] #Get the second substring (locus tag 1)
                second = spaces[4] #Get the fifth substring (locus tag 2)
                first_split = first.split('_') #Split the locus tag 1
                second_split = second.split('_') #Split the locus tag 2

                g1 = first_split[0] #Get strain 1
                g2 = second_split[0] #Get strain 2
                g1 = g1.replace('K2W83', 'DSMZ').replace('LDX55', 'IBH001') #Fix strain name 1
                g2 = g2.replace('K2W83', 'DSMZ').replace('LDX55', 'IBH001') #Fix strain name 2
                
                if g1.startswith('H'): #If the strain 1 name starts with H
                    g1 = g1[:4] + '-' + g1[4:] #Add the -
                if g2.startswith('H'): #Same for strain 2
                    g2 = g2[:4] + '-' + g2[4:]

            if i.startswith('t='): #If the line starts with t=
                line = i.rstrip() #Remove trailing whitespaces
                line1 = re.sub('=', '', line) #Remove the equal symbol
                line2 = re.sub('\s+', '\t', line1) #Replace any number of whitespaces with a tabulation
                tabs = line2.split('\t') #Split the line by tabulations

                dnds = float(tabs[7]) #Retrieve the dN/dS value
                dn = float(tabs[9]) #Retrieve the dN value
                ds = float(tabs[11]) #Retrieve the dS value
                
                pair = tuple(sorted((g1, g2))) #Create a tuple with the strain names in alphabetical order
                print(pair) #Print the tuple
                
                ##Add dN and dS filters here if wanted
                
                if pair not in new_pairwise_dict.keys(): #If the pair is not in the output dictionary
                    new_pairwise_dict[pair] = [[dn], [ds], [dnds]] #Add the pair as a key and the dN, dS and dN/dS as lists
                else: #If the pair is already in the dictionary
                    new_pairwise_dict[pair][0].append(dn) #Append the values to the preexisting lists
                    new_pairwise_dict[pair][1].append(ds)
                    new_pairwise_dict[pair][2].append(dnds)

    return new_pairwise_dict #Return the dictionary with pairwise values
                
# =============================================================================
# 3. Apply the function to input files
# =============================================================================

locus_pairs = {} #Create an empty dictionary
folders = [folder for folder in os.listdir(inpath) if os.path.isdir(f'{inpath}/{folder}') and folder.startswith('A1401')] #Retrieve the folders with core gene alignments
for folder in folders: #Loop through the folders
    print(f'{inpath}/{folder}/{folder}.txt') #Print the path to the CodeML file
    locus_pairs = parse_codeml_output(f'{inpath}/{folder}/{folder}.txt', locus_pairs) #Apply the function to update the input dictionary

        
# =============================================================================
# 4. Save statistics to files
# =============================================================================

with open(finalfile, 'w') as stats_file: #Open the output statistics file
    stats_file.write('strain1\tstrain2\tmean_dN\tmedian_dN\tmean_dS\tmedian_dS\tmean_w\tmedian_w\n') #Write file headers
    for fp in sorted(locus_pairs.keys()): #Loop through the strain pairs
        mean_dN = mean(locus_pairs[fp][0]) #Calculate the mean dN
        median_dN = median(locus_pairs[fp][0]) #Calculate the median dN
        mean_dS = mean(locus_pairs[fp][1]) #Calculate the mean dS
        median_dS = median(locus_pairs[fp][1]) #Calculate the median dS
        mean_w = mean(locus_pairs[fp][2]) #Calculate the mean dN/dS
        median_w = median(locus_pairs[fp][2]) #Calculate the median dN/dS
        stats_file.write(f'{fp[0]}\t{fp[1]}\t{mean_dN:.4f}\t{median_dN:.4f}\t{mean_dS:.4f}\t') #Write the strain names and statistics to file
        stats_file.write(f'{median_dS:.4f}\t{mean_w:.4f}\t{median_w:.4f}\n') #Write the remaining statistics and add a linebreak at the end
