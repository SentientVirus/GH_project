#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 15:46:33 2025

Script to parse ClonalFrameML results and count the total lengths of the 
genomes in the dataset.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os

# =============================================================================
# 2. Define path to inputs and phylogroup divisions
# =============================================================================

workdir = os.path.expanduser('~') + '/GH_project' #Path to the working directory
phylogroups = ['A', 'B-C'] #Phylogroup designations

# =============================================================================
# 2. Create a function to parse the .xmfa files to calculate the final genome
# lengths used by ClonalFrameML.
# =============================================================================

def get_lengths(infile):
    '''
    Function to calculate the total length of each genome in an alignment,
    excluding gaps.

    Parameters
    ----------
    infile : str
        Path to the input file.

    Returns
    -------
    length_dict : dict
        Dictionary with strain names (str) as keys and genome lengths (int) as
        values.
    '''
    
    length_dict = {} #Create empty dictionary
    with open(infile) as handle: #Open input file
        #Loop through the lines in the file, avoiding the ones that don't contain sequence data
        file_reader = [line for line in handle if not line.startswith('#') and not line.startswith('=')]
        for l in file_reader: #Loop through the retained lines
            if l[0] == '>': #If the line starts with >
                strain = l.split()[-1].replace('.fna', '') #Retrieve the strain name
            if strain not in length_dict: #If the strain name is not in the dictionary
                length_dict[strain] = 0 #Include it with an initial length of 0
            elif l[0] != '>': #If the line doesn't start with >
                length_dict[strain] += len(l.replace('-', '').strip()) #Add the length of the segment to the total for the strain
    return length_dict #Return the dictionary

# =============================================================================
# 3. Loop through the phylogroups and save the lengths to files.
# =============================================================================

for i in phylogroups: #Loop through phylogroups
    infile = f'{workdir}/xmfa_sslcb500_genome_len/phylogroup{i}/AK38gh_phylogr{i}_sslcb500_mod.fna.xmfa' #Path to the input file
    outfile = infile.replace('_sslcb500_mod.fna.xmfa', '_lengths.tab') #Path to the output file

    length_dict = get_lengths(infile) #Apply the function to get the lengths
    with open(outfile, 'w') as handle: #Open the output file
        handle.write('Strain\tGenome_length(bp)\n') #Write headers
        #Write strain names and their corresponding lengths
        [handle.write(f'{strain.replace("DSM", "DSMZ")}\t{length_dict[strain]}\n') for strain in sorted(length_dict.keys())]