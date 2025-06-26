#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:18:40 2021

This script calculates general statistics for all core genes based
on the files where pairwise substitution rates were calculated.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os
from numpy import mean, median
import pandas as pd

# =============================================================================
# 1. Define input variables and paths
# =============================================================================

codeml_dir = os.path.expanduser('~') + '/GH_project/all_core/results' #Input folder
outfile = f'{codeml_dir}/global/core_general_stats.tsv'

dN_list = [] #Create empty list to store dN values
dS_list = [] #Create empty list to store dS values
w_list = [] #Create empty list to store dN/dS values
folders = [folder for folder in os.listdir(codeml_dir) if os.path.isdir(f'{codeml_dir}/{folder}') and folder.startswith('A1401')]
for folder in folders:    
    for file in os.listdir(f'{codeml_dir}/{folder}'): #Loop through files in input folder
        if file.endswith('dNdS.tsv'): #If the files have dN and dS data
            dNdS = pd.read_csv(f'{codeml_dir}/{folder}/{file}', sep = '\t') #Read the files as dataframes
            dN_list += list(dNdS['dN']) #Add the remaining rows to the dN list
            dS_list += list(dNdS['dS']) #Add the remaining rows to the dS list
            w_list += list(dNdS['w']) #Add the remaining rows to the dN/dS list
    
    with open(outfile, 'w') as stats_file: #Open the output file
        full_list = [dN_list, dS_list, w_list] #Create a list of lists with the parameters to write to the file
        rate_name = ['dN', 'dS', 'w'] #Create a list with the names of the rates (column headers)
        for l in range(len(full_list)): #Loop through the length of the lists
            mean_stat = mean(full_list[l]) #Get the mean of each statistic
            median_stat = median(full_list[l]) #Get the median of each statistic
            stats_file.write(f'Mean_{rate_name[l]}\t{mean_stat}\n') #Write the mean of the statistic to the output file
            stats_file.write(f'Median_{rate_name[l]}\t{median_stat}\n') #Write the median of the statistic to the output file

        
        

