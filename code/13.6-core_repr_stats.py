#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:18:40 2021

This script calculares general statistics for all core genes based
on the files where pairwise substitution rates were calculated.

@author: Marina Mota-Merlo
"""
import os
from numpy import mean, median
import pandas as pd

# =============================================================================
# 1. Define input variables and paths
# =============================================================================
genes = ['GS1', 'GS2', 'BRS', 'NGB'] #Genes whose statistics will be summarized
to_exclude = ['A1001', 'A1003', 'A1404'] #Strains to exclude from GS1
        
codeml_dir = os.path.expanduser('~') + '/GH_project/results/core_genes/GH_repr' #Input folder

for gene in genes: #Loop through gene types
    dN_list = [] #Create empty list to store dN values
    dS_list = [] #Create empty list to store dS values
    w_list = [] #Create empty list to store dN/dS values
    
    for file in os.listdir(codeml_dir): #Loop through files in input folder
        if file.startswith(gene) and file.endswith('dNdS.tsv'): #If the files have dN and dS data
            dNdS = pd.read_csv(f'{codeml_dir}/{file}', sep = '\t') #Read the files as dataframes
            if gene == 'GS1': #If the gene is GS1
                for strain in to_exclude: #Loop through strains to exclude
                    index_out = dNdS[dNdS['locus1'].str.contains(strain) | dNdS['locus1'].str.contains(strain)].index #Get the index of rows with strains to exclude
                    dNdS.drop(index_out, inplace = True) #Remove these rows
            dN_list += list(dNdS['dN']) #Add the remaining rows to the dN list
            dS_list += list(dNdS['dS']) #Add the remaining rows to the dS list
            w_list += list(dNdS['w']) #Add the remaining rows to the dN/dS list
    
    with open(f'{codeml_dir}/global_stats/core_global_stats_{gene}.tsv', 'w') as stats_file: #Open the output file
        full_list = [dN_list, dS_list, w_list] #Create a list of lists with the parameters to write to the file
        rate_name = ['dN', 'dS', 'w'] #Create a list with the names of the rates (column headers)
        for l in range(len(full_list)): #Loop through the length of the lists
            mean_stat = mean(full_list[l]) #Get the mean of each statistic
            median_stat = median(full_list[l]) #Get the median of each statistic
            stats_file.write(f'Mean_{rate_name[l]}\t{mean_stat}\n') #Write the mean of the statistic to the output file
            stats_file.write(f'Median_{rate_name[l]}\t{median_stat}\n') #Write the median of the statistic to the output file

        
        

