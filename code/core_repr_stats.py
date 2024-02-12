#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:18:40 2021

@author: marina
"""
import os
from numpy import mean, median
import pandas as pd


genes = ['GS1', 'GS2', 'BRS', 'NGB']

        
codeml_dir = 'results/core_genes/GH_repr'
#os.chdir(codeml_dir)

for gene in genes:
    dN_list = []
    dS_list = []
    w_list = []
    
    for file in os.listdir(codeml_dir):
        if file.startswith(gene) and file.endswith('dNdS.tsv'):
            dNdS = pd.read_csv(f'{codeml_dir}/{file}', sep = '\t')
            dN_list += list(dNdS['dN'])
            dS_list += list(dNdS['dS'])
            w_list += list(dNdS['w'])
    
    with open(f'{codeml_dir}/global_stats/core_global_stats_{gene}.tsv', 'w') as stats_file:
        full_list = [dN_list, dS_list, w_list]
        rate_name = ['dN', 'dS', 'w']
        for l in range(len(full_list)):
            mean_stat = mean(full_list[l])
            median_stat = median(full_list[l])
            stats_file.write(f'Mean_{rate_name[l]}\t{mean_stat}\n')
            stats_file.write(f'Median_{rate_name[l]}\t{median_stat}\n')

# os.chdir('codeml_out')
# for folder in os.listdir():
#     if os.path.isdir(folder) and ('GS' in folder or 'BRS' in folder or 'short' in folder or 'NCB' in folder):
# #    if os.path.isdir(folder) and 'GH32_' in folder:
#         dNdS = pd.read_csv(f'{folder}/dNdS.tsv', sep = '\t')
#         dN_list += list(dNdS['dN'])
#         dS_list += list(dNdS['dS'])
#         w_list += list(dNdS['w'])
        
        
# with open('global_stats.tsv', 'w') as stats_file:
#     full_list = [dN_list, dS_list, w_list]
#     rate_name = ['dN', 'dS', 'w']
#     for l in range(len(full_list)):
#         mean_stat = mean(full_list[l])
#         median_stat = median(full_list[l])
#         stats_file.write(f'Mean_{rate_name[l]}\t{mean_stat}\n')
#         stats_file.write(f'Median_{rate_name[l]}\t{median_stat}\n')
        
        

