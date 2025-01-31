#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:18:40 2021

From the CodeML output, this script creates a tab file
with the dN and dS values for all the comparisons of
core genes, and it calculates statistics for each strain
pair (not general statistics for all core genes).

@author: Marina Mota Merlo
"""
import os
from numpy import mean, median
import re

outpath = 'results/core_genes'
finalfile = 'core_pairwise_metrics.tsv'

##Parse codeml output (.txt)

def parse_codeml_output(infile, pairwise_dict):

    status = 0

    new_pairwise_dict = pairwise_dict.copy()

    # outfile = open(outfile1, 'w')

    # outfile.write('locus1\tlocus2\tdN\tdS\tw\n')

    results = open(infile, 'r')

    print('results: ', results)
    
    for i in results:

        if i.startswith('pairwise comparison, codon frequencies'):

            status = 1



        if status == 1:

            if i[0].isdigit():

                line = i.rstrip()

                line2 = re.sub('\(', '', line)
                
                line3 = re.sub('\)', '', line2)

                spaces = line3.split(' ')

                first = spaces[1]

                second = spaces[4]
                
                first_split = first.split('_')

                second_split = second.split('_')

                g1 = first_split[0]

                g2 = second_split[0]




            if i.startswith('t='):

                line = i.rstrip()

                line1 = re.sub('=', '', line)

                line2 = re.sub('\s+', '\t', line1)

                tabs = line2.split('\t')

                dnds = float(tabs[7])

                dn = float(tabs[9])

                ds = float(tabs[11])
                
                pair = tuple(sorted((g1, g2)))
                print(pair)
                
                ##Add dN and dS filters here if wanted
                
                if pair not in new_pairwise_dict.keys():
                    new_pairwise_dict[pair] = [[dn], [ds], [dnds]]
                else:
                    # print("It's aliiiive!\n")
                    new_pairwise_dict[pair][0].append(dn)
                    new_pairwise_dict[pair][1].append(ds)
                    new_pairwise_dict[pair][2].append(dnds)

    return new_pairwise_dict
                
                


##Read files

locus_pairs = {}
for folder in os.listdir(os.chdir(outpath)):
    if os.path.isdir(folder) and 'A1401' in folder and '12720' not in folder:
        locus_pairs = parse_codeml_output(f'{folder}/{folder}.txt', locus_pairs)

        
##Get mean and median values and save them to files
                
with open(finalfile, 'w') as stats_file:
    stats_file.write('strain1\tstrain2\tmean_dN\tmedian_dN\tmean_dS\tmedian_dS\tmean_w\tmedian_w\n')
    for fp in sorted(locus_pairs.keys()):
        mean_dN = mean(locus_pairs[fp][0])
        median_dN = median(locus_pairs[fp][0])
        mean_dS = mean(locus_pairs[fp][1])
        median_dS = median(locus_pairs[fp][1])
        mean_w = mean(locus_pairs[fp][2])
        median_w = median(locus_pairs[fp][2])
        stats_file.write(f'{fp[0]}\t{fp[1]}\t{mean_dN:.4f}\t{median_dN:.4f}\t{mean_dS:.4f}\t')
        stats_file.write(f'{median_dS:.4f}\t{mean_w:.4f}\t{median_w:.4f}\n')
        

