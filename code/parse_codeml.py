#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Created on Tue June 21 2022

Script from https://github.com/faylward/dnds, slightly modified by Julia 
Pedersen and Marina Mota. It parses the CodeML output and calculates some
statistics for each subset of genes.

"""



import re, logging, traceback

from numpy import mean, median

# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# Load inputs from Snakefile
# =============================================================================

infile = snakemake.input['txt'] #'codeml_out/a_kunkeei_GH70/a_kunkeei_GH70.txt'
outfile1 = snakemake.output['dNdS'] #'codeml_out/a_kunkeei_GH70/new_dNdS.tsv'
outfile2 = snakemake.output['stats'] #'codeml_out/a_kunkeei_GH70/new_stats.tsv'

# =============================================================================
# Define function to parse CodeML output
# =============================================================================

def parse_codeml_output(infile, outfile1, outfile2):

    dn_list = []
    ds_list = []
    dnds_list = []

    status = 0



    outfile = open(outfile1, 'w')

    outfile.write('locus1\tlocus2\tdN\tdS\tw\n') #Add header to output file

    results = open(infile, 'r')
    
    for i in results:

        if i.startswith('pairwise comparison, codon frequencies'):

            status = 1



        if status == 1:

            if i[0].isdigit():

                line = i.rstrip()

                line2 = re.sub('\(', '', line)

                line3 = re.sub('\)', '', line2)

                spaces = line3.split(" ")

                first = spaces[1]

                second = spaces[4]


            if i.startswith('t='):

                line = i.rstrip()

                line1 = re.sub('=', '', line)

                line2 = re.sub('\s+', '\t', line1)

                tabs = line2.split('\t')

                dnds = float(tabs[7])
                dnds_list.append(dnds)

                dn = float(tabs[9])
                dn_list.append(dn)

                ds = float(tabs[11])
                ds_list.append(ds)



#                if ds < 50 and ds > 0.01 and dnds < 99:

#                    outfile.write(f'{first}\t{second}\t{dn}\t{ds}\t{dnds}\n')
                outfile.write(f'{second}\t{first}\t{dn}\t{ds}\t{dnds}\n')
                print(f'Added dN, dS and w for pair {first} vs {second}.')
    
    with open(outfile2, 'w') as outno2:
        filtered_ds = [ds for ds in ds_list if ds >= 0]
        filtered_dn = [dn_list[i] for i in range(len(ds_list)) if ds_list[i] >= 0]
        filtered_dnds = [dnds_list[j] for j in range(len(ds_list)) if ds_list[j] >= 0]                
        outno2.write(f'Mean_dN\t{mean(filtered_dn):.4f}\n')
        outno2.write(f'Median_dN\t{median(filtered_dn):.4f}\n')
        outno2.write(f'Mean_dS\t{mean(filtered_ds):.4f}\n')
        outno2.write(f'Median_dS\t{median(filtered_ds):.4f}\n')
        outno2.write(f'Mean_w\t{mean(filtered_dnds):.4f}\n')
        outno2.write(f'Median_w\t{median(filtered_dnds):.4f}\n')

# =============================================================================
# Run function in list comprehension
# =============================================================================

[parse_codeml_output(infile[i], outfile1[i], outfile2[i]) for i in range(len(infile))]

