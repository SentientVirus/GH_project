#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 17 2022

This is a script to retrieve the Interproscan annotations of those genes that
have a GH70 or GH32 domain in the A. kunkeei genomes.

author: Marina Mota-Merlo

"""
# =============================================================================
# Script to get GH70 domain sequences and GH32 protein sequences
# =============================================================================

import os, sys    #import necessary packages
from pathlib import Path
import csv
import logging, traceback

# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ])
                  )

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# Defining paths and loading parameters from Snakemake
# =============================================================================

direct = 'interproscan'
indir = os.path.dirname(snakemake.input[0])
outdir = os.path.dirname(snakemake.output[0])
gene_types = snakemake.output
file_list = snakemake.input
strains = [Path(file).stem for file in file_list] #Extract strain names from files
strains = [strain[:-2] for strain in strains] #Remove _1 from strain names

# =============================================================================
# Function to get GH70 domain tags from Interproscan annotations (in same folder)
# =============================================================================

def get_GH70(directory, strains, domain_annot):
    GH70_gene_names = [] #Create empty list
    for filename in sorted(strains, key=lambda v: v.upper()): #Loop through files in folder
        if filename == 'Fhon2':
            filename = 'fhon2'

        with open(f'{directory}/{filename}.tsv') as annot: #Open annotations
            file = csv.reader(annot) #Read file
            for line in file: #Loop through file
                line = line[0].split('\t') #Get list with information in the line
                if domain_annot in line: #If it is a GH70 domain
                    loctag = line[0].replace('-', '')
                    loctag = loctag.replace('fhon2', 'FHON2')
                    print(loctag)
                    if loctag not in GH70_gene_names: #Retrieve domain position
                        GH70_gene_names.append(loctag)
            
    return GH70_gene_names #Returns dictionary locus_tag: (domain start, domain end)


# =============================================================================
# Same as the function above, but gets GH32 full protein locus tags
# =============================================================================

def get_GH32(directory, strains, domain_annot):
    GH32_gene_names = []
    for filename in sorted(strains, key=lambda v: v.upper()):
        with open(f'{directory}/{filename}.tsv') as annot:
            file = csv.reader(annot)
            for line in file:
                line = line[0].split('\t')
                val = line[0].split('_')[1]
                if val.startswith('RS'):
                    val = val[2:]
                if domain_annot in line and (('RS' in line[0] and int(val) > 10**3) or ('RS' not in line[0] and int(line[0].split('_')[1]) > 10**4)): #Take only locus tags > 10000
                    loctag = line[0].replace('-', '') 
                    if loctag not in GH32_gene_names: #Retrieve domain position
                        GH32_gene_names.append(loctag)
                        
    return GH32_gene_names

# =============================================================================
# Function to get only the domain sequence from a protein
# =============================================================================

def save_annot(directory, strains, outdir, gene_names, prefix = 'GH70'):
    outfile = f'{outdir}/{prefix}_interpro.tsv'
    with open(outfile, 'w+') as handle:
        handle.write('')
        
    for filename in sorted(strains, key=lambda v: v.upper()): #Loop through files in folder
        if filename == 'Fhon2':
            filename = 'fhon2'

        with open(f'{directory}/{filename}.tsv') as annot: #Open annotations
            file = csv.reader(annot) #Read file
            for line in file: #Loop through file
                MP2_check = False
                loctag = line[0].split('\t')[0].replace('-', '')
                if loctag.startswith('H') or loctag.startswith('A') or loctag.startswith('G'):
                    full_loctag = 'AKU' + loctag
                elif 'LDX55' in loctag:
                    full_loctag = loctag
                    loctag = full_loctag.replace('LDX55', 'IBH001')
                elif 'K2W83_RS' in loctag:
                    full_loctag = loctag
                    loctag = full_loctag.replace('K2W83_RS', 'DSMZ_')
                elif loctag == 'MP2_13350':
                    print("MP2's GS1 retrieved")
                    full_loctag = 'APS55_RS03850'
                    loctag = full_loctag.replace('APS55_RS', 'MP2_')
                    MP2_check = True
                elif loctag == 'MP2_13360':
                    full_loctag = 'APS55_RS03845'
                    loctag = full_loctag.replace('APS55_RS', 'MP2_')
                    MP2_check = True
                elif loctag == 'MP2_14250':
                    full_loctag = 'APS55_RS03400'
                    loctag = full_loctag.replace('APS55_RS', 'MP2_')
                    MP2_check = True
                line_text = line[0].split('\t')[1:]
                line_text = [loctag, full_loctag] + line_text
                line_text = '\t'.join(line_text)
                if loctag in gene_names or full_loctag in gene_names or MP2_check:
                    with open(outfile, 'a') as handle:
                        handle.write(f'{line_text}\n')

# =============================================================================
# Implementation of the functions
# =============================================================================
for gene_type in gene_types: #Loop through output filenames
    gtype = Path(gene_type).stem.split('_')[0] #Get gene type: GH70 or GH32
    print(gtype)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir) #Create outdir if it does not exist
        
    if gtype == 'GH70':
        domannot = 'Glycosyl hydrolase family 70'
        GH70_list = get_GH70(direct, strains, domannot)
        save_annot(direct, strains, outdir, GH70_list, gtype)
        
    elif gtype == 'GH32':
        domannot = 'SSF75005'
        GH32_list = get_GH32(direct, strains, domannot)
        save_annot(direct, strains, outdir, GH32_list, gtype)
        
