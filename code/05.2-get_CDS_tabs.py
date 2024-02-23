# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 21:57:29 2021

This script retrieves positional information about the CDS in and around the
deleted region so that this information can be later plotted.

@author: Marina Mota Merlo
"""
import os, logging, traceback
import re

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
# Input/output definitions
# =============================================================================

#homedir = '/home/marina'
outdir = 'plots/tabfiles'
phylo_file = snakemake.input.tree
inputs = snakemake.input.gbff

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
# =============================================================================
# Class definition
# =============================================================================

class CDS:
    def __init__(self, name, start, end, strand, locus):
        self.name = name #3-letter code for the gene
        self.start = start #Start position of the protein gene
        self.end = end #End position of the gene
        self.strand = strand #Strand where the gene is (+ or -)
        self.locus = locus #Name of the file (genome)
    def __str__(self):
        return f'Protein: {self.name} in genome {self.locus}, {self.start}-{self.end}, strand {self.strand}'

# =============================================================================
# Dictionary to plot strands following the order in the phylogeny
# =============================================================================

phylogeny = {}
n = 1
with open(phylo_file) as f:
    lines = f.readlines()
    for line in lines:
        line = line.replace(' ', '')
        line = line.replace('\n', '')
        no = str(n)
        while len(no) < 3:
            no = '0' + no
        n += 1
        phylogeny[line] = no

# =============================================================================
# Create tab files with CDS positions
# =============================================================================

#omit = ['Fhon13', 'H3B2-03M-01', 'H3B2-03M-04'] #Or use shadow: minimal in Snakemake

for filename in inputs:
    if filename.endswith('.gbk'):
        f = (row for row in open(filename))
        locus = next(f)[12:23] #Take locus name (OX...)
            
        for line in f:
            if 'strain=' in line:
                strain = ''.join(e for e in line if (e.isalnum() or e == '-'))[6:]
                out_file = f'{outdir}/{phylogeny[strain]}_{strain}_gpr.tab'
                print(out_file)
                with open(out_file, 'w') as gn: #Start writing tab file
                    gn.write('name\tstart\tend\tstrand\tregion\n') #Add headers to file
            if line[5:8] == 'CDS':
                check = False
                if 'complement' in line: #If the gene is in the complementary strand
                    line = re.sub(r'[(complement)]', '', line) #Remove the word "complement"
                    strand = '-' #Save strand
                else:
                    strand = '+'
                position = line[21:-1].split('..') #Save position of CDS
                for i in range(8):
                    gene = next(f) #Go to the next line
                    if '/gene' in gene and 'gtfC' not in gene: #If the next line includes gene name
                        gene = gene[:-1].replace(' ', '') #Remove spaces
                        gene = gene.replace('/gene=', '') #Remove other text
                        gene = re.sub(r'[("")]', '', gene) #Remove ""
                        check = True
                        break
                        #print(gene)
                    if '/locus_tag' in gene: #If the next line includes gene name
                        locus_tag = gene[:-1].replace(' ', '') #Remove spaces
                        locus_tag = locus_tag.replace('/locus_tag=', '') #Remove other text
                        locus_tag = re.sub(r'[("")]', '', locus_tag) #Remove ""
                        print(locus_tag)
                if check == False:
                    gene = locus_tag
                start = position[0].replace('ji', '')
                start = start.replace('>', '')
                start = start.replace('<', '')
                end = position[1].replace(',1', '')
                end = end.replace('>', '')
                end = end.replace('<', '')
                p = CDS(gene, start, end, strand, locus) #Save info in object

                basename = filename.split('/')[1].split('_')[0]
                out_file = f'{outdir}/{phylogeny[strain]}_{basename}_gpr.tab'
                with open(out_file, 'a') as gn: #Add object info to file
                    gn.write(f'{p.name}\t{p.start}\t{p.end}\t{p.strand}\t{p.locus}\n')
