#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 17 2022

This is a script to retrieve the complete protein and nucleotide sequence of
the genes that contain GH70 domains in A. kunkeei based on previous
Interproscan annotations.

author: Julia Pedersen and Marina Mota

"""
# =============================================================================
# Script to get GH70 domain sequences and GH32 protein sequences
# =============================================================================

import os, sys    #import necessary packages
from pathlib import Path
from Bio import SeqIO
import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import as_fasta
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
indir = 'gbks'
gene_types = snakemake.output
file_list = snakemake.input
strains = [Path(file).stem for file in file_list] #Extract strain names from files
strains = [strain[:-2] for strain in strains] #Remove _1 from strain names


# =============================================================================
# Same as the function above, but gets GH70 full protein locus tags
# =============================================================================

def get_GH70(directory, strains, domain_annot):
    gene_names = []
    for filename in sorted(strains, key=lambda v: v.upper()):
        with open(f'{directory}/{filename}.tsv') as annot:
            file = csv.reader(annot)
            for line in file:
                line = line[0].split('\t')
                if domain_annot in line:
                    loctag = line[0].replace('-', '')
                    loctag = loctag.replace('fhon2', 'FHON2')
                    if filename == 'MP2':
                        if loctag == 'MP2_13350':
                            loctag = 'APS55_RS03850'
                        elif loctag == 'MP2_13360':
                            loctag = 'APS55_RS03845'
                        elif loctag == 'MP2_14250':
                            loctag = 'APS55_RS03400'
                    gene_names.append(loctag)
    return gene_names

def get_GH32(directory, strains, domain_annot):
    gene_names = []
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
                    gene_names.append(loctag)
    return gene_names
 
# =============================================================================
# Function to parse nucleotide and amino acid sequences from GenBank files:
# =============================================================================
                                           
def parse_GH(indir, outdir, gene_domains, out_prefix):
    locus_tags = gene_domains                   #locus_tags of interest
    
    
    outfile_faa = f'{outdir}/complete_{out_prefix}.faa' #Name outputs
    outfile_fna = f'{outdir}/complete_{out_prefix}.fna'
    
    with open(outfile_faa, 'w') as faa, open(outfile_fna, 'w') as fna:
        faa.write('')
        fna.write('')
    
    for gbk_file in sorted(strains, key=lambda v: v.upper()):
        print(gbk_file)
        with open(f'{indir}/{gbk_file}_1.gbk') as new_gbk:

            #parse faa and fna sequences
            with open(outfile_faa, 'a') as out_faa, open(outfile_fna, 'a') as out_fna:                #open output file where all faa and fna sequences will be saved.
                for record in SeqIO.parse(new_gbk, 'genbank'):      #for each record in gbk file
                    if record.features:                             #if record.features
                        for feature in record.features:             #loop thorugh each feature
                            if feature.type == 'CDS':                 #if there is a CDS feature type
                                if feature.qualifiers['locus_tag'][0][:3] == 'AKU':
                                    locus_tag = feature.qualifiers['locus_tag'][0][3:]  #retrieve locus tag as locus_tag
                                else:
                                    locus_tag = feature.qualifiers['locus_tag'][0]
                                print(locus_tag)
                                if locus_tag in locus_tags: #if locus_tag is present in locus_tags
                                    prot_record = SeqRecord(Seq(feature.qualifiers['translation'][0]), locus_tag, locus_tag, '')
                                    out_faa.write(as_fasta(prot_record))     #write >locus_tag and amino acid sequence to output file
                                    gene_record = SeqRecord(Seq(feature.location.extract(record).seq), locus_tag, locus_tag, '')
                                    out_fna.write(as_fasta(gene_record))    #write >locus_tag and nucleotide sequence to output file
                                    print(f'{out_prefix} with locus tag {locus_tag} added.\n')

# =============================================================================
# Implementation of the functions
# =============================================================================

for gene_type in gene_types: #Loop through output filenames
    gtype = Path(gene_type).stem.replace('complete_', '') #Get gene type: GH70
    outdir = f'data/fasta/{gtype}'
    if not os.path.exists(outdir):
        os.makedirs(outdir) #Create outdir if it does not exist
        
    if gtype == 'GH70':
        domannot = 'Glycosyl hydrolase family 70'
        GH70_prot = get_GH70(direct, strains, domannot)
        parse_GH(indir, outdir, GH70_prot, gtype)
    elif gtype == 'GH32':
        domannot = 'SSF75005'
        GH32_prot = get_GH32(direct, strains, domannot)
        parse_GH(indir, outdir, GH32_prot, gtype)
        
