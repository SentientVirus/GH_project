#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:19:15 2022

@author: marina
"""

import logging, traceback
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import as_fasta

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
# Setting inputs and outputs
# =============================================================================

path2files = os.path.dirname(snakemake.input[0])

strains = [os.path.basename(snk).split('_')[0] for snk in snakemake.input]

outfile = snakemake.output[0]

#This has to be defined manually, positions of truncated GH70 genes
positions = {'A1401': (1321040, 1325531), 'G0102': (1343637, 1347958), 
             'H3B1-11M': (1324383, 1327787), 'H3B2-09X': (1379690, 1385310), 
             'Fhon2': (1396935, 1402745)}

# =============================================================================
# Function to take a substring from a string
# =============================================================================

def take_substring(strinput, start, end):
    newstr = strinput[start:end]
    return newstr


# =============================================================================
# Retrieve truncated genes in fna format
# =============================================================================

for strain in strains:
    infile = f'{strain}_genomic.fna'
    with open(f'{path2files}/{infile}') as fna: #Read input genome
        genome = SeqIO.parse(fna, 'fasta')
        for record in genome: #Loop through chromosomes/plasmids
            if 'chromosome' in record.description: #Genes are in the chromosome
                new_seq = take_substring(record.seq, positions[strain][0],
                                            positions[strain][1]) #Take sequence
                new_seq = Seq(new_seq[::-1]).complement() #Convert to right strand
                new_record = SeqRecord(Seq(new_seq), record.id, record.name, 
                                       record.description) #Convert to record
                with open(outfile, 'a+') as out_file: #Write to output
                    out_file.write(as_fasta(new_record))
                    