#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:51:51 2023

This is a script used to create tabs including information about the
glycosyl hydrolase genes under study.

@author: Marina Mota Merlo

"""

# =============================================================================
# Import packages
# =============================================================================

from Bio import GenBank
import logging, traceback
import os

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
# Define objects to store GenBank information and sorting function
# =============================================================================

def correct_loctag(locus_tag):
    corrected = locus_tag.replace('APS55_RS', 'MP2_').replace('K2W83_RS', 'DSM_').replace('LDX55', 'IBH001').replace('AKU', '')
    return corrected

class gbk_entry:
    def __init__(self, gene_name, ID, accession = 'None', strand = '-', start = 0, end = 0, representative = False):
        self.gene_name = gene_name
        self.ID = ID
        self.strain = ID.split('_')[0]
        if self.strain == 'APS55':
            self.strain = 'MP2'
        elif 'FHON' in self.strain:
            self.strain = self.strain.replace('HON', 'hon')
        elif 'K2W83' == self.strain:
            self.strain = 'DSMZ12361'
        elif 'LDX55' == self.strain:
            self.strain = 'IBH001'
        elif self.strain[0] == 'H':
            self.strain = self.strain[:4] + '-' + self.strain[4:]
        self.locus_tag = correct_loctag(self.ID)
        if self.strain != 'MP2' and self.strain != 'DSMZ12361' and self.strain != 'IBH001':
            self.locus_tag = 'AKU' + self.ID
        self.ID = ID
        self.accession = accession
        self.strand = strand
        self.start = start
        self.end = end
        self.representative = representative
    def update_repr(self, rep_list):
        if self.strain in rep_list:
            self.representative = True
    def update_pos_acc(self, location, strand, accession):
        self.start = int(location[0])
        self.end = int(location[1])
        self.strand = strand
        self.accession = accession
    def __str__(self):
        return f'<{self.ID} ({self.gene_name}) at {(self.start, self.end)} ({self.strand}) with accession {self.accession}>'

# =============================================================================
# Function to retrieve position and accession id from locus tag
# =============================================================================

def search_gbk(filename, loc_tag, path):
    gbk_file = f'{path}/{filename}'
    with open(gbk_file) as handle:
        for record in GenBank.parse(handle):
            for feature in record.features:
                if feature.key == 'CDS':
                    for qual in feature.qualifiers:
                        if 'locus_tag' in qual.key:
                            locus = qual.value
                            if loc_tag in locus:
                                if 'complement' in feature.location:
                                    strand = '-'
                                    loc = feature.location.replace('complement(', '').strip(')')
                                else:
                                    strand = '+'
                                    loc = feature.location
                                loc = loc.split('..')
                                
                                for qual2 in feature.qualifiers:
                                    if 'protein_id' in qual2.key:
                                        acc = qual2.value.strip('\"')
                                
                                return loc, strand, acc

# =============================================================================
# Read input parameters from Snakemake
# =============================================================================

GS1 = snakemake.params.GS1
GS2 = snakemake.params.GS2
GS3 = snakemake.params.GS3
GS4 = snakemake.params.GS4
BRS = snakemake.params.BRS
NGB = snakemake.params.NGB
short = snakemake.params.short
S1 = snakemake.params.S1
S2a = snakemake.params.S2a
S2b = snakemake.params.S2b
S3 = snakemake.params.S3

representatives = snakemake.params.representatives

gbk_dir = os.path.expanduser('~') + '/Akunkeei_files/gbff'
workdir = os.path.expanduser('~') + '/GH_project'

# =============================================================================
# Make strain names consistent with locus tags
# =============================================================================

repres = [rep.strip('-') for rep in representatives]

# =============================================================================
# Join separated GS2_BRS into one locus_tag, and remove extra locus tags
# =============================================================================

GS2_BRS = [GS[:-2] for GS in GS2 if '_2' in GS]
GS2 = [GS for GS in GS2 if '_2' not in GS]
BRS = [BS for BS in BRS if BS[:-2] not in GS2_BRS]


# =============================================================================
# To access variables by name
# =============================================================================

myVars = globals()

# =============================================================================
# Define list of genes to loop through
# =============================================================================

type_list = ['GS1', 'GS2', 'GS3', 'GS4', 'BRS', 'GS2_BRS', 'NGB', 'short', 'S1', 'S2a', 'S2b', 'S3']

# =============================================================================
# Save info into gbk_entry object and write to file
# =============================================================================

for GH_type in type_list:
    var = myVars[GH_type]
    out_file = f'{workdir}/tables/{GH_type}_suppl.tab'
    with open(out_file, 'w') as out_tab:
        out_tab.write('Strain\tGene_name\tLocus_tag\tLocus_id\tProtein_id\tStrand\tStart\tEnd\tRepresentative\n')
    for locus_tag in var:
        gbk = gbk_entry(GH_type, locus_tag)
        print(gbk.ID)
        gbk.update_repr(repres)
        location, strand, acc = search_gbk(f'{gbk.strain}_genomic.gbff', gbk.ID, gbk_dir)
        gbk.update_pos_acc(location, strand, acc)
        print(f'Get info for gene {gbk.ID}')
        with open(out_file, 'a') as out_tab:
            out_tab.write(f'{gbk.strain}\t{gbk.gene_name}\t{gbk.locus_tag}\t')
            out_tab.write(f'{gbk.ID}\t{gbk.accession}\t{gbk.strand}\t{gbk.start}\t')
            out_tab.write(f'{gbk.end}\t{gbk.representative}\n')
            
