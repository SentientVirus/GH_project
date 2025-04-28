#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:51:51 2023

This is a script used to create tabs including information about the
glycosyl hydrolase genes under study.

@author: Marina Mota-Merlo

"""

# =============================================================================
# 0. Import packages
# =============================================================================

from Bio import GenBank
import logging, traceback
import os
import csv

# =============================================================================
# 0. Logging
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
# 1. Define objects to store CDS information and function to fix locus tags
# =============================================================================

def correct_loctag(locus_tag):
    '''Function to convert the locus tags into the version that is displayed
    in most figures and tables.
    Input: locus_tag (the locus tag to be modified)
    Output: corrected (the modified locus tag)'''
    corrected = locus_tag.replace('APS55_RS', 'MP2_').replace('K2W83_RS', 'DSM_').replace('LDX55', 'IBH001').replace('AKU', '')
    return corrected

class gbk_entry:
    def __init__(self, gene_name, ID, accession = 'None', strand = '-', start = 0, end = 0, GH_len = '', representative = False):
        #Assign properties to the object
        self.gene_name = gene_name
        self.ID = ID
        self.strain = ID.split('_')[0]
        if self.strain == 'APS55': #Conditional statement to assign the right strain name
            self.strain = 'MP2'
        elif 'FHON' in self.strain:
            self.strain = self.strain.replace('HON', 'hon')
        elif 'K2W83' == self.strain:
            self.strain = 'DSMZ12361'
        elif 'LDX55' == self.strain:
            self.strain = 'IBH001'
        elif self.strain[0] == 'H':
            self.strain = self.strain[:4] + '-' + self.strain[4:]
        self.locus_tag = correct_loctag(self.ID) #Apply function to correct the locus tag
        if self.strain != 'MP2' and self.strain != 'DSMZ12361' and self.strain != 'IBH001': #If the strain is from Dyrhage et al. (2022)
            self.locus_tag = 'AKU' + self.ID #Add AKU at the beginning of the locus tag
        self.ID = ID
        self.accession = accession
        self.strand = strand
        self.start = start
        self.end = end
        self.length = end - start
        self.GH_len = GH_len
        self.representative = representative
    def update_repr(self, rep_list): #Function to determine if the gene is in a representative strain
        if self.strain in rep_list:
            self.representative = True
    def update_pos_acc(self, location, strand, accession): #Function to update the location of the CDS
        self.start = int(location[0]) #Start position
        self.end = int(location[1]) #End position
        self.strand = strand #Strand
        self.accession = accession #NCBI accession
        self.length = self.end-self.start #CDS length
    def update_domains(self, domain_list): #Function to add the no. of domains present in the protein from list
        self.SP = domain_list[0] #The first element of the list is the presence/absence of the signal peptide
        self.GH = domain_list[1] #The second element is the no. of GH domains
        if domain_list[3] > 0: #If the fourth element is different to 0
            self.GB = f'{domain_list[2]}+{domain_list[3]}' #Add it together with the third element as glucan-bindind domains
        else: self.GB = domain_list[2] #Otherwise, the third element is the no. of glucan-binding domains
    def __str__(self): #Function to determine what print(object) returns
        return f'<{self.ID} ({self.gene_name}) at {(self.start, self.end)} ({self.strand}) with accession {self.accession}>'

# =============================================================================
# 2. Define the function to retrieve position and accession id from locus tag
# =============================================================================

def search_gbk(filename, loc_tag, path):
    '''Function that reads a file and searches for a locus tag to retrieve
    positional information.
    Inputs:
        - filename: Name of the file to be searched
        - loc_tag: Locus tag to be searched
        - path: Path to the file
    Outputs:
        - loc: Start and end position of the locus tag.
        - strand: Strand where the locus tag is present (+/-)
        - acc: Protein accession number of the CDS'''
    gbk_file = f'{path}/{filename}' #Full path to the GenBank file
    with open(gbk_file) as handle: #Open the GenBank file
        for record in GenBank.parse(handle): #Loop through records in the file
            for feature in record.features: #Loop through features in the records
                if feature.key == 'CDS': #If the feature is a CDS
                    for qual in feature.qualifiers: #Loop through qualifiers in the feature
                        if 'locus_tag' in qual.key: #If locus_tag is in the qualifier name
                            locus = qual.value #Retrieve the value of the qualifier
                            if loc_tag in locus: #If the locus tag to be searched is in the retrieved locus tag
                                if 'complement' in feature.location: #If complement is in the string with the location
                                    strand = '-' #Set the strand to the reverse
                                    loc = feature.location.replace('complement(', '').strip(')') #Remove complement() from the string
                                else: #Otherwise
                                    strand = '+' #Set the strand to forward
                                    loc = feature.location #Don't modify the location text
                                loc = loc.split('..') #Split the location text into a list where the first element is the start and the second, the end
                                
                                for qual2 in feature.qualifiers: #Loop through qualifiers again
                                    if 'protein_id' in qual2.key: #If protein_id is in the name of the qualifier
                                        acc = qual2.value.strip('\"') #Set the accession to the value of the qualifier, removing special characters
                                
                                return loc, strand, acc #Return positional information

# =============================================================================
# Read input parameters from Snakemake
# =============================================================================
#Lists with the locus tags of each gene type
GS1 = snakemake.params.GS1
GS2 = snakemake.params.GS2
GS3 = snakemake.params.GS3
GS4 = snakemake.params.GS4
BRS = snakemake.params.BRS
NGB = snakemake.params.NGB
short = snakemake.params.short
S1 = snakemake.params.S1 + ['G0417_12900', 'G0417_12910']
S2a = snakemake.params.S2a
S2b = snakemake.params.S2b
S3 = snakemake.params.S3
outgroup_file = snakemake.params.outgroup_file

representatives = snakemake.params.representatives

myVars = globals() #Variable that allows to access other variables by name

gbk_dir = os.path.expanduser('~') + '/Akunkeei_files/gbff' #Directory with GenBank files
workdir = os.path.expanduser('~') + '/GH_project' #Working directory
interpro_dir = f'{workdir}/interproscan' #Directory with Interproscan annotations

type_list = ['GS1', 'GS2', 'GS3', 'GS4', 'BRS', 'GS2_BRS', 'NGB', 'short', 
             'S1', 'S2a', 'S2b', 'S3'] #List of genes to loop through

repres = [rep.strip('-') for rep in representatives] #Make strain names consistent with locus tags

#Join separated GS2_BRS into one locus_tag, and remove extra locus tags
GS2_BRS = [GS[:-2] for GS in GS2 if '_2' in GS]
GS2 = [GS for GS in GS2 if '_2' not in GS]
BRS = [BS for BS in BRS if BS[:-2] not in GS2_BRS]


# =============================================================================
# Save info into gbk_entry object and write to file
# =============================================================================

for GH_type in type_list: #Loop through gene types
    var = myVars[GH_type] #Retrieve the variable with the right gene type
    out_file = f'{workdir}/tables/{GH_type}_suppl.tab' #Define the path to the output file
    with open(out_file, 'w') as out_tab: #Create or empty the output file
        #Add headers
        out_tab.write('Strain\tGene_name\tLocus_tag\tLocus_id\t')
        out_tab.write('Protein_id\tStrand\tStart\tEnd\tLength\t')
        out_tab.write('Signal_peptide\tGH_domains\tGH_domain_length\t')
        out_tab.write('GB_domains\tRepresentative\n')
        
    for locus_tag in var: #Loop through locus tags in the list
        gbk = gbk_entry(GH_type, locus_tag) #Create an empty object
        gbk.update_repr(repres) #Add representative status
        location, strand, acc = search_gbk(f'{gbk.strain}_genomic.gbff', gbk.ID, gbk_dir) #Get the positional information of the CDS
        gbk.update_pos_acc(location, strand, acc) #Add the information to the object
        print(f'Get info for gene {gbk.ID}') #Print progress
        
        domain_dict = {} #Create an empty dictionary to store domains
        with open(f'{interpro_dir}/{gbk.strain}.tsv') as interpro: #Open the Interproscan file
            annot = csv.reader(interpro, delimiter ='\t') #Read the file
            
            for row in annot: #Loop through rows in the file
                to_compare = row[0].upper().replace('-', '') #Modify locus tags to remove - symbols

                #Manually provide the right locus tags for MP2
                if to_compare == 'MP2_13360':
                    to_compare = 'APS55_RS03845'
                elif to_compare == 'MP2_13350':
                    to_compare = 'APS55_RS03850'
                elif to_compare == 'MP2_14250':
                    to_compare = 'APS55_RS03400'
                if to_compare == gbk.ID: #If the current locus tag corresponds with the object ID
                    if to_compare not in domain_dict.keys(): #If the locus tag is not in the dictionary yet
                        domain_dict[to_compare] = [False, 0, 0, 0] #Create a list indicating total absence of domains
                    if row[3] == 'Pfam': #If the annotation is a Pfam annotation
                        if 'Choline-binding' in row[5]: #If the domain is annotated as choline binding
                            domain_dict[to_compare][2] += 1 #Add +1 to the glucan-binding domain count
                        elif 'cell wall binding' in row[5]: #If the domain is annotated as cell wall-binding
                            domain_dict[to_compare][3] += 1 #Add +1 to the cell wall-binding domain count
                        elif 'Glycosyl hydrolase family 70' in row[5]: #If the domain is GH70
                            domain_dict[to_compare][1] += 1 #Add +1 to the GH domain count
                            gbk.GH_len += f'{int(row[7])-int(row[6])}+' #Assign GH domain length
                    elif row[3] == 'SMART' and 'glyco_32' in row[5]: #If the annotation comes from SMART and the domain is GH32
                        domain_dict[to_compare][1] += 1 #Add +1 to the GH32 domain count
                        gbk.GH_len = f'{int(row[7])-int(row[6])}+' #Assign GH domain length
                    if 'signal peptide' in row[5]: #If the domain is a signal peptide
                        domain_dict[to_compare][0] = True #Set the signal peptide presence to True
        gbk.update_domains(domain_dict[gbk.ID]) #Add the domains to the object
        gbk.GH_len = gbk.GH_len[:-1]
        if gbk.GH_len == '':
            gbk.GH_len = '0'
        
        if gbk.strain in ['DSMZ12361', 'IBH001', 'MP2']: #If the strain is not from Dyrhage et al. (2022)
            #Change the locus tag to the ID, and vice-versa
            locus_tag = gbk.ID
            locus_ID = gbk.locus_tag
            gbk.locus_tag = locus_tag
            gbk.ID = locus_ID
                        
        with open(out_file, 'a') as out_tab: #Open output file
            #Write the information stored in the object to the file
            out_tab.write(f'{gbk.strain}\t{gbk.gene_name}\t{gbk.locus_tag}\t')
            out_tab.write(f'{gbk.ID}\t{gbk.accession}\t{gbk.strand}\t{gbk.start}\t')
            out_tab.write(f'{gbk.end}\t{gbk.length}\t{gbk.SP}\t{gbk.GH}\t')
            out_tab.write(f'{gbk.GH_len}\t{gbk.GB}\t{gbk.representative}\n')
