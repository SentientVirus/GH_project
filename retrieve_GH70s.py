#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 17 2022

author: Julia Pedersen and Marina Mota

"""
# =============================================================================
# Script to get GH70 domain sequences and GH32 protein sequences
# =============================================================================

import os                                           #import necessary packages
from Bio import SeqIO
import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import as_fasta

direct = 'interproscan'
indir = 'gbks'
gene_types = snakemake.input

def get_domain_pos(directory, domain_annot):
    gene_names = {}
    for filename in sorted(os.listdir(directory), key=lambda v: v.upper()):
        with open(f'{directory}/{filename}') as annot:
            file = csv.reader(annot)
            for line in file:
                line = line[0].split('\t')
                if domain_annot in line:
                    if line[0] not in gene_names.keys():
                        gene_names[line[0]] = (int(line[6]), int(line[7]))
                    else:
                        gene_names[f'{line[0]}_2'] = (int(line[6]), int(line[7]))
    return gene_names

def get_GH32(directory, domain_annot):
    gene_names = []
    for filename in sorted(os.listdir(directory), key=lambda v: v.upper()):
        with open(f'{directory}/{filename}') as annot:
            file = csv.reader(annot)
            for line in file:
                line = line[0].split('\t')
                if domain_annot in line and int(line[0].split('_')[1]) > 10**4:
                    gene_names.append(line[0])
    return gene_names

def get_domains(sequence, start, end):
    domain = sequence[start:end]
    return domain



#Function to parse nucleotide and amino acid sequence from genbank:

def parse_faa_fna(indir, outdir, gene_domains, out_prefix):
    locus_tags = list(gene_domains.keys())                   #locus_tags of interest
    
    
    outfile_faa = f'{outdir}/{out_prefix}.faa' #Name outputs
    outfile_fna = f'{outdir}/{out_prefix}.fna'
    
    with open(outfile_faa, 'w') as faa, open(outfile_fna, 'w') as fna:
        faa.write('')
        fna.write('')
    
    for gbk_file in sorted(os.listdir(indir)):
        if gbk_file.endswith('_1.gbk'): #Define inputs
            with open(f'{indir}/{gbk_file}') as new_gbk:

                #parse faa and fna sequences
                with open(outfile_faa, 'a') as out_faa, open(outfile_fna, 'a') as out_fna:                #open output file where all faa and fna sequences will be saved.
                    for record in SeqIO.parse(new_gbk, 'genbank'):      #for each record in gbk file
                        if record.features:                             #if record.features
                            for feature in record.features:             #loop thorugh each feature
                                if feature.type == 'CDS':                 #if there is a CDS feature type
                                    locus_tag = feature.qualifiers['locus_tag'][0]      #retrieve locus tag as locus_tag
                                    if locus_tag in locus_tags: #if locus_tag is present in locus_tags
                                        if f'{locus_tag}_2' in locus_tags:
                                            locus_tag_list = [locus_tag, f'{locus_tag}_2']
                                        else:
                                            locus_tag_list = [locus_tag]
                                        for tag in locus_tag_list:
                                            domain = get_domains(feature.qualifiers['translation'][0], gene_domain[locus_tag][0], gene_domain[locus_tag][1])
                                            prot_record = SeqRecord(Seq(domain), tag, tag, '')
                                            out_faa.write(as_fasta(prot_record))     #write >locus_tag and amino acid sequence to output file
                                            gene_dom = get_domains(feature.location.extract(record).seq, gene_domain[locus_tag][0]*3, gene_domain[locus_tag][1]*3)
                                            gene_record = SeqRecord(Seq(gene_dom), tag, tag, '')
                                            out_fna.write(as_fasta(gene_record))    #write >locus_tag and nucleotide sequence to output file
                                            return print(f'GH70 domain with locus tag {tag} added.\n')
                                            
def parse_GH32(indir, outdir, gene_domains, out_prefix):
    locus_tags = gene_domains                   #locus_tags of interest
    
    
    outfile_faa = f'{outdir}/{out_prefix}.faa' #Name outputs
    outfile_fna = f'{outdir}/{out_prefix}.fna'
    
    with open(outfile_faa, 'w') as faa, open(outfile_fna, 'w') as fna:
        faa.write('')
        fna.write('')
    
    for gbk_file in sorted(os.listdir(indir)):
        if gbk_file.endswith('_1.gbk'): #Define inputs
            with open(f'{indir}/{gbk_file}') as new_gbk:

                #parse faa and fna sequences
                with open(outfile_faa, 'a') as out_faa, open(outfile_fna, 'a') as out_fna:                #open output file where all faa and fna sequences will be saved.
                    for record in SeqIO.parse(new_gbk, 'genbank'):      #for each record in gbk file
                        if record.features:                             #if record.features
                            for feature in record.features:             #loop thorugh each feature
                                if feature.type == 'CDS':                 #if there is a CDS feature type
                                    locus_tag = feature.qualifiers['locus_tag'][0]      #retrieve locus tag as locus_tag
                                    if locus_tag in locus_tags: #if locus_tag is present in locus_tags
                                        prot_record = SeqRecord(Seq(feature.qualifiers['translation'][0]), locus_tag, locus_tag, '')
                                        out_faa.write(as_fasta(prot_record))     #write >locus_tag and amino acid sequence to output file
                                        gene_record = SeqRecord(Seq(feature.location.extract(record).seq), locus_tag, locus_tag, '')
                                        out_fna.write(as_fasta(gene_record))    #write >locus_tag and nucleotide sequence to output file
                                        return print(f'GH32 with locus tag {locus_tag} added.\n')

for gene_type in gene_types:
    outdir = f'data/fasta/{gene_type}'
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)
        
    if gene_type == 'GH70':
        domannot = 'Glycosyl hydrolase family 70'
        gene_domain = get_domain_pos(direct, domannot)
        parse_faa_fna(indir, outdir, gene_domain, gene_type)
        
    elif gene_type == 'GH32':
        domannot = 'SSF75005'
        GH32_domain = get_GH32(direct, domannot)
        parse_GH32(indir, outdir, GH32_domain, gene_type)
        
