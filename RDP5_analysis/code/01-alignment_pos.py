#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 13:31:36 2022

This script concatenates two regions of the genome of A. kunkeei and,
keeping all the gaps, generates an alignment of the region in the
representative subset of 38 strains or the subset of 34 strains without 
transposons.

@author: Marina Mota-Merlo
"""
# =============================================================================
# 0. Import required modules
# =============================================================================

import logging, traceback, sys
from Bio import SeqIO
import os
import subprocess

# =============================================================================
# 0. Logging
# =============================================================================

log = os.path.expanduser('~') + '/GH_project/RDP5_analysis/logs/get_A-D_alignments.log'

if os.path.exists(log):
    os.remove(log)

if not os.path.exists(os.path.dirname(log)):
    os.makedirs(os.path.dirname(log))

logger = logging.getLogger()

logging.basicConfig(filename = log, level = logging.INFO,
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

sys.stdout = open(log, 'a')

# =============================================================================
# 0. Input definition
# =============================================================================

#Dictionary with the locus tags in each strain for each of the genes included in the aligned region
#Basal strains A1001 and A1404 are not considered
bcrA_loctag = {'A0901': 'AKUA0901_13170', 'A1003': 'AKUA1003_12430',
               'A1202': 'AKUA1202_13410', 'A1401': 'AKUA1401_12650', #'A1404': 'AKUA1404_13340',
               'A1805': 'AKUA1805_12670', 'DSMZ12361': 'K2W83_RS06120',
               'Fhon2': 'AKUFHON2_13430', 'G0101': 'AKUG0101_12670',
               'G0403': 'AKUG0403_13000', 'H1B1-04J': 'AKUH1B104J_12890',
               'H1B1-05A': 'AKUH1B105A_12180', 'H1B3-02M': 'AKUH1B302M_12740', 
               'H3B1-01A': 'AKUH3B101A_13140', 'H3B1-04J': 'AKUH3B104J_12880', 
               'H3B1-04X': 'AKUH3B104X_13100', 'H3B2-02X': 'AKUH3B202X_12750', 
               'H3B2-03J': 'AKUH3B203J_13270', 'H3B2-03M': 'AKUH3B203M_12350', 
               'H3B2-06M': 'AKUH3B206M_12730', 'H3B2-09X': 'AKUH3B209X_13230',
               'H4B1-11J': 'AKUH4B111J_13460', 'H4B2-02J': 'AKUH4B202J_12780',
               'H4B2-04J': 'AKUH4B204J_13230', 'H4B2-05J': 'AKUH4B205J_12870',
               'H4B2-06J': 'AKUH4B206J_13300', 'H4B2-11M': 'AKUH4B211M_12880', 
               'H4B4-02J': 'AKUH4B402J_12490', 'H4B4-05J': 'AKUH4B405J_13250',
               'H4B4-06M': 'AKUH4B406M_13350', 'H4B4-12M': 'AKUH4B412M_13120', 
               'H4B5-01J': 'AKUH4B501J_12770', 'H4B5-03X': 'AKUH4B503X_12560', 
               'H4B5-04J': 'AKUH4B504J_13360', 'H4B5-05J': 'AKUH4B505J_12770',
               'IBH001': 'LDX55_06275', 'MP2': 'APS55_RS03910'}

ydiL_loctag = {'A0901': 'AKUA0901_13260', 'A1003': 'AKUA1003_12530', 
               'A1202': 'AKUA1202_13510', 'A1401': 'AKUA1401_12740', #'A1404': 'AKUA1404_13400',
               'A1805': 'AKUA1805_12760', 'DSMZ12361': 'K2W83_RS06170',
               'Fhon2': 'AKUFHON2_13520', 'G0101': 'AKUG0101_12770',
               'G0403': 'AKUG0403_13090', 'H1B1-04J': 'AKUH1B104J_12980',
               'H1B1-05A': 'AKUH1B105A_12290', 'H1B3-02M': 'AKUH1B302M_12870', 
               'H3B1-01A': 'AKUH3B101A_13230', 'H3B1-04J': 'AKUH3B104J_12980',
               'H3B1-04X': 'AKUH3B104X_13190', 'H3B2-02X': 'AKUH3B202X_12840', 
               'H3B2-03J': 'AKUH3B203J_13360', 'H3B2-03M': 'AKUH3B203M_12460', 
               'H3B2-06M': 'AKUH3B206M_12820', 'H3B2-09X': 'AKUH3B209X_13320',
               'H4B1-11J': 'AKUH4B111J_13550', 'H4B2-02J': 'AKUH4B202J_12870',
               'H4B2-04J': 'AKUH4B204J_13320', 'H4B2-05J': 'AKUH4B205J_12960',
               'H4B2-06J': 'AKUH4B206J_13390', 'H4B2-11M': 'AKUH4B211M_12970',
               'H4B4-02J': 'AKUH4B402J_12580', 'H4B4-05J': 'AKUH4B405J_13340',
               'H4B4-06M': 'AKUH4B406M_13440', 'H4B4-12M': 'AKUH4B412M_13210', 
               'H4B5-01J': 'AKUH4B501J_12860', 'H4B5-03X': 'AKUH4B503X_12650', 
               'H4B5-04J': 'AKUH4B504J_13450', 'H4B5-05J': 'AKUH4B505J_12860',
               'IBH001': 'LDX55_06320', 'MP2': 'APS55_RS03855'}

wzx_loctag =  {'A0901': 'AKUA0901_13400', 'A1003': 'AKUA1003_12580',
               'A1202': 'AKUA1202_13560', 'A1401': 'AKUA1401_12810',
               'A1805': 'AKUA1805_12860', 'DSMZ12361': 'K2W83_RS06205',
               'Fhon2': 'AKUFHON2_13630', 'G0101': 'AKUG0101_12870',
               'G0403': 'AKUG0403_13160', 'H1B1-04J': 'AKUH1B104J_13100',
               'H1B1-05A': 'AKUH1B105A_12340', 'H1B3-02M': 'AKUH1B302M_12940', 
               'H3B1-01A': 'AKUH3B101A_13300', 'H3B1-04J': 'AKUH3B104J_13060',
               'H3B1-04X': 'AKUH3B104X_13290', 'H3B2-02X': 'AKUH3B202X_12900', 
               'H3B2-03J': 'AKUH3B203J_13430', 'H3B2-03M': 'AKUH3B203M_12530', 
               'H3B2-06M': 'AKUH3B206M_12880', 'H3B2-09X': 'AKUH3B209X_13430',
               'H4B1-11J': 'AKUH4B111J_13610', 'H4B2-02J': 'AKUH4B202J_12930',
               'H4B2-04J': 'AKUH4B204J_13410', 'H4B2-05J': 'AKUH4B205J_13060',
               'H4B2-06J': 'AKUH4B206J_13530', 'H4B2-11M': 'AKUH4B211M_13070',
               'H4B4-02J': 'AKUH4B402J_12670', 'H4B4-05J': 'AKUH4B405J_13400', 
               'H4B4-06M': 'AKUH4B406M_13500', 'H4B4-12M': 'AKUH4B412M_13340', 
               'H4B5-01J': 'AKUH4B501J_12960', 'H4B5-03X': 'AKUH4B503X_12770', 
               'H4B5-04J': 'AKUH4B504J_13520', 'H4B5-05J': 'AKUH4B505J_12980',
               'IBH001': 'LDX55_06355', 'MP2': 'APS55_RS03825'}

tagU_plus2_loctags = {'A0901': 'AKUA0901_13530', 'A1003': 'AKUA1003_12720',
               'A1202': 'AKUA1202_13690', 'A1401': 'AKUA1401_12940',
               'A1805': 'AKUA1805_12990', 'DSMZ12361': 'K2W83_RS06270',
               'Fhon2': 'AKUFHON2_13760', 'G0101': 'AKUG0101_13000',
               'G0403': 'AKUG0403_13290', 'H1B1-04J': 'AKUH1B104J_13230',
               'H1B1-05A': 'AKUH1B105A_12470', 'H1B3-02M': 'AKUH1B302M_13070', 
               'H3B1-01A': 'AKUH3B101A_13430', 'H3B1-04J': 'AKUH3B104J_13190',
               'H3B1-04X': 'AKUH3B104X_13420', 'H3B2-02X': 'AKUH3B202X_13030', 
               'H3B2-03J': 'AKUH3B203J_13560', 'H3B2-03M': 'AKUH3B203M_12690', 
               'H3B2-06M': 'AKUH3B206M_13010', 'H3B2-09X': 'AKUH3B209X_13560',
               'H4B1-11J': 'AKUH4B111J_13740', 'H4B2-02J': 'AKUH4B202J_13060',
               'H4B2-04J': 'AKUH4B204J_13540', 'H4B2-05J': 'AKUH4B205J_13190',
               'H4B2-06J': 'AKUH4B206J_13660', 'H4B2-11M': 'AKUH4B211M_13200',
               'H4B4-02J': 'AKUH4B402J_12800', 'H4B4-05J': 'AKUH4B405J_13530',
               'H4B4-06M': 'AKUH4B406M_13630', 'H4B4-12M': 'AKUH4B412M_13470', 
               'H4B5-01J': 'AKUH4B501J_13090', 'H4B5-03X': 'AKUH4B503X_12900', 
               'H4B5-04J': 'AKUH4B504J_13650', 'H4B5-05J': 'AKUH4B505J_13110',
               'IBH001': 'LDX55_06420', 'MP2': 'APS55_RS03760'}

# =============================================================================
# 0. Create a gene object class
# =============================================================================
# Maybe not needed, I can use CDS objects
class geneObj:
    '''A class to store CDS information, can be initialized as an empty object
    with default values'''
    def __init__(self, seq = '', start = 0, end = 0, strand = '+', name = 'unk.', annot = 'hypothetical protein', locus_tag = '', strain = '', organism = 'Apilactobacillus kunkeei'):
        self.seq = seq
        self.start = start
        self.end = end
        self.strand = strand
        self.annot = annot
        self.name = name
        self.locus_tag = locus_tag
        self.strain = strain
        self.organism = organism
    def __str__(self): # What print(class) returns
        return f'<{self.locus_tag} ({self.name}) from {self.organism} {self.strain} at position ({self.start}, {self.end}), strand {self.strand}>'
    
# =============================================================================
# 1. Create paths and input files
# =============================================================================

#Paths
gbk_dir = os.path.expanduser('~') + '/Akunkeei_files/gbff'
fna_dir = os.path.expanduser('~') + '/Akunkeei_files/fna'
infile_suffix = 'genomic.gbff'
infile_suffix2 = 'genomic.fna'
pos_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab'
out_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/fna'
threads = 48 #Threads to run mafft
padding = 0 #5000
prefixes = ['A', 'B', 'C', 'D']

if not os.path.exists(pos_dir):
    os.makedirs(pos_dir)
    
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

for prefix in prefixes:
    tabdir = f'{pos_dir}/{prefix}'
    tab_out = f'{tabdir}/{prefix}.tab'
    aln_dir = f'{out_dir}/{prefix}'
    aln_file = f'{aln_dir}/{prefix}.mafft.fasta'
    
    if not os.path.exists(tabdir):
        os.makedirs(tabdir)
        
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    
    with open(tab_out, 'w+') as gpos:
        gpos.write('locus_tag\tstrain\tgene_name\tstart\tend\tstrand\n')
        
    for i in [1,2]:
        with open(f'{aln_dir}/{prefix}{i}.fasta', 'w+') as out_file:
            out_file.write('')
    
    full_record = {}
    dict_pos = {}
    all_genes = {}
    genes_in_segment = {}
    loop =  list(bcrA_loctag.keys())

    for strain in loop:
        segment1 = [0, 0]
        segment2 = [0, 0]

        with open(f'{gbk_dir}/{strain}_{infile_suffix}') as gbk_file:
            gbk = SeqIO.parse(gbk_file, 'genbank')
            genes_in_segment[strain] = []
            all_genes[strain] = []
            for genomic_element in gbk:
                for CDS in genomic_element.features:
                    if 'locus_tag' in CDS.qualifiers and CDS.type == 'CDS' and 'translation' in CDS.qualifiers:
                        gene_obj = geneObj(seq = CDS.qualifiers['translation'][0], 
                                        start = int(CDS.location.start), 
                                        end = int(CDS.location.end),
                                        strand = CDS.location.strand,
                                        annot = CDS.qualifiers['product'][0],
                                        locus_tag = CDS.qualifiers['locus_tag'][0],
                                        strain = strain)
                        
                        if strain == 'MP2':
                            gene_obj.strand *= -1
                        
                        if 'gene' in CDS.qualifiers:
                            gene_obj.name = CDS.qualifiers['gene'][0]
                        
                        all_genes[strain].append(gene_obj)
                        
                        if gene_obj.locus_tag == bcrA_loctag[strain]:
                            if strain != 'MP2':
                                segment1[0] = gene_obj.start - padding
                            else:
                                segment2[1] = gene_obj.end + padding
                        elif gene_obj.locus_tag == ydiL_loctag[strain]:
                            if strain != 'MP2':
                                segment1[1] = gene_obj.end
                            else:
                                segment2[0] = gene_obj.start
                        elif gene_obj.locus_tag == wzx_loctag[strain]:
                            if strain != 'MP2':
                                segment2[0] = gene_obj.start
                            else:
                                segment1[1] = gene_obj.end
                        elif gene_obj.locus_tag == tagU_plus2_loctags[strain]:
                            if strain != 'MP2':
                                segment2[1] = gene_obj.end + padding
                            else:
                                segment1[0] = gene_obj.start - padding
                          
        dict_pos[strain] = [segment1, segment2]
            
        for gene in all_genes[strain]:
            check1 = False
            check2 = False
            padding = 0
            if prefix == 'A':
                padding = 1
            if dict_pos[strain][0][0] <= gene.start < dict_pos[strain][0][1]:
                if strain == 'MP2':
                    gene.start = gene.start - dict_pos[strain][0][0] + 1
                    gene.end = gene.end - dict_pos[strain][0][0]
                else:
                    segment_length = dict_pos[strain][1][1] - dict_pos[strain][1][0] + 1 + padding
                    gene_start = gene.start
                    gene.start = dict_pos[strain][0][1] - gene.end + segment_length + 1
                    gene.end = dict_pos[strain][0][1] - gene_start + segment_length 
                
                genes_in_segment[strain].append(gene)
                check1 = True
        
            elif dict_pos[strain][1][0] <= gene.start < dict_pos[strain][1][1]:
                segment_length = dict_pos[strain][0][1] - dict_pos[strain][0][0]
                if strain == 'MP2':
                    gene.start = gene.start - dict_pos[strain][1][0] + segment_length + 1
                    gene.end = gene.end - dict_pos[strain][1][0] + segment_length
                else:
                    gene_start = gene.start
                    gene.start = dict_pos[strain][1][1] - gene.end + 1
                    gene.end = dict_pos[strain][1][1] - gene_start
                genes_in_segment[strain].append(gene)
                check2 = True

            if check1 or check2:
                with open(tab_out, 'a+') as pos_file:
                    pos_file.write(f'{gene.locus_tag}\t{strain}\t{gene.name}\t{gene.start}\t{gene.end}\t{gene.strand}\n')
                        
                        
        with open(f'{fna_dir}/{strain}_{infile_suffix2}') as fna:
            fna_read = SeqIO.parse(fna, 'fasta')
            for record in fna_read:
                if strain != 'MP2':
                    segment = record.seq[dict_pos[strain][0][0]:dict_pos[strain][0][1]] + '!'
                    middle_point = len(segment)-1
                    segment += record.seq[dict_pos[strain][1][0]:dict_pos[strain][1][1]]
                else: # Almost fixed, there is now a single 1 kb gap in MP2, possibly transposons?
                    genome_len = len(record.seq)
                    start_s1 = genome_len - dict_pos[strain][1][1]
                    end_s1 = genome_len - dict_pos[strain][1][0] + 1
                    start_s2 = genome_len - dict_pos[strain][0][1]
                    end_s2 = genome_len - dict_pos[strain][0][0] + 1
                    segment = record.seq[start_s1:end_s1] + '!'
                    middle_point = len(segment)-1
                    segment += record.seq[start_s2:end_s2]
                    
                record.id = strain
                record.description = ''
                record.seq = segment
                
                # if prefix == 'B' or prefix == 'D':
                #     record.seq = record.seq.replace('!', '')

                full_record[strain] = record.seq.split('!')
                    
                ind = 3
                for seq in full_record[strain]:
                    ind -= 1
                    record.seq = seq.reverse_complement()
                    record.id = f'{strain}_{ind}'
                    with open(f'{aln_dir}/{prefix}{ind}.fasta', 'a+') as global_out:
                        if prefix != 'C' and prefix != 'D':
                            print(f'Saving record {record.id} to file {aln_dir}/{prefix}{ind}.fasta')
                            SeqIO.write(record, global_out, 'fasta')
                        elif strain not in ['H1B1-05A', 'H1B3-02M', 'H3B2-03M', 'MP2']:
                            with open(f'{aln_dir}/{prefix}{ind}.fasta', 'a+') as global_out:
                                print(f'Saving record {record.id} to file {aln_dir}/{prefix}{ind}.fasta')
                                SeqIO.write(record, global_out, 'fasta')
                                    
                break
           
    for ind in [1, 2]:
        raw_data = f'{aln_dir}/{prefix}{ind}.fasta'
        alignment = raw_data.replace('fasta', 'mafft.fasta')
        subprocess.run(f'mafft --auto --thread {threads} {raw_data} > {alignment} 2>> {log}', shell = True)

    with open(f'{aln_dir}/{prefix}1.mafft.fasta') as seq1, open(f'{aln_dir}/{prefix}2.mafft.fasta') as seq2:
        records1 = list(SeqIO.parse(seq1, 'fasta'))
        records2 = list(SeqIO.parse(seq2, 'fasta'))
        with open(aln_file, 'w+') as last_outfile:
            last_outfile.write('')
        for j in range(0, len(records1)):
            if records1[j].id.split('_')[0] == records2[j].id.split('_')[0] and records1[j].seq != '' and records2[j].seq != '':
                if prefix == 'A' or prefix == 'C':
                    new_record = records1[j] + '!'
                else: new_record = records1[j]
                new_record.seq += records2[j].seq
                new_record.id = new_record.id.replace('_1', '')
                new_record.description = ''
                with open(aln_file, 'a+') as last_outfile:
                    SeqIO.write(new_record, last_outfile, 'fasta')
                    
