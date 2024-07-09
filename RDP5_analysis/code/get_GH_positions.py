#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 13:31:36 2022

@author: marina
"""

import logging, traceback, sys
from Bio import SeqIO
import os
import subprocess

# =============================================================================
# 0. Logging
# =============================================================================
log = os.path.expanduser('~') + '/GH_project/RDP5_analysis/logs/get_alignments.log'

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
# Paths
gbk_dir = os.path.expanduser('~') + '/Akunkeei_files/gbff'
fna_dir = os.path.expanduser('~') + '/Akunkeei_files/fna'
infile_suffix = 'genomic.gbff'
infile_suffix2 = 'genomic.fna'
pos_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab'
out_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/fna'
start_padding = 25000
end_padding = 35000
threads = 16 #Threads to run mafft

# Strains and comparisons of interest
strain_groups = {'H3B2-03M': 0, 'H4B4-02J': 0, 'H4B5-03X': 0, 'H1B3-02M': 1,
                 'H3B2-09X': 1, 'H4B5-05J': 1, 'MP2': 2, 'H3B2-02X': 2}

group_names = {0: 'root_GS1_S2-3_subset', 1: 'GS1-2_BRS', 2: 'only_GS1+GS2'}

GS1_loctag = {'H3B2-03M': 'AKUH3B203M_12480', 'H4B4-02J': 'AKUH4B402J_12600',
              'H4B5-03X': 'AKUH4B503X_12670', 'H1B3-02M': 'AKUH1B302M_12880',
              'H4B5-05J': 'AKUH4B505J_12880', 'H3B2-09X': 'AKUH3B209X_13340',
              'H3B2-02X': 'AKUH3B202X_12850', 'MP2': 'APS55_RS03850'}

if not os.path.exists(pos_dir):
    os.makedirs(pos_dir)
    
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
for group in group_names.values():
    with open(f'{pos_dir}/{group}_positions.tab', 'w+') as pos_file:
        pos_file.write('locus_tag\tstrain\tgene_name\tstart\tend\tstrand\n')
    with open(f'{out_dir}/{group}_seqs.fasta', 'w+') as out_file:
        out_file.write('')


dict_pos = {}
genes_in_segment = {}
for strain in strain_groups.keys():
    outfile = f'{group_names[strain_groups[strain]]}_seqs.fasta'
    with open(f'{gbk_dir}/{strain}_{infile_suffix}') as gbk_file:
        gbk = SeqIO.parse(gbk_file, 'genbank')
        start = 0
        end = 0
        strand = 0
        genes_in_segment[strain] = []
        for genomic_element in gbk:
            for CDS in genomic_element.features:
                if start == 0 and 'locus_tag' in CDS.qualifiers and CDS.qualifiers['locus_tag'][0] == GS1_loctag[strain] and CDS.type == 'CDS':
                     start = int(CDS.location.start)
                     end = int(CDS.location.end)
                     strand = CDS.location.strand
                     dict_pos[strain] = [start - start_padding, end + end_padding, strand]
                     if strain == 'MP2':
                         dict_pos[strain] = [start - end_padding, end + start_padding, strand*-1]
            if start != 0:
                for gene in genomic_element.features:
                    if 'locus_tag' in gene.qualifiers and int(gene.location.start) >= dict_pos[strain][0] and int(gene.location.end) <= dict_pos[strain][1] and gene.type == 'CDS':
                        gene_tag = gene.qualifiers['locus_tag'][0]
                        if 'gene' in gene.qualifiers:
                            gene_name = gene.qualifiers['gene'][0]
                        else: gene_name = '-'
                        gene_start = int(gene.location.start) - dict_pos[strain][0]
                        gene_end = int(gene.location.end) - dict_pos[strain][0]
                        gene_strand = gene.location.strand
                        if strain == 'MP2':
                            full_length = dict_pos[strain][1] - dict_pos[strain][0]
                            gene_strand*=-1
                            gene_start = dict_pos[strain][1] - int(gene.location.end)
                            gene_end = dict_pos[strain][1] - int(gene.location.start)
                            
                            
                        with open(f'{pos_dir}/{group_names[strain_groups[strain]]}_positions.tab', 'a+') as pos_file:
                            pos_file.write(f'{gene_tag}\t{strain}\t{gene_name}\t{gene_start}\t{gene_end}\t{gene_strand}\n')
                        # info_dict = {gene_tag: [gene_name, gene_start, gene_end, gene_strand]}
                        # genes_in_segment[strain].append(info_dict)
                    
                
    with open(f'{out_dir}/{outfile}', 'a+') as out_file:
        with open(f'{fna_dir}/{strain}_{infile_suffix2}') as fna:
            fna_read = SeqIO.parse(fna, 'fasta')
            for record in fna_read:
                segment = record.seq[dict_pos[strain][0]:dict_pos[strain][1]]
                record.id = strain
                record.description = ''
                record.seq = segment
                if strain == 'MP2':
                    record.seq = record.seq.reverse_complement()
                SeqIO.write(record, out_file, 'fasta')
                break
                
                  
for group in group_names.values():
    raw_data = f'{out_dir}/{group}_seqs.fasta'
    alignment = raw_data.replace('fasta', 'mafft.fasta')
    subprocess.run(f'mafft-linsi --thread {threads} {raw_data} > {alignment} 2> {log}', shell = True)            
        
    

# start = 1270000 #Chosen so that reference gene is included and distance is 40000
# end = 1301000 #Original 1310000
# comp_dir = '../tabfiles'
# genome_dir = '../data/fna'
# outdir = '../RDP5/unaligned'
# gbk_dir = '../annotated_gbks'
# GH_dir = '../GH_sequences'
# GH_types = ['GH70', 'GH32']
# GH_dir2 = '../codeml_analysis/codeml_input'
# GH_dir3 = '../codeml_analysis/inputs'
# outdir = 'GH_positions'

# GH_loci = []
# with open(f'{GH_dir2}/region_{GH_types[0]}_prot.fasta', 'r') as GH70, open(f'{GH_dir3}/a_kunkeei_all{GH_types[1]}_prot.fasta', 'r') as GH32:
#     for record in SeqIO.parse(GH70, 'fasta'):
#         GH_loci.append(record.id)
#     for record in SeqIO.parse(GH32, 'fasta'):
#         GH_loci.append(record.id)


# with open(f'{GH_dir}/same_domain_sequences.txt', 'r') as GH70, open(f'{GH_dir}/same_{GH_types[1]}_proteins_sequences.txt', 'r') as GH32:
#     for line in GH70:
#         if not line.startswith('The'):
#             GH70_list = line.split(' ')
#             for GH70_locus in GH70_list:
#                 GH70_locus = GH70_locus.rstrip()
#                 if GH70_locus not in GH_loci and len(GH70_locus) > 0:
#                     if GH70_list[0].rstrip() in GH_loci:
#                         GH_loci.append(GH70_locus)
#     for line in GH32:
#         if not line.startswith('The'):
#             GH32_list = line.split(' ')
#             for GH32_locus in GH32_list:
#                 GH32_locus = GH32_locus.rstrip()
#                 if GH32_locus not in GH_loci and len(GH32_locus) > 0:
#                     if GH32_list[0].rstrip() in GH_loci:
#                         GH_loci.append(GH32_locus)

# GH_loci.sort()

# pos_dict = {}
# prot_dict = {}
# #In this section, I fill a dictionary with start/end positions of the deleted region
# for leaf in list(strain_groups.keys()):
#     for filename in sorted(os.listdir(comp_dir)):
#         if leaf in filename and filename.endswith('gpr.tab'):
#             check = False
#             with open(f'{comp_dir}/{filename}', 'r') as gfile:
#                 genes = csv.reader(gfile, delimiter='\t')
#                 for gene in genes:
#                     if 'name' not in gene: #Skip first line
#                         if gene[0] == 'ohrR' and not check: #Gene used to center the plotted fragments in the right position
#                             start = int(gene[1])+5000
#                             if leaf == 'MP2' or leaf == 'H3B2-03M':
#                                 start += 1400
#                             elif leaf == 'H1B1-05A':
#                                 start += 1300
#                             elif leaf == 'H1B3-02M':
#                                 start += 1450
#                             end = start + 30500 #Original 40000
#                             check = True
#                     if check:
#                         pos_dict[leaf] = (start, end)
                        
#     filename = f'{gbk_dir}/{leaf}.gbk'
#     f = (row for row in open(filename))

#     for line in f: #Loop through lines in file
#         if line[22:27] == 'locus': #Get gene names
#             gene = line[:-1].replace(' ', '') #Remove spaces
#             gene = gene.replace('/locus_tag=', '') #Remove other text
#             gene = re.sub(r'[("")]', '', gene) #Remove ""
#             if any(x in line for x in GH_loci): #If it is the gene(s) of interest
#                 line2 = next(f)
#                 prot_pos = line2
#                 prot_pos = prot_pos.split('..')
#                 prot_pos = [int(re.sub('\D', '', prot_pos[i])) for i in range(len(prot_pos))]
#                 if pos_dict[leaf][1] >= prot_pos[0] >= pos_dict[leaf][0]:
#                     prot_dict[gene] = (prot_pos[0] - pos_dict[leaf][0],
#                                         prot_pos[1] - pos_dict[leaf][0])
                    
# #with open(outfile, 'w') as out:
# #    out.write('locus_tag\tstart\tend\n')
                    
# gap_prot_dict = {}
# #Before I finish this I need to add the gaps to the original position(s).
# for group in group_names.values():
#     with open(f'{outdir}/{group}_{outdir}.tab', 'w') as out:
#         out.write('locus_tag\tstart\tend\n')
#     with open(f'aligned/{group}_regions.mafft.fasta', 'r') as alignment:
#         for record in SeqIO.parse(alignment, 'fasta'):
#             for locus in prot_dict.keys():
#                 count1 = 0
#                 count2 = 0
#                 if record.id in locus:
#                     start, end = prot_dict[locus]
#                     #print(group, locus)
                    
#                 ##OBS! This won't work, I need to take into account identical proteins
#                     for i in range(len(record.seq)):
#                         if record.seq[i] != '-':
#                             count2 += 1
#                         count1 += 1
#                         if count2 == start:
#                             gapstart = count1
#                         elif count2 == end:
#                             gapend = count1
#                             gap_prot_dict[locus] = (gapstart, gapend)
#                             break
                    
#                     with open(f'{outdir}/{group}_{outdir}.tab', 'a') as out:
#                         out.write(f'{locus}\t{gap_prot_dict[locus][0]}\t{gap_prot_dict[locus][1]}\n')
