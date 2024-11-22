# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 00:34:10 2021

This is a script that can be used to retrieve the sequences of the genes
that are in the deleted region but are not glycosyl hydrolases.

@author: Marina Mota Merlo
"""
import os
import logging, traceback
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
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
# Define objects to store GenBank information and sorting function
# =============================================================================

class gbk_entry:
    def __init__(self, name, locus_tag, translation, gene_seq, description = ''):
        self.name = name
        self.locus_tag = locus_tag
        self.translation = translation
        self.gene_seq = gene_seq
        self.length = len(self.translation)
        self.description = description
    def __str__(self):
        return f'<{self.name} at {self.locus_tag} with length {self.length}>'
    
def sort_key(text:str):
    return(text[0].upper(), text[1:])

# =============================================================================
# Set paths, parameters and define dictionary to store information
# =============================================================================

current_dir = '/home/marina/GH_project'
outdir = current_dir + '/' + os.path.dirname(snakemake.output[0]) #'/home/marina/GH_project/data/fasta/other_genes'
gbk_dir = os.path.dirname(snakemake.input[0])

gene_dict = {}

if not os.path.exists(outdir):
    os.makedirs(outdir)

genes = snakemake.params.gene_names
genes = [f'{gene}_I' for gene in genes]

# gtypes = ['GH70', 'GH32']

# =============================================================================
# Loop through gene type and strain
# =============================================================================

os.chdir(gbk_dir)
file_list = sorted(os.listdir(os.getcwd()), key=sort_key)
for filename in file_list:
    for gene_n in genes:
        with open(filename) as handle:
            print(filename, gene_n)
            for record in SeqIO.parse(handle, 'genbank'):
                for feature in record.features:
                    # Retrieve locus tag and sequence
                    if 'gene' in feature.qualifiers.keys() and 'translation' in feature.qualifiers.keys():
                        gene_str = feature.qualifiers['gene'][0]
                        locus_tag = feature.qualifiers['locus_tag'][0][3:]
                        # Filtering for genes that appear several times
                        if gene_str == gene_n:
                            #gene_n = gene_n.split('_')[0]
                            print(gene_n, 'check 3')
                            file_out = f'{outdir}/a_kunkeei_{gene_n[:-2]}.faa'
                            file_fna = file_out.replace('.faa', '.fna')
                            if gene_n[:-2] not in gene_dict.keys():
                                gene_dict[gene_n[:-2]] = []
                                with open(file_out, 'w') as gn, open(file_fna, 'w') as gn2:
                                    gn.write('')
                                    gn2.write('')
                            
                            # Create GenBank object
                            gbkobj = gbk_entry(gene_n[:-2], locus_tag, 
                                               feature.qualifiers['translation'][0], 
                                               feature.location.extract(record).seq,
                                               feature.qualifiers['product'][0])
                            
                            gene_dict[gene_n[:-2]].append(gbkobj)
                            
                            # Save entries in fasta format
                            with open(file_out, 'a') as outfile:
                                sequence = Seq(gbkobj.translation)
                                entry = SeqRecord(sequence, gbkobj.locus_tag, 
                                                  gbkobj.name, description = gbkobj.description)
                                outfile.write(as_fasta(entry))
                            with open(file_fna, 'a') as outfile2:
                                gene_entry = SeqRecord(gbkobj.gene_seq, gbkobj.locus_tag, 
                                                  gbkobj.name, description = gbkobj.description)
                                outfile2.write(as_fasta(gene_entry))
                            print(gbkobj)

# for gene_n in genes:
#     os.chdir(gbk_dir)
#     file_list = sorted(os.listdir(os.getcwd()), key=sort_key)
#     for filename in file_list:
#         with open(filename) as handle:
#             for record in SeqIO.parse(handle, 'genbank'):
#                 for feature in record.features:
#                     # Retrieve locus tag and sequence
#                     if 'gene' in feature.qualifiers.keys() and 'translation' in feature.qualifiers.keys():
#                         gene_str = feature.qualifiers['gene'][0]
#                         locus_tag = feature.qualifiers['locus_tag'][0][3:]
#                         # Filtering for genes that appear several times
#                         if gene_str == gene_n:
#                             gene_n = gene_n.split('_')[0]
#                             print(gene_n, 'check 3')
#                             file_out = f'{outdir}/a_kunkeei_{gene_n}.faa'
#                             file_fna = file_out.replace('.faa', '.fna')
#                             if gene_n not in gene_dict.keys():
#                                 gene_dict[gene_n] = []
#                                 with open(file_out, 'w') as gn, open(file_fna, 'w') as gn2:
#                                     gn.write('')
#                                     gn2.write('')
                            
#                             # Create GenBank object
#                             gbkobj = gbk_entry(gene_n, locus_tag, 
#                                                feature.qualifiers['translation'][0], 
#                                                feature.location.extract(record).seq,
#                                                feature.qualifiers['product'][0])
                            
#                             gene_dict[gene_n].append(gbkobj)
                            
#                             # Save entries in fasta format
#                             with open(file_out, 'a') as outfile:
#                                 sequence = Seq(gbkobj.translation)
#                                 entry = SeqRecord(sequence, gbkobj.locus_tag, 
#                                                   gbkobj.name, description = gbkobj.description)
#                                 outfile.write(as_fasta(entry))
#                             with open(file_fna, 'a') as outfile2:
#                                 gene_entry = SeqRecord(gbkobj.gene_seq, gbkobj.locus_tag, 
#                                                   gbkobj.name, description = gbkobj.description)
#                                 outfile2.write(as_fasta(gene_entry))
#                             print(gbkobj)
                                    
                                        
