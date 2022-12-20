# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 00:34:10 2021

@author: usuario
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
    def __init__(self, name, subset, locus_tag, translation, gene_seq, description = ''):
        self.name = name
        self.subset = subset
        self.locus_tag = locus_tag
        self.translation = translation
        self.gene_seq = gene_seq
        self.length = len(self.translation)
        self.description = description
    def __str__(self):
        return f'<{self.name} at {self.locus_tag} with length {self.length}, {self.subset} subset>'
    
def sort_key(text:str):
    return(text[0].upper(), text[1:])

# =============================================================================
# Set paths, parameters and define dictionary to store information
# =============================================================================

current_dir = '/home/marina/GH_project'
outdir = current_dir + '/' + os.path.dirname(snakemake.output[0]) #'/home/marina/GH_project/data/fasta/other_genes'
gbk_dir = os.path.dirname(snakemake.input[0]) #'/home/marina/Akunkeei_files/gbff'

gene_dict = {}

if not os.path.exists(outdir):
    os.makedirs(outdir)

genes = snakemake.params.gene_names
genes = list(map(lambda x: x.replace('yhdG_2','yhdG'), genes))
# genes = ['ohrR', 'yifK', 'yhdG', 'GH39', 'nox', 'hpcG', 'oppA', 'sbnD', 'epsE', 
#          'epsL', 'ywqE', 'ywqD', 'ywqC', 'tagU'] #Run without GH39

gtypes = snakemake.params.gtypes
# gtypes = ['GH70', 'GH32']

# =============================================================================
# Loop through gene type and strain
# =============================================================================

for gtype in gtypes:
    if gtype == 'GH70':
        GH70_subset = snakemake.params.GH70_subset
        # GH70_subset = ['G0403', 'H1B3-02M', 'H3B2-03J', 'H3B1-04J', 'H3B1-04X', 
        #               'H4B2-02J', 'H4B2-04J', 'H4B5-04J', 'H4B5-05J']
    else:
        GH70_subset = snakemake.params.GH32_subset
        # GH70_subset = ['A1001', 'A1805', 'H1B1-04J', 'H3B2-03M',
        #               'H4B4-02J', 'H4B4-12M', 'H4B5-03X']
    for gene_n in genes:
        os.chdir(gbk_dir)
        file_list = sorted(os.listdir(os.getcwd()), key=sort_key)
        for filename in file_list:
            sname = filename.split('_')[0]
            if filename.endswith('.gbff') and sname in GH70_subset:
                with open(filename) as handle:
                    for record in SeqIO.parse(handle, 'genbank'):
                        for feature in record.features:
                            # Retrieve locus tag and sequence
                            if 'gene' in feature.qualifiers.keys() and 'translation' in feature.qualifiers.keys():
                                gene_str = feature.qualifiers['gene'][0] 
                                locus_tag = feature.qualifiers['locus_tag'][0][3:]
                                # Filtering for genes that appear several times
                                if (gene_str == gene_n and (int(locus_tag.split('_')[1]) > 12000)) or ('mhpD' in gene_str and gene_n == 'hpcG'):
                                    if not (gtype == 'GH70' and gene_n == 'yhdG' and (int(locus_tag.split('_')[1]) < 12200)):
                                        # Create empty output files and account for genes that appear several times
                                        if (gene_n, gtype) not in gene_dict.keys():
                                            gene_dict[(gene_n, gtype)] = []
                                            file_out = f'{outdir}/a_kunkeei_{gene_n}_{gtype}.faa'
                                            file_fna = file_out.replace('.faa', '.fna')
                                            if gene_n != 'yhdG':
                                                with open(file_out, 'w') as gn:
                                                    gn.write('')
                                        elif (sname.replace('-', '') in [gbk.locus_tag.split('_')[0] for gbk in gene_dict[(gene_n, gtype)]]):
                                            gene_n = f'{gene_n}_2'
                                            if (gene_n, gtype) not in gene_dict.keys() and gene_n != 'sbnD_2' and gene_n != 'GH39_2':
                                                gene_dict[(gene_n, gtype)] = []
                                                file_out = f'{outdir}/a_kunkeei_{gene_n}_{gtype}.faa' #Add .fna files to code
                                                file_fna = file_out.replace('.faa', '.fna')
                                                with open(file_out, 'w') as gn, open(file_fna, 'w') as gn2:
                                                    gn.write('')
                                                    gn2.write('')
                                        
                                        # Create GenBank object
                                        gbkobj = gbk_entry(gene_n, gtype, locus_tag, 
                                                           feature.qualifiers['translation'][0], 
                                                           feature.location.extract(record).seq,
                                                           feature.qualifiers['product'][0])
                                        
                                        # Filter out variants of the gene outside target region, save info to dictionary
                                        if gbkobj.name != 'sbnD_2' and gbkobj.name != 'GH39_2':
                                            gene_dict[(gene_n, gtype)].append(gbkobj)
                                        gene_n = gene_n.split('_')[0]
                                        
                                        # Save entries in fasta format
                                        file_out = f'{outdir}/a_kunkeei_{gbkobj.name}_{gtype}.faa'
                                        file_fna = file_out.replace('.faa', '.fna')
                                        if (gbkobj.name != 'sbnD_2') and (gbkobj.name != 'yhdG') and ('GH39' not in gbkobj.name or (gbkobj.name == 'GH39' and gbkobj.length > 500)):
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
                                        