# -*- coding: utf-8 -*-
"""
This is a script that counts how many instances of each gene are in each genome.
Needed to make the CDS plot of the dynamic region in ete3.

@author: Marina Mota-Merlo
"""
import os, logging, traceback
import csv
from Bio import SeqIO

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
# Define inputs from Snakemake
# =============================================================================

output_file = snakemake.output.file
outdir = os.path.expanduser('~') + '/GH_project/plots/counts'
if not os.path.exists(outdir):
    os.makedirs(outdir) #Create output directory if it doesn't exist

#Read list of locus tags per type
GS1 = snakemake.params.GS1
GS2 = snakemake.params.GS2
GS3 = snakemake.params.GS3
GS4 = snakemake.params.GS4
BRS = snakemake.params.BRS
short = snakemake.params.short
NGB = snakemake.params.NGB
S1 = snakemake.params.S1
S2a = snakemake.params.S2a
S2b = snakemake.params.S2b
S2 = S2a + S2b
S3 = snakemake.params.S3
phylo = snakemake.input.tree
GH70 = snakemake.input.GH70

#Create a dictionary to calculate the value of the counts for GH70 genes
value_dict = {}
with open(GH70) as lengths:
    for record in SeqIO.parse(lengths, 'fasta'):
        if len(record.seq) >= 750: #Add only 0.5 to the count for proteins that are not full-length
            value_dict[record.id] = 1
        elif record.id not in short:
            value_dict[record.id] = 0.5


all_GH70 = GS1 + GS2 + GS3 + GS4 + BRS + short + NGB

all_GH32 = S1 + S2a + S2b + S3

count = [0, 0, 0, 0, 0, 0, 0, 0]
with open(output_file, 'w') as tab:
    tab.write('strain\tGS1\tGS2\tBRS\tNGB\tshort_gene\ttotal_GH70\tS1\tS2\tS3\ttotal_GH32\tloci\n')

with open(phylo) as phy: #Order by phylogeny
    p = (tip for tip in phy)
    for strain in p:
        loci = ''
        strain = strain.strip()
        if len(strain) > 0:
            strain = strain.replace('-', '')
            count[0] = sum(value_dict[s] for s in GS1 if (strain.upper() in s) or ('APS55' in s and strain == 'MP2') or ('LDX55' in s and strain == 'IBH001') or ('K2W83' in s and strain == 'DSMZ12361')) #Account for locus tags with a different format for MP2
            count[1] = sum(value_dict[s] for s in GS2 if (strain.upper() in s) or ('APS55' in s and strain == 'MP2') or ('LDX55' in s and strain == 'IBH001') or ('K2W83' in s and strain == 'DSMZ12361'))
            count[2] = sum(value_dict[s] for s in BRS if (strain.upper() in s) or ('APS55' in s and strain == 'MP2') or ('LDX55' in s and strain == 'IBH001') or ('K2W83' in s and strain == 'DSMZ12361'))
            count[3] = sum([1 for s in NGB if (strain.upper() in s) or ('APS55' in s and strain == 'MP2') or ('LDX55' in s and strain == 'IBH001') or ('K2W83' in s and strain == 'DSMZ12361')])
            count[4] = sum([1 for s in short if (strain.upper() in s) or ('APS55' in s and strain == 'MP2') or ('LDX55' in s and strain == 'IBH001') or ('K2W83' in s and strain == 'DSMZ12361')])
            count = [int(num) for num in count] #Convert floats to ints
            total = sum(count)
            for GH70 in all_GH70: #Add loci names
                if strain.upper() in GH70 or ('APS55' in GH70 and strain == 'MP2') or ('LDX55' in GH70 and strain == 'IBH001') or ('K2W83' in GH70 and strain == 'DSMZ12361'):
                    loci += f'{GH70}, '
            count[5] = sum([1 for s in S1 if strain in s or ('K2W83' in s and strain == 'DSMZ12361')])
            count[6] = sum([1 for s in S2 if strain in s or ('K2W83' in s and strain == 'DSMZ12361')])
            count[7] = sum([1 for s in S3 if strain in s or ('K2W83' in s and strain == 'DSMZ12361')])
            total2 = sum(count[5:8])
            #Add also cell wall-binding domains
            '''with open(f'../Alignment/interproscan/{strain}.tsv') as interpro:
                annot = csv.reader(interpro, delimiter ='\t')
                cb = 0
                cwb = 0
                for row in annot:
                    if row[3] == 'Pfam':
                        if 'Choline-binding' in row[5]:
                            cb += 1
                        elif 'cell wall binding' in row[5]:
                            cwb += 1
            count[5] = cb
            count[6] = cwb
            count[7] = cb + cwb'''
            #Write info to file
            with open(output_file, 'a') as tab:
                tab.write(f'{strain}\t{count[0]}\t{count[1]}\t{count[2]}\t'+
                          f'{count[3]}\t{count[4]}\t{total}\t{count[5]}\t'+
                          f'{count[6]}\t{count[7]}\t{total2}\t{loci[:-2]}\n')
        count = [0, 0, 0, 0, 0, 0, 0, 0]
