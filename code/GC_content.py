#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 16:51:25 2024

Script to calculate the GC content at the third codon position of genes,
applied to glycosyl hydrolases

@author: Marina Mota Merlo
"""

# Import modules
import os
import logging, traceback
from Bio.SeqUtils import GC123 #Calculates GC content at each codon position
from Bio import SeqIO
#from pyfields import fields, init_fields

# =============================================================================
# 0. Logging
# =============================================================================

# logfile = 'logs/GC_content.log' #Change to snakemake.log[0]
# logging.basicConfig(filename = logfile, level = logging.INFO,
#                     format = '%(asctime)s %(message)s',
#                     datefmt = '%Y-%m-%d %H:%M:%S')

# def handle_exception(exc_type, exc_value, exc_traceback):
#     if issubclass(exc_type, KeyboardInterrupt):
#         sys.__excepthook__(exc_type, exc_value, exc_traceback)
#         return

#     logger.error(''.join(["Uncaught exception: ",
#                           *traceback.format_exception(exc_type, exc_value, exc_traceback)
#                           ]))

# sys.excepthook = handle_exception

# sys.stdout = open(logfile, 'a')


# =============================================================================
# 1. Declare inputs
# =============================================================================
inpath = '/home/marina/Akunkeei_files/fna'
infolder = 'data/fasta/'
GH_files = [f'../{infolder}/GH70/complete_GH70_all.fna', 
            f'../{infolder}/GH70/GH70_all.fna',
            f'../{infolder}/GH32/GH32_all.fna'] + [f'{inpath}/{f}' for f in os.listdir(inpath) if os.path.isfile(os.path.join(inpath, f))]

outfile = '../GC_content/GH_GC.tab'
outdir = os.path.dirname(outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir)

type_dict = {'A1001_12310': 'GS1', 'A1003_12540': 'GS1', 'A1202_13520': 'GS1', 
 'A1401_12750': 'GS1', 'A1805_12820': 'GS1', 'A2002_12300': 'GS1', 
 'A2003_12280': 'GS1', 'A2101_12830': 'GS1', 'A2103_12520': 'GS1', 
 'FHON2_13540': 'GS1', 'G0101_12800': 'GS1', 'G0102_12710': 'GS1', 
 'G0102_12720': 'GS1', 'G0103_12710': 'GS1', 'G0401_12720': 'GS1', 
 'G0402_12760': 'GS1', 'G0403_13100': 'GS1', 'G0404_12740': 'GS1', 
 'G0405_12720': 'GS1', 'G0406_13110': 'GS1', 'G0407_12710': 'GS1', 
 'G0408_12740': 'GS1', 'G0410_13100': 'GS1', 'G0412_13100': 'GS1', 
 'G0414_12740': 'GS1', 'G0415_12740': 'GS1', 'G0417_12930': 'GS1', 
 'G0420_13320': 'GS1', 'G0601_12750': 'GS1', 'G0602_12700': 'GS1', 
 'G0702_12720': 'GS1', 'G0801_12740': 'GS1', 'G0802_12690': 'GS1', 
 'G0803_12710': 'GS1', 'G0804_12720': 'GS1', 'H1B104J_13010': 'GS1', 
 'H1B105A_12300': 'GS1', 'H1B302M_12880': 'GS1', 'H2B105J_12840': 'GS1', 
 'H2B105J_12850': 'GS1', 'H3B101A_13240': 'GS1', 'H3B101J_12980': 'GS1', 
 'H3B101X_12830': 'GS1', 'H3B101X_12840': 'GS1', 'H3B102A_13470': 'GS1', 
 'H3B102X_13340': 'GS1', 'H3B102X_13350': 'GS1', 'H3B103J_12780': 'GS1', 
 'H3B103M_13200': 'GS1', 'H3B103X_12820': 'GS1', 'H3B103X_12830': 'GS1', 
 'H3B104J_12990': 'GS1', 'H3B104X_13200': 'GS1', 'H3B107A_13230': 'GS1', 
 'H3B109M_12560': 'GS1', 'H3B110M_12750': 'GS1', 'H3B111A_13170': 'GS1', 
 'H3B111M_12560': 'GS1', 'H3B202M_12850': 'GS1', 'H3B202M_12860': 'GS1', 
 'H3B202X_12850': 'GS1', 'H3B203J_13370': 'GS1', 'H3B203M_12480': 'GS1', 
 'H3B204J_13240': 'GS1', 'H3B205J_13230': 'GS1', 'H3B206M_12830': 'GS1', 
 'H3B206M_12840': 'GS1', 'H3B207X_12780': 'GS1', 'H3B208X_13260': 'GS1', 
 'H3B209X_13340': 'GS1', 'H4B101A_12820': 'GS1', 'H4B101A_12830': 'GS1', 
 'H4B102A_13620': 'GS1', 'H4B102A_13630': 'GS1', 'H4B103J_12820': 'GS1', 
 'H4B103J_12830': 'GS1', 'H4B104A_13560': 'GS1', 'H4B104A_13570': 'GS1', 
 'H4B111J_13560': 'GS1', 'H4B111J_13570': 'GS1', 'H4B114J_12820': 'GS1', 
 'H4B114J_12830': 'GS1', 'H4B116J_12780': 'GS1', 'H4B116J_12790': 'GS1', 
 'H4B202J_12880': 'GS1', 'H4B203M_12880': 'GS1', 'H4B203M_12890': 'GS1', 
 'H4B204J_13330': 'GS1', 'H4B205J_12990': 'GS1', 'H4B210M_13560': 'GS1', 
 'H4B210M_13570': 'GS1', 'H4B211M_13000': 'GS1', 'H4B303J_13010': 'GS1', 
 'H4B303J_13020': 'GS1', 'H4B402J_12600': 'GS1', 'H4B403J_12780': 'GS1', 
 'H4B403J_12790': 'GS1', 'H4B404J_12780': 'GS1', 'H4B404J_12790': 'GS1', 
 'H4B405J_13350': 'GS1', 'H4B405J_13360': 'GS1', 'H4B406M_13450': 'GS1', 
 'H4B406M_13460': 'GS1', 'H4B410M_13340': 'GS1', 'H4B410M_13350': 'GS1', 
 'H4B411M_12830': 'GS1', 'H4B411M_12840': 'GS1', 'H4B412M_13240': 'GS1', 
 'H4B501J_12890': 'GS1', 'H4B502X_13260': 'GS1', 'H4B503X_12670': 'GS1', 
 'H4B504J_13460': 'GS1', 'H4B505J_12880': 'GS1', 'H4B507J_12770': 'GS1', 
 'H4B507X_13710': 'GS1', 'H4B508X_13230': 'GS1', 'APS55_RS03850': 'GS1', 
 'LDX55_06325': 'GS1', 'K2W83_RS06180': 'GS1', 'FHON2_13560': 'GS2', 
 'FHON2_13570': 'GS2', 'G0403_13120': 'GS2', 'G0406_13130': 'GS2', 
 'G0420_13340': 'GS2', 'H1B302M_12900': 'GS2', 'H3B101A_13260': 'GS2', 
 'H3B101J_13010_2': 'GS2', 'H3B102A_13490': 'GS2', 'H3B103M_13220': 'GS2', 
 'H3B104J_13020_2': 'GS2', 'H3B104X_13220': 'GS2', 'H3B107A_13250': 'GS2', 
 'H3B109M_12580': 'GS2', 'H3B111A_13190': 'GS2', 'H3B111M_12590': 'GS2', 
 'H3B202X_12860': 'GS2', 'H3B203J_13390': 'GS2', 'H3B204J_13260': 'GS2', 
 'H3B205J_13250': 'GS2', 'H3B208X_13280': 'GS2', 'H3B209X_13360': 'GS2', 
 'H3B209X_13370': 'GS2', 'H4B202J_12890_2': 'GS2', 'H4B204J_13350_2': 'GS2', 
 'H4B502X_13280': 'GS2', 'H4B504J_13480': 'GS2', 'H4B505J_12900': 'GS2', 
 'H4B507X_13730': 'GS2', 'H4B508X_13250': 'GS2', 'APS55_RS03845': 'GS2', 
 'LDX55_06335_2': 'GS2', 'K2W83_RS06185': 'GS2', 'A1401_12760': 'BRS',
 'FHON2_13550': 'BRS', 'G0403_13110': 'BRS', 'G0406_13120': 'BRS', 
 'G0420_13330': 'BRS', 'H1B302M_12890': 'BRS', 'H3B101A_13250': 'BRS', 
 'H3B101J_12990': 'BRS', 'H3B101J_13000': 'BRS', 'H3B101J_13010': 'BRS', 
 'H3B102A_13480': 'BRS', 'H3B103M_13210': 'BRS', 'H3B104J_13000': 'BRS', 
 'H3B104J_13010': 'BRS', 'H3B104J_13020': 'BRS', 'H3B104X_13210': 'BRS', 
 'H3B107A_13240': 'BRS', 'H3B109M_12570': 'BRS', 'H3B111A_13180': 'BRS', 
 'H3B111M_12570': 'BRS', 'H3B111M_12580': 'BRS', 'H3B203J_13380': 'BRS', 
 'H3B204J_13250': 'BRS', 'H3B205J_13240': 'BRS', 'H3B208X_13270': 'BRS', 
 'H3B209X_13350': 'BRS', 'H4B202J_12890': 'BRS', 'H4B204J_13340': 'BRS', 
 'H4B204J_13350': 'BRS', 'H4B502X_13270': 'BRS', 'H4B504J_13470': 'BRS', 
 'H4B505J_12890': 'BRS', 'H4B507X_13720': 'BRS', 'H4B508X_13240': 'BRS', 
 'LDX55_06330': 'BRS', 'LDX55_06335': 'BRS', 'K2W83_RS06185_2': 'BRS', 
 'A0901_14250': 'partial', 'A1001_13210': 'partial', 'A1002_14250': 'partial', 
 'A1003_13400': 'partial', 'A1201_14240': 'partial', 'A1202_14500': 'partial', 
 'A1401_13710': 'partial', 'A1404_14290': 'partial', 'A1802_14260': 'partial', 
 'A1803_14250': 'partial', 'A1805_13670': 'partial', 'A2001_14380': 'partial', 
 'A2002_13200': 'partial', 'A2003_13180': 'partial', 'A2101_13680': 'partial', 
 'A2102_13600': 'partial', 'A2103_13380': 'partial', 'FHON2_14510': 'partial', 
 'G0101_13740': 'partial', 'G0102_13660': 'partial', 'G0103_13650': 'partial', 
 'G0401_13660': 'partial', 'G0402_13700': 'partial', 'G0403_13980': 'partial', 
 'G0404_13680': 'partial', 'G0405_13700': 'partial', 'G0406_13990': 'partial', 
 'G0407_13680': 'partial', 'G0408_13680': 'partial', 'G0410_14060': 'partial', 
 'G0412_14060': 'partial', 'G0414_13680': 'partial', 'G0415_13680': 'partial', 
 'G0417_13870': 'partial', 'G0420_14200': 'partial', 'G0601_13690': 'partial', 
 'G0602_13640': 'partial', 'G0702_13770': 'partial', 'G0801_13680': 'partial', 
 'G0802_13630': 'partial', 'G0803_13650': 'partial', 'G0804_13660': 'partial', 
 'H1B104J_13940': 'partial', 'H1B105A_13110': 'partial', 'H1B302M_13810': 'partial', 
 'H2B105J_13740': 'partial', 'H3B101A_14180': 'partial', 'H3B101J_13900': 'partial', 
 'H3B101X_13730': 'partial', 'H3B102A_14410': 'partial', 
 'H3B102X_14140': 'partial', 'H3B103J_13710': 'partial', 'H3B103M_14110': 'partial', 
 'H3B103X_13720': 'partial', 'H3B104J_13910': 'partial', 'H3B104X_14110': 'partial', 
 'H3B107A_14170': 'partial', 'H3B109M_13490': 'partial', 'H3B110M_13680': 'partial', 
 'H3B111A_14080': 'partial', 'H3B111M_13500': 'partial', 'H3B202M_13750': 'partial', 
 'H3B202X_13780': 'partial', 'H3B203J_14310': 'partial', 'H3B203M_13420': 'partial', 
 'H3B204J_14170': 'partial', 'H3B204M_13620': 'partial', 'H3B205J_14170': 'partial', 
 'H3B206M_13730': 'partial', 'H3B207X_13710': 'partial', 'H3B208X_14200': 'partial', 
 'H3B209X_14340': 'partial', 'H4B101A_13720': 'partial', 'H4B102A_14430': 'partial', 
 'H4B103J_13720': 'partial', 'H4B104A_14460': 'partial', 'H4B111J_14460': 'partial', 
 'H4B114J_13720': 'partial', 'H4B116J_13680': 'partial', 'H4B202J_13700': 'partial',
 'H4B203M_13780': 'partial', 'H4B204J_14250': 'partial', 'H4B205J_13920': 'partial',
 'H4B206J_14400': 'partial', 'H4B210M_14460': 'partial', 'H4B211M_13940': 'partial',
 'H4B303J_13910': 'partial', 'H4B402J_13530': 'partial', 'H4B403J_13680': 'partial',
 'H4B404J_13680': 'partial', 'H4B405J_14150': 'partial', 'H4B406M_14380': 'partial',
 'H4B410M_14140': 'partial', 'H4B411M_13710': 'partial', 'H4B412M_14170': 'partial',
 'H4B501J_13790': 'partial', 'H4B502X_14190': 'partial', 'H4B503X_13610': 'partial', 
 'H4B504J_14270': 'partial', 'H4B505J_13820': 'partial', 'H4B507J_13700': 'partial', 
 'H4B507X_14660': 'partial', 'H4B508X_14160': 'partial', 'APS55_RS03400': 'partial',
 'LDX55_06780': 'partial', 'K2W83_RS06625': 'partial', 'A0901_05380': 'NGB', 
 'A1002_05380': 'NGB', 'A1003_04750': 'NGB', 'A1201_05380': 'NGB', 
 'A1202_05530': 'NGB', 'A1401_04720': 'NGB', 'A1802_05380': 'NGB', 
 'A1803_05380': 'NGB', 'A1805_04920': 'NGB', 'A2001_05400': 'NGB', 
 'A2101_04930': 'NGB', 'A2102_04710': 'NGB', 'A2103_04760': 'NGB', 
 'FHON2_04830': 'NGB', 'G0101_04800': 'NGB', 'G0102_04760': 'NGB', 
 'G0103_04750': 'NGB', 'G0401_04750': 'NGB', 'G0402_04760': 'NGB', 
 'G0404_04770': 'NGB', 'G0405_04750': 'NGB', 'G0407_04750': 'NGB', 
 'G0408_04780': 'NGB', 'G0410_04950': 'NGB', 'G0412_04950': 'NGB', 
 'G0414_04750': 'NGB', 'G0415_04750': 'NGB', 'G0417_04850': 'NGB', 
 'G0601_04760': 'NGB', 'G0602_04740': 'NGB', 'G0702_04760': 'NGB', 
 'G0801_04780': 'NGB', 'G0802_04740': 'NGB', 'G0803_04750': 'NGB', 
 'G0804_04760': 'NGB', 'H1B105A_04750': 'NGB', 'H1B302M_04960': 'NGB', 
 'H2B105J_05060': 'NGB', 'H3B101A_04720': 'NGB', 'H3B101X_05080': 'NGB', 
 'H3B102A_04730': 'NGB', 'H3B102X_04940': 'NGB', 'H3B103J_04760': 'NGB', 
 'H3B103M_04730': 'NGB', 'H3B103X_05060': 'NGB', 'H3B104X_04750': 'NGB', 
 'H3B107A_04720': 'NGB', 'H3B109M_04710': 'NGB', 'H3B110M_04750': 'NGB', 
 'H3B111A_04730': 'NGB', 'H3B111M_04690': 'NGB', 'H3B202M_05080': 'NGB', 
 'H3B202X_04770': 'NGB', 'H3B203J_04720': 'NGB', 'H3B203M_04710': 'NGB', 
 'H3B204J_04730': 'NGB', 'H3B205J_04710': 'NGB', 'H3B206M_05060': 'NGB', 
 'H3B207X_04750': 'NGB', 'H3B208X_04740': 'NGB', 'H3B209X_04660': 'NGB', 
 'H4B101A_04950': 'NGB', 'H4B102A_04970': 'NGB', 'H4B103J_04940': 'NGB', 
 'H4B104A_04920': 'NGB', 'H4B111J_04920': 'NGB', 'H4B114J_05070': 'NGB', 
 'H4B116J_04920': 'NGB', 'H4B202J_04590': 'NGB', 'H4B203M_04920': 'NGB', 
 'H4B204J_04710': 'NGB', 'H4B205J_04990': 'NGB', 'H4B210M_05070': 'NGB', 
 'H4B211M_04810': 'NGB', 'H4B303J_05100': 'NGB', 'H4B403J_04940': 'NGB', 
 'H4B404J_04920': 'NGB', 'H4B405J_04950': 'NGB', 'H4B406M_05270': 'NGB', 
 'H4B410M_04950': 'NGB', 'H4B411M_04940': 'NGB', 'H4B502X_04740': 'NGB', 
 'H4B503X_04680': 'NGB', 'H4B504J_05510': 'NGB', 'H4B505J_04980': 'NGB', 
 'H4B507J_04760': 'NGB', 'H4B507X_04770': 'NGB', 'H4B508X_04720': 'NGB', 
 'K2W83_RS02365': 'NGB', 'G0101_12790': 'S1', 'G0102_12700': 'S1', 
 'G0103_12700': 'S1', 'G0401_12710': 'S1', 'G0402_12750': 'S1', 
 'G0404_12730': 'S1', 'G0405_12710': 'S1', 'G0407_12700': 'S1', 
 'G0408_12730': 'S1', 'G0410_13090': 'S1', 'G0412_13090': 'S1', 
 'G0414_12730': 'S1', 'G0415_12730': 'S1', 'G0417_12920': 'S1', 
 'G0601_12740': 'S1', 'G0602_12690': 'S1', 'G0702_12710': 'S1', 
 'G0801_12730': 'S1', 'G0802_12680': 'S1', 'G0803_12700': 'S1', 
 'G0804_12710': 'S1', 'H3B103J_12770': 'S1', 'H3B110M_12740': 'S1', 
 'H3B207X_12770': 'S1', 'H4B205J_12980': 'S1', 'H4B211M_12990': 'S1', 
 'H4B501J_12880': 'S1', 'H4B507J_12760': 'S1', 'K2W83_RS06175': 'S1', 
 'A0901_13270': 'S2a', 'A1001_12300': 'S2a', 'A1002_13270': 'S2a', 
 'A1201_13260': 'S2a', 'A1802_13280': 'S2a', 'A1803_13270': 'S2a',
 'A1805_12810': 'S2a', 'A2001_13390': 'S2a', 'A2002_12290': 'S2a',
 'A2003_12270': 'S2a', 'A2101_12820': 'S2a', 'A2102_12620': 'S2a', 
 'H1B104J_13000': 'S2a', 'H3B203M_12470': 'S2a', 'H3B204M_12620': 'S2a', 
 'H4B206J_13400': 'S2a', 'H4B402J_12590': 'S2a', 'H4B412M_13230': 'S2a',
 'H4B503X_12660': 'S2a', 'A1805_12800': 'S2b', 'A2101_12810': 'S2b', 
 'H1B104J_12990': 'S2b', 'H4B412M_13220': 'S2b', 'A0901_13360': 'S3',
 'A1001_12360': 'S3', 'A1002_13360': 'S3', 'A1201_13350': 'S3', 
 'A1404_13450': 'S3', 'A1802_13370': 'S3', 'A1803_13360': 'S3', 
 'A2001_13480': 'S3', 'A2002_12350': 'S3', 'A2003_12330': 'S3', 
 'A2102_12710': 'S3', 'H3B203M_12520': 'S3', 'H3B204M_12710': 'S3',
 'H4B206J_13490': 'S3', 'H4B402J_12660': 'S3', 'H4B412M_13310': 'S3', 
 'H4B503X_12760': 'S3'}





#genomes = [f'{inpath}/{f}' for f in os.listdir(inpath) if os.path.isfile(os.path.join(inpath, f))]

# =============================================================================
# 2. Create class GC entry
# =============================================================================
class GC_entry:
    def __init__(self, seqid, seq, description, strain = '', gene_type = ''):
        self.id = seqid # Locus tag of the gene or locus id of the chromosome/plasmid
        self.seq = seq # Sequence of the gene, not used to print
        self.info = description
        self.strain = strain # Strain where the gene is present
        self.gene_type = gene_type # Type of gene/genomic element (glycosyl hydrolase, chromosome, plasmid)
    def get_gene_type(self, type_dict):
        if '.' not in self.id:
            self.gene_type = type_dict[self.id]
        elif 'chromosome' in self.info:
            self.gene_type = 'chromosome'
        elif 'plasmid' in self.info:
            self.gene_type = 'plasmid'
        else:
            self.gene_type = 'genomic'
    def get_GC123(self):
        GC = GC123(self.seq)
        self.GC_total = GC[0] # Mean GC content (%) of a sequence
        self.GC1 = GC[1] # GC content (%) at the first codon position
        self.GC2 = GC[2] # GC content (%) at the second codon position
        self.GC3 = GC[3] # GC content (%) at the third codon position
    def __str__(self): # What print(class) returns
        return f'{self.id}\t{self.strain}\t{self.gene_type}\t{self.GC_total:.2f}\t{self.GC1:.2f}\t{self.GC2:.2f}\t{self.GC3:.2f}\n'


# =============================================================================
# 3. Calculate GC contents
# =============================================================================
# GC_entries = []

no_entries = 0
with open(outfile, 'w') as f:
    f.write('id\tstrain\tseqtype\tGC\tGC1\tGC2\tGC3\n')
for file in GH_files:
    with open(file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            GC_obj = GC_entry(record.id, record.seq, record.description)
            no_entries += 1
            if 'genomic' not in file:
                strain = record.id.split('_')[0]
                if strain == 'APS55':
                    strain = 'MP2'
                elif strain == 'LDX55':
                    strain = 'IBH001'
                elif strain == 'K2W83':
                    strain = 'DSM12361'
            
            else:
                strain = os.path.basename(file).split('_')[0]
            if strain.startswith('H') and '-' not in strain:
                strain = strain[:4] + '-' + strain[4:]
            GC_obj.strain = strain
            GC_obj.get_GC123()
            GC_obj.get_gene_type(type_dict)
            if 'complete' in file:
                GC_obj.gene_type += '_complete'
            with open(outfile, 'a') as fout:
                fout.write(str(GC_obj))
            # GC_entries.append(GC_obj)
        
#genome_subset = [genome for genome in genomes if any(strain in genome for strain in strains)]

# =============================================================================
# 3. Create a text file with a list of genes, type of GH and the GC content at 
# each codon position, including the GC content of the chromosome where the 
# strains are found
# =============================================================================

# =============================================================================
# 3. Function to calculate GC content of the entire genome
# =============================================================================
