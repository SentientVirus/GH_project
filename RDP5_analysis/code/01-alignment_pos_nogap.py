#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 13:31:36 2022

@author: marina
"""

import logging, traceback, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import subprocess

# =============================================================================
# 0. Logging
# =============================================================================
log = os.path.expanduser('~') + '/GH_project/RDP5_analysis/logs/get_alignments.log'

# if os.path.exists(log):
#     os.remove(log)

# if not os.path.exists(os.path.dirname(log)):
#     os.makedirs(os.path.dirname(log))

# logger = logging.getLogger()

# logging.basicConfig(filename = log, level = logging.INFO,
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

# sys.stdout = open(log, 'a')

# =============================================================================
# 0. Input definition
# =============================================================================
# Paths
gbk_dir = os.path.expanduser('~') + '/Akunkeei_files/gbff'
fna_dir = os.path.expanduser('~') + '/Akunkeei_files/fna'
infile_suffix = 'genomic.gbff'
infile_suffix2 = 'genomic.fna'
pos_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/tab/all_subsets_nogap'
out_dir = os.path.expanduser('~') + '/GH_project/RDP5_analysis/files/fna/all_subsets_nogap'
threads = 48 #Threads to run mafft

# Strains and comparisons of interest
# strain_groups = {'H3B2-03M': 0, 'H4B4-02J': 0, 'H4B5-03X': 0, 'H4B4-12M': 0, 
#                  'H4B4-06M': 0, 'H1B1-04J': 0, 'A0901': 0, 'H1B3-02M': 0,
#                  'H3B2-09X': 1, 'H4B5-05J': 1, 'H3B1-04X': 1, 'H4B5-04J': 1,
#                  'MP2': 2, 'H3B2-02X': 2, 'H3B2-03J': 2, 'G0403': 2}

group_names = {0: 'root_GS1_S2-3_subset', 1: 'GS1-2_BRS', 2: 'only_GS1+GS2'}

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

yhdG_loctag = {'A0901': 'AKUA0901_13230', 'A1003': 'AKUA1003_12490', 
               'A1202': 'AKUA1202_13470', 'A1401': 'AKUA1401_12710', #'A1404': 'AKUA1404_13400',
               'A1805': 'AKUA1805_12730', 'DSMZ12361': 'K2W83_RS06150',
               'Fhon2': 'AKUFHON2_13490', 'G0101': 'AKUG0101_12730',
               'G0403': 'AKUG0403_13060', 'H1B1-04J': 'AKUH1B104J_12950',
               'H1B1-05A': 'AKUH1B105A_12240', 'H1B3-02M': 'AKUH1B302M_12800', 
               'H3B1-01A': 'AKUH3B101A_13200', 'H3B1-04J': 'AKUH3B104J_12940',
               'H3B1-04X': 'AKUH3B104X_13160', 'H3B2-02X': 'AKUH3B202X_12810', 
               'H3B2-03J': 'AKUH3B203J_13330', 'H3B2-03M': 'AKUH3B203M_12410', 
               'H3B2-06M': 'AKUH3B206M_12790', 'H3B2-09X': 'AKUH3B209X_13290',
               'H4B1-11J': 'AKUH4B111J_13520', 'H4B2-02J': 'AKUH4B202J_12840',
               'H4B2-04J': 'AKUH4B204J_13290', 'H4B2-05J': 'AKUH4B205J_12930',
               'H4B2-06J': 'AKUH4B206J_13360', 'H4B2-11M': 'AKUH4B211M_12940',
               'H4B4-02J': 'AKUH4B402J_12550', 'H4B4-05J': 'AKUH4B405J_13310',
               'H4B4-06M': 'AKUH4B406M_13410', 'H4B4-12M': 'AKUH4B412M_13180', 
               'H4B5-01J': 'AKUH4B501J_12830', 'H4B5-03X': 'AKUH4B503X_12620', 
               'H4B5-04J': 'AKUH4B504J_13420', 'H4B5-05J': 'AKUH4B505J_12830',
               'IBH001': 'LDX55_06305', 'MP2': 'APS55_RS03880'}

CDS8_loctag = {'A0901': 'AKUA0901_13240', 'A1003': 'AKUA1003_12500', 
               'A1202': 'AKUA1202_13480', 'A1401': 'AKUA1401_12720', #'A1404': 'AKUA1404_13400',
               'A1805': 'AKUA1805_12740', 'DSMZ12361': 'K2W83_RS06155',
               'Fhon2': 'AKUFHON2_13500', 'G0101': 'AKUG0101_12740',
               'G0403': 'AKUG0403_13070', 'H1B1-04J': 'AKUH1B104J_12960',
               'H1B1-05A': 'AKUH1B105A_12270', 'H1B3-02M': 'AKUH1B302M_12830', 
               'H3B1-01A': 'AKUH3B101A_13210', 'H3B1-04J': 'AKUH3B104J_12950',
               'H3B1-04X': 'AKUH3B104X_13170', 'H3B2-02X': 'AKUH3B202X_12820', 
               'H3B2-03J': 'AKUH3B203J_13340', 'H3B2-03M': 'AKUH3B203M_12440', 
               'H3B2-06M': 'AKUH3B206M_12800', 'H3B2-09X': 'AKUH3B209X_13300',
               'H4B1-11J': 'AKUH4B111J_13530', 'H4B2-02J': 'AKUH4B202J_12850',
               'H4B2-04J': 'AKUH4B204J_13300', 'H4B2-05J': 'AKUH4B205J_12940',
               'H4B2-06J': 'AKUH4B206J_13370', 'H4B2-11M': 'AKUH4B211M_12950',
               'H4B4-02J': 'AKUH4B402J_12560', 'H4B4-05J': 'AKUH4B405J_13320',
               'H4B4-06M': 'AKUH4B406M_13420', 'H4B4-12M': 'AKUH4B412M_13190', 
               'H4B5-01J': 'AKUH4B501J_12840', 'H4B5-03X': 'AKUH4B503X_12630', 
               'H4B5-04J': 'AKUH4B504J_13430', 'H4B5-05J': 'AKUH4B505J_12840',
               'IBH001': 'LDX55_06310', 'MP2': 'APS55_RS03865'}

nox_loctag = {'A0901': 'AKUA0901_13250', 'A1003': 'AKUA1003_12520', 
               'A1202': 'AKUA1202_13500', 'A1401': 'AKUA1401_12730', #'A1404': 'AKUA1404_13400',
               'A1805': 'AKUA1805_12750', 'DSMZ12361': 'K2W83_RS06165',
               'Fhon2': 'AKUFHON2_13510', 'G0101': 'AKUG0101_12760',
               'G0403': 'AKUG0403_13080', 'H1B1-04J': 'AKUH1B104J_12970',
               'H1B1-05A': 'AKUH1B105A_12280', 'H1B3-02M': 'AKUH1B302M_12840', 
               'H3B1-01A': 'AKUH3B101A_13220', 'H3B1-04J': 'AKUH3B104J_12970',
               'H3B1-04X': 'AKUH3B104X_13180', 'H3B2-02X': 'AKUH3B202X_12830', 
               'H3B2-03J': 'AKUH3B203J_13350', 'H3B2-03M': 'AKUH3B203M_12450', 
               'H3B2-06M': 'AKUH3B206M_12810', 'H3B2-09X': 'AKUH3B209X_13310',
               'H4B1-11J': 'AKUH4B111J_13540', 'H4B2-02J': 'AKUH4B202J_12860',
               'H4B2-04J': 'AKUH4B204J_13310', 'H4B2-05J': 'AKUH4B205J_12950',
               'H4B2-06J': 'AKUH4B206J_13380', 'H4B2-11M': 'AKUH4B211M_12960',
               'H4B4-02J': 'AKUH4B402J_12570', 'H4B4-05J': 'AKUH4B405J_13330',
               'H4B4-06M': 'AKUH4B406M_13430', 'H4B4-12M': 'AKUH4B412M_13200', 
               'H4B5-01J': 'AKUH4B501J_12850', 'H4B5-03X': 'AKUH4B503X_12640', 
               'H4B5-04J': 'AKUH4B504J_13440', 'H4B5-05J': 'AKUH4B505J_12850',
               'IBH001': 'LDX55_06315', 'MP2': 'APS55_RS03860'}

CDS7_loctag = {'A0901': 'AKUA0901_13260', 'A1003': 'AKUA1003_12530', 
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

CDS5_loctag =  {'A0901': 'AKUA0901_13400', 'A1003': 'AKUA1003_12580',
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

ywqC_loctag = {'A0901': 'AKUA0901_13500', 'A1003': 'AKUA1003_12680',
               'A1202': 'AKUA1202_13660', 'A1401': 'AKUA1401_12910', 
               'A1805': 'AKUA1805_12960', 'DSMZ12361': 'K2W83_RS06255',
               'Fhon2': 'AKUFHON2_13730', 'G0101': 'AKUG0101_12970',
               'G0403': 'AKUG0403_13260', 'H1B1-04J': 'AKUH1B104J_13200',
               'H1B1-05A': 'AKUH1B105A_12440', 'H1B3-02M': 'AKUH1B302M_13040',
               'H3B1-01A': 'AKUH3B101A_13400', 'H3B1-04J': 'AKUH3B104J_13160',
               'H3B1-04X': 'AKUH3B104X_13390', 'H3B2-02X': 'AKUH3B202X_13000',
               'H3B2-03J': 'AKUH3B203J_13530', 'H3B2-03M': 'AKUH3B203M_12630',
               'H3B2-06M': 'AKUH3B206M_12980', 'H3B2-09X': 'AKUH3B209X_13530', 
               'H4B1-11J': 'AKUH4B111J_13710', 'H4B2-02J': 'AKUH4B202J_13030', 
               'H4B2-04J': 'AKUH4B204J_13510', 'H4B2-05J': 'AKUH4B205J_13160',
               'H4B2-06J': 'AKUH4B206J_13630', 'H4B2-11M': 'AKUH4B211M_13170',
               'H4B4-02J': 'AKUH4B402J_12770', 'H4B4-05J': 'AKUH4B405J_13500',
               'H4B4-06M': 'AKUH4B406M_13600', 'H4B4-12M': 'AKUH4B412M_13440',
               'H4B5-01J': 'AKUH4B501J_13060', 'H4B5-03X': 'AKUH4B503X_12870',
               'H4B5-04J': 'AKUH4B504J_13620', 'H4B5-05J': 'AKUH4B505J_13080',
               'IBH001': 'LDX55_06405', 'MP2': 'APS55_RS03775'}

tagU_loctag = {'A0901': 'AKUA0901_13510', 'A1003': 'AKUA1003_12690',
               'A1202': 'AKUA1202_13670', 'A1401': 'AKUA1401_12920', 
               'A1805': 'AKUA1805_12970', 'DSMZ12361': 'K2W83_RS06260',
               'Fhon2': 'AKUFHON2_13740', 'G0101': 'AKUG0101_12980',
               'G0403': 'AKUG0403_13270', 'H1B1-04J': 'AKUH1B104J_13210',
               'H1B1-05A': 'AKUH1B105A_12450', 'H1B3-02M': 'AKUH1B302M_13050',
               'H3B1-01A': 'AKUH3B101A_13410', 'H3B1-04J': 'AKUH3B104J_13170',
               'H3B1-04X': 'AKUH3B104X_13400', 'H3B2-02X': 'AKUH3B202X_13010',
               'H3B2-03J': 'AKUH3B203J_13540', 'H3B2-03M': 'AKUH3B203M_12660',
               'H3B2-06M': 'AKUH3B206M_12990', 'H3B2-09X': 'AKUH3B209X_13540', 
               'H4B1-11J': 'AKUH4B111J_13720', 'H4B2-02J': 'AKUH4B202J_13040', 
               'H4B2-04J': 'AKUH4B204J_13520', 'H4B2-05J': 'AKUH4B205J_13170',
               'H4B2-06J': 'AKUH4B206J_13640', 'H4B2-11M': 'AKUH4B211M_13180',
               'H4B4-02J': 'AKUH4B402J_12780', 'H4B4-05J': 'AKUH4B405J_13510',
               'H4B4-06M': 'AKUH4B406M_13610', 'H4B4-12M': 'AKUH4B412M_13450',
               'H4B5-01J': 'AKUH4B501J_13070', 'H4B5-03X': 'AKUH4B503X_12880',
               'H4B5-04J': 'AKUH4B504J_13630', 'H4B5-05J': 'AKUH4B505J_13090',
               'IBH001': 'LDX55_06410', 'MP2': 'APS55_RS03770'}

tagU_plus2_loctag = {'A0901': 'AKUA0901_13530', 'A1003': 'AKUA1003_12720',
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
if not os.path.exists(pos_dir):
    os.makedirs(pos_dir)
    
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
with open(f'{pos_dir}/all_subsets_positions.tab', 'w+') as gpos:
    gpos.write('locus_tag\tstrain\tgene_name\tstart\tend\tstrand\n')
    
for ind in range(1, 6):
    rawdata = f'{out_dir}/all_subsets_seqs{ind}.fasta'
    with open(rawdata, 'w+') as seq_no:
        seq_no.write('')

full_record = {}
dict_pos = {}
all_genes = {}
genes_in_segment = {}
loop =  list(bcrA_loctag.keys())
# loop.remove('H1B3-02M')
for strain in loop:
    segment1 = [0, 0]
    segment2 = [0, 0]
    segment3 = [0, 0]
    segment4 = [0, 0]
    segment5 = [0, 0]
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
                            segment1[0] = gene_obj.start
                        else:
                            segment5[1] = gene_obj.end
                    elif gene_obj.locus_tag == yhdG_loctag[strain]:
                        if strain != 'MP2':
                            segment1[1] = gene_obj.end
                        else:
                            segment5[0] = gene_obj.start
                    elif gene_obj.locus_tag == CDS8_loctag[strain]:
                        if strain != 'MP2':
                            segment2[0] = gene_obj.start
                        else:
                            segment4[1] = gene_obj.end
                    elif gene_obj.locus_tag == nox_loctag[strain]:
                        if strain != 'MP2':
                            segment2[1] = gene_obj.end
                        else:
                            segment4[0] = gene_obj.start
                    elif gene_obj.locus_tag == CDS7_loctag[strain]:
                        segment3[0] = gene_obj.start
                        segment3[1] = gene_obj.end
                    elif gene_obj.locus_tag == CDS5_loctag[strain]:
                        if strain != 'MP2':
                            segment4[0] = gene_obj.start
                        else:
                            segment2[1] = gene_obj.end
                    elif gene_obj.locus_tag == ywqC_loctag[strain]:
                        if strain != 'MP2':
                            segment4[1] = gene_obj.end
                        else:
                            segment2[0] = gene_obj.start
                    elif gene_obj.locus_tag == tagU_loctag[strain]:
                        if strain != 'MP2':
                            segment5[0] = gene_obj.start
                        else:
                            segment1[1] = gene_obj.end
                    elif gene_obj.locus_tag == tagU_plus2_loctag[strain]:
                        if strain != 'MP2':
                            segment5[1] = gene_obj.end
                        else:
                            segment1[0] = gene_obj.start
                      
    dict_pos[strain] = [segment5, segment4, segment3, segment2, segment1]
     
    cummulength = 0
    # for segment in dict_pos[strain]:
    #     cummulength += segment[1]-segment[0]
    for i in range(0, len(dict_pos[strain])):
        j = len(dict_pos[strain])-1-i
        if i > 0:
            cummulength += dict_pos[strain][i-1][1]-dict_pos[strain][i-1][0]
        for gene in all_genes[strain]:
            check = False
            if strain == 'MP2' and dict_pos[strain][j][0] <= gene.start < dict_pos[strain][j][1]:
                gene.start = gene.start - dict_pos[strain][j][0] + 1
                gene.end = gene.end - dict_pos[strain][j][0] + 1
                check = True
                genes_in_segment[strain].append(gene)
            elif strain != 'MP2' and dict_pos[strain][i][0] <= gene.start < dict_pos[strain][i][1]:
                    
                
                print(i, cummulength)
                gene_start = gene.start
                gene.start = dict_pos[strain][i][1] - gene.end + cummulength + 2
                gene.end = dict_pos[strain][i][1] - gene_start + cummulength + 2
                check = True
                genes_in_segment[strain].append(gene)
        
            if check:
                with open(f'{pos_dir}/all_subsets_positions.tab', 'a+') as pos_file:
                    pos_file.write(f'{gene.locus_tag}\t{strain}\t{gene.name}\t{gene.start}\t{gene.end}\t{gene.strand}\n')
                        
                        
        with open(f'{fna_dir}/{strain}_{infile_suffix2}') as fna:
            fna_read = SeqIO.parse(fna, 'fasta')
            count = 0
            for record in fna_read:
                if count > 0:
                    break
                if strain != 'MP2':
                    segment = record.seq[dict_pos[strain][j][0]:dict_pos[strain][j][1]+1] #+ '!'
                else: # Almost fixed, there is now a single 1 kb gap in MP2, possibly transposons?
                    genome_len = len(record.seq)
                    start_s = genome_len - dict_pos[strain][i][1]
                    end_s = genome_len - dict_pos[strain][i][0] + 1
                    segment = record.seq[start_s:end_s]
                    
                record.id = strain
                record.description = ''
                record.seq = segment
                full_record[strain] = record.seq

                record.seq = record.seq.reverse_complement()
                record.id = f'{strain}_{i+1}'
                with open(f'{out_dir}/all_subsets_seqs{i+1}.fasta', 'a+') as nogap:
                    print(f'Saving record {record.id} to file {out_dir}/all_subsets_seqs{i+1}.fasta')
                    SeqIO.write(record, nogap, 'fasta')
                    count += 1
       
with open(f'{out_dir}/all_subsets_seqs.mafft.fasta', 'w+') as last_outfile:
    last_outfile.write('')
    
record_dict = {}
for ind in range(1, 6):
    raw_data = f'{out_dir}/all_subsets_seqs{ind}.fasta'
    alignment = raw_data.replace('fasta', 'mafft.fasta')
    subprocess.run(f'mafft --auto --thread {threads} {raw_data} > {alignment} 2>> {log}', shell = True)
    
    
    with open(f'{out_dir}/all_subsets_seqs{ind}.mafft.fasta') as seq:
        records = list(SeqIO.parse(seq, 'fasta'))
        for k in range(0, len(records)):
            strain = records[k].id.split('_')[0]
            if strain not in record_dict:
                record_dict[strain] = ''
            record_dict[strain] += records[k].seq

for strain in record_dict.keys():
    seqrec = SeqRecord(record_dict[strain], id = strain, name = strain, 
                       description = '')
    with open(f'{out_dir}/all_subsets_seqs.mafft.fasta', 'a+') as last_outfile:
        SeqIO.write(seqrec, last_outfile, 'fasta')
        print(f'Writing concatenated alignment of strain {strain} to file {out_dir}/all_subsets_seqs.mafft.fasta')
                    
