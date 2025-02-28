#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:28:42 2023

This is a script to generate tab files to later
use PyGenomeViz to plot the region of the genome where
the NGB gene is in the representative set of strains.

@author: Marina Mota Merlo
"""
from Bio.Blast.Applications import NcbiblastpCommandline as cline_blast
import os

# =============================================================================
# In this section, we run Blast between pairs of strains following the 
# phylogeny.
# =============================================================================

# phylo_order = {1:'G0101', 2:'H4B2-05J', 3:'H4B2-11M', 4:'H1B1-04J', 
#               5:'H4B4-12M', 6:'H4B5-01J', 7:'H3B1-01A', 8:'H3B2-03J', 
#               9:'H3B1-04X', 10:'H3B2-02X', 11:'A1003', 12:'H4B2-04J', 
#               13:'H4B5-03X', 14:'H4B4-02J', 15:'H1B3-02M', 16:'Fhon2', 
#               17:'DSMZ12361', 18:'H4B2-02J', 19:'H1B1-05A', 20:'H3B1-04J', 
#               21:'A1401', 22:'A1202', 23:'H4B5-05J', 24:'H3B2-09X', 
#               25:'H4B5-04J', 26:'A1805', 27:'MP2', 28:'G0403', 29:'H3B2-06M', 
#               30:'H4B4-05J', 31:'H4B1-11J', 32:'H4B4-06M', 33:'A0901', 
#               34:'IBH001', 35:'H4B2-06J', 36:'H3B2-03M', 37:'A1404', 38:'A1001'}

# phylo_order = {1:'G0101', 2:'H4B2-05J', 3:'H4B2-11M', 4:'H4B5-01J', 
#               5:'H4B4-12M', 6:'H1B1-04J', 7:'H3B1-01A', 8:'H3B2-03J', 
#               9:'H4B2-04J', 10:'H3B1-04X', 11:'H3B2-02X', 12:'A1003',
#               13:'H1B3-02M', 14:'H4B4-02J', 15:'H4B5-03X', 16:'Fhon2', 
#               17:'H1B1-05A', 18:'H4B2-02J', 19:'DSMZ12361', 20:'H3B1-04J', 
#               21:'A1202', 22:'A1401', 23:'H4B5-05J', 24:'H3B2-09X', 
#               25:'H4B5-04J', 26:'A1805', 27:'MP2', 28:'G0403', 29:'H4B4-05J', 
#               30:'H3B2-06M', 31:'H4B1-11J', 32:'H4B4-06M', 33:'IBH001', 
#               34:'H4B2-06J', 35:'A0901', 36:'H3B2-03M', 37:'A1404', 38:'A1001'}

phylo_order = {1:'H4B4-12M', 2:'H4B5-01J', 3:'G0101', 4:'H4B2-05J', 
              5:'H4B2-11M', 6:'H1B1-04J', 7:'H4B5-03X', 8:'H4B4-02J', 
              9:'H1B3-02M', 10:'H3B2-02X', 11:'H3B1-04X', 12:'A1003',
              13:'H3B1-01A', 14:'H3B2-03J', 15:'H4B2-04J', 16:'Fhon2', 
              17:'H4B2-02J', 18:'H1B1-05A', 19:'DSMZ12361', 20:'H3B1-04J', 
              21:'H3B2-06M', 22:'H4B4-05J', 23:'H4B1-11J', 24:'H4B4-06M', 
              25:'H3B2-03M', 26:'A0901', 27:'H4B2-06J', 28:'IBH001', 
              29:'H3B2-09X', 30:'H4B5-05J', 31:'A1401', 32:'A1202', 
              33:'H4B5-04J', 34:'MP2', 35:'A1805', 36:'G0403', 
              37:'A1404', 38:'A1001'}


outpath = os.path.expanduser('~') + '/GH_project/blast_tabs'
folder_path = os.path.expanduser('~') + '/Akunkeei_files/fna'

def get_tabs(phylo_dict, folder_path, outpath):
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for i in range(1, len(phylo_dict)):
        strain1_file = f'{folder_path}/{phylo_dict[i]}_genomic.fna'
        strain2_file = f'{folder_path}/{phylo_dict[i+1]}_genomic.fna'
        outfile = f'{outpath}/NGB_{i}.tab'
        cline_input = cline_blast(cmd = 'blastn', query = strain2_file, subject = strain1_file, remote = False, out = outfile, outfmt = 7)
        os.system(str(cline_input));

get_tabs(phylo_order, folder_path, outpath)
