#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 17:51:55 2025

Code to retrieve protein sequences from other A. kunkeei species from NCBI,
add the species name to the accession, and save the sequences to a single
fasta file.
Environment: ncbi_download.yml

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import subprocess
import os
from Bio import SeqIO

# =============================================================================
# 1. Set path to output
# =============================================================================

workdir = os.path.expanduser('~') + '/GH_project/add_species/files'
outfile = f'{workdir}/blastp_hits.faa'
finalfile = f'{workdir}/blastp_formatted.faa'

# =============================================================================
# 2. Define the sequences to be downloaded
# =============================================================================
acc_list = ['WP_353322098.1', 'WP_353317313.1', 'WP_249510840.1', 
            'WP_105964111.1', 'WP_105957247.1', 'WP_290047116.1',
            'WP_290047118.1', 'WP_290046731.1', 'WP_290047117.1',
            'WP_105987841.1', 'WP_217304049.1', 'WP_353485634.1', 
            'WP_217304214.1', 'WP_285289173.1', 'WP_220729248.1', 
            #All Apilactobacillus sequences 
            'WP_248720560.1', 'WP_089940200.1', 'WP_248720543.1',
            'WP_252796734.1', 'WP_203617611.1', 'WP_338347939.1',
            'WP_203617613.1', 'WP_412129441.1', 'WP_412142894.1',
            'WP_172692410.1', 'WP_261912693.1', 'WP_036092298.1',
            'WP_111303383.1', 'WP_273750349.1', 'WP_367114384.1',
            'WP_080732772.1', 'WP_010690266.1', 'WP_094511690.1',
            'WP_273774597.1', 'WP_415619990.1', 'WP_415621139.1',
            'WP_056961808.1', 'WP_260116902.1', 'WP_260116903.1',
            'WP_260117220.1', 'WP_186432386.1', 'WP_186432140.1',
            'WP_324249989.1', 'WP_135411726.1', 'WP_071707849.1',
            'WP_231458450.1', 'WP_205143387.1', 'WP_205144328.1',
            'WP_286135565.1', #All Lactobacillaceae sequences
            'WP_080984265.1', 'WP_032672932.1', 'WP_323728076.1', #Sequences from other searches
            'WP_053002791.1', 'WP_013990284.1', 'WP_002891310.1',
            'WP_002884027.1', #Non-Lactobacillaceae sequences
            'AAU08001.1', 'AAU08004.1', 'AAU08015.1',
            'AAU08011.1', 'AAD10952.1', 'CAB65910.2',
            'BAA14241.1', 'AKE50934.1', 'CDX66896.1',
            'CDX65123.1', 'WP_051592287.1', #Sequences from Meng et al. (2018)
            'KDB00248.1' #BrsD
            ]


desc_dict = {'AAU08001.1': 'DSR', 'AAU08004.1': 'MSR', 'AAU08015.1': 'RSR',
             'AAU08011.1': 'DSR',  'AAD10952.1': 'DSR', 'CAB65910.2': 'ASR',
             'BAA14241.1': 'MSR', 'AKE50934.1': 'DSR', 'CDX66896.1': 'a1,2_BRS',
             'CDX65123.1': 'a1,3_BRS', 'WP_051592287.1': 'a1,2_BRS',
             'KDB00248.1': 'a1,2_BRS'}

# acc_list = ['WP_353322098.1', 'WP_353317313.1', 'WP_249510840.1', 
#             'WP_105964111.1', 'WP_105957247.1', 'WP_290047116.1',
#             'WP_290047118.1', 'WP_290046731.1', 'WP_290047117.1',
#             'WP_105987841.1', 'WP_217304049.1', 'WP_353485634.1', 
#             'WP_217304214.1', 'WP_285289173.1', 'WP_220729248.1', 
#             'WP_248720560.1', 'WP_089940200.1', 'WP_248720543.1', 
#             'WP_252795570.1', 'WP_203617611.1', 'WP_338347939.1', 
#             'WP_203617613.1', 'WP_412129441.1', 'WP_412142894.1', 
#             'WP_172692410.1', 'WP_261912693.1', 'WP_036092298.1', 
#             'WP_111303383.1', 'WP_273750349.1', 'WP_367114384.1', 
#             'WP_080732772.1', 'WP_010690266.1', 'WP_094511690.1',
#             'WP_273774597.1', 'WP_415619990.1', 'WP_415621139.1', 
#             'WP_056961808.1', 'WP_260116902.1', 'WP_260116903.1', 
#             'WP_260117220.1', 'WP_186432386.1', 'WP_186432140.1',
#             'WP_231458450.1', 'WP_205143387.1', 'WP_205144328.1',
#             'WP_286135565.1', 'WP_324249989.1', 'WP_135411726.1', 
#             'WP_071707849.1', 'WP_252767194.1']
#            'WP_411292356.1', 'WP_374915322.1', 'WP_196780219.1',
#            'WP_278251525.1']

# =============================================================================
# 3. Create a string with the command to be run
# =============================================================================

list_command = ''
for acc in acc_list:
    list_command += f'{acc} '
list_command = list_command.strip()

full_command = f'ncbi-acc-download --molecule protein --out {outfile} {list_command}'

# =============================================================================
# 4. Run the command
# =============================================================================

subprocess.run(full_command, shell = True)

# =============================================================================
# 5. Read the output file and modify the fasta headers
# =============================================================================
    
records = []
with open(outfile) as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        species = record.description.split('[')[1]
        genus = species.split(' ')[0]
        epithet = species.split(' ')[1][:-1]
        if len(species.split(' ')) > 2 and epithet == 'sp':
            epithet += f'_{species.split(" ")[2][:-1]}'
        elif len(species.split(' ')) > 2:
            epithet += f'i_{species.split(" ")[2][:-1]}'
        elif epithet == 'sp.':
            epithet = epithet[:-1] + '_G418'
        gtype = ''
        if record.id in desc_dict.keys():
            gtype = f'{desc_dict[record.id]}_'
        new_description = f'{record.id}_{gtype}{genus}_{epithet}'
        record.description = new_description
        record.id = new_description
        print(new_description)
        records.append(record)
        
with open(finalfile, 'w') as handle:
    SeqIO.write(records, handle, 'fasta')
        
