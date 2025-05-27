#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 11:45:45 2023

Script that was used for the Interproscan annotations of strains DSMZ 12361 and
IBH001. It changes the column with protein IDs in the Interproscan output to 
a column with locus tags, so that these files can be navigated using the same
scripts as for the A. kunkeei annotations made by Dyrhage et al. (2022).

@author: Marina Mota Merlo
"""

# =============================================================================
# Script to change protein IDs in Interproscan annotations to locus tags
# =============================================================================

import pandas as pd
from Bio import SeqIO

path = "../Akunkeei_files/gbff"
new_strains = ["DSMZ12361", "IBH001"]
prot2loctag = {}
for s in new_strains:
    with open(f"{path}/{s}_genomic.gbff") as genb:
        info = SeqIO.parse(genb, "genbank")
        for record in info:
            for feature in record.features:
                tag = feature.qualifiers.get('locus_tag')
                prot_id = feature.qualifiers.get('protein_id')
                if tag and prot_id:
                    tag = tag[0]
                    prot_id = prot_id[0]
                    prot2loctag[prot_id] = tag
    with open(f"{s}_id.tsv") as df_file:
        df = pd.read_csv(df_file, sep = "\t", header = None)
        protein_ids = df[0]
        df2 = df.replace({0: prot2loctag})
        df2.to_csv(f"{s}.tsv", sep = "\t", header = False, index = False)
