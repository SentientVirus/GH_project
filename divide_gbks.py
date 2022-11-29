#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 10:01:24 2022

@author: marina
"""

from Bio import SeqIO
import os

directory = 'analyses/data/genbanks'

for file in os.listdir(directory):
    for record in SeqIO.parse(file, 'genbank'):
        print(record.id)
