#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:44:38 2023

@author: marina
"""

import phylonetwork as pn

infile = 'data/fasta/GH70/trees/BRS_repset.mafft.faa.treefile'

str_tree = ''
with open(infile) as tree:
    for line in tree:
        str_tree += line
        
net = pn.PhylogeneticNetwork(eNewick = str_tree)
net.nodes(data = False)
net.draw()