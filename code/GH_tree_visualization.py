# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 16:04:16 2021

@author: usuario
"""
import os
import csv
#Idea for this script:
#   1. Take all locus ids from tree file (NO, filter them in Interproscan)
#   2. Create protein objects based on Interproscan info
#       ·Include protein length and domain subobjects with length, chosen shape, color, etc.
#   3. Write to file with | separators
class domain_obj:
    def __init__(self, locus, start, end, dom_type, shape, color):
        self.locus = locus
        self.start = int(start)
        self.end = int(end)
        self.length = max(self.start, self.end) - min(self.start, self.end)
        self.type = dom_type
        self.shape = shape
        self.color = color
    def __str__(self):
        return f'<{self.type} domain of length {self.length} at locus {self.locus}>'
    
class protein_obj:
    def __init__(self, locus, length, domains = []):
        self.locus = locus
        self.length = int(length)
        self.domains = domains
    def __str__(self):
        return f'<{len(self.domains)}-domain protein of length {self.length} at locus {self.locus}>'

types = {'SP': ['Signal peptide', 'DI', '#00197F'], 'GH70': ['Glycosyl hydrolase family 70', 'EL', '#FF006F'], 'DUF': ['DUF5776', 'TR', '#7D8292'],'GH32': ['Glycosyl hydrolase family 32', 'EL', '#EE3900'], 'CB': ['Choline-binding', 'RE', '#FFAB00'], 'CWB': ['Cell-wall binding', 'RE', '#FF8900']}
indir = '../../interproscan'

GH_types = ['GH70', 'GH32']
for GH_type in GH_types:
    ref_tree = f'../data/fasta/{GH_type}/trees/{GH_type}_functional_repset.mafft.faa.treefile'
    if GH_type == 'GH32':
        ref_tree = ref_tree.replace('functional_', '')
    dom_out = f'../data/tabs/{GH_type}_domain_file.txt'
    
    proteins = {}
    locus_list = []
        
    with open(ref_tree) as ref: #Use the tree with all domains as a reference
        ref_text = ref.read()
        print(ref_text)
    
    for filename in os.listdir(indir):
        with open(f'{indir}/{filename}') as file:
            domain_tab = csv.reader(file, delimiter='\t')
            for line in domain_tab:
                sig_pep = False
                locus = line[0].replace('-', '').upper()
                if locus == 'MP2_13360':
                    locus = 'APS55_RS03845'
                elif locus == 'MP2_13350':
                    locus = 'APS55_RS03850'
                
                if locus in ref_text:
                    print(locus)
                    if ('GH32' in ref_tree) and ((line[3] == 'Pfam' and 'Choline' in line[5]) or ('A1404' in locus and line[3] == 'Phobius' and 'Signal peptide region' in line[5]) or (line[3] == 'Pfam' and 'signal peptide' in line[5]) or (line[3] == 'TIGRFAM' and 'signal peptide' in line[5]) or (line[3] == 'SMART' and 'glyco_32' in line[5])):
                        if line[3] == 'TIGRFAM' and sig_pep == False:
                            dom_type = 'SP'
                            sig_pep = True
                        elif line[3] == 'SMART':
                            dom_type = 'GH32'
                        elif line[3] == 'Pfam':
                            if 'Choline' in line[5]:
                                dom_type = 'CB'
                            elif sig_pep == False: 
                                dom_type = 'SP'
                                sig_pep = True
                        if line[3] == 'Phobius' and sig_pep == False:
                            dom_type = 'SP'
                            sig_pep = False
                        
                        domain = domain_obj(locus, line[6], line[7], dom_type, types[dom_type][1], types[dom_type][2])
                        if locus not in locus_list:
                            protein = protein_obj(locus, line[2], [domain])
                            proteins[locus] = protein
                            locus_list.append(locus)
                        if domain not in proteins[locus].domains and domain.locus == locus:
                            proteins[locus].domains.append(domain)
                            
                    elif ('GH70' in ref_tree) and (line[3] == 'Pfam' or (line[3] == 'TIGRFAM' and 'signal peptide' in line[5])):
                        if line[3] == 'TIGRFAM' or (line[3] == 'Pfam' and 'signal peptide' in line[5]) and sig_pep == False:
                            dom_type = 'SP'
                            sig_pep = True
                        elif 'family 70' in line[5]:
                            dom_type = 'GH70'
                        elif 'Choline' in line[5]:
                            dom_type = 'CB'
                        elif 'cell wall' in line[5]:
                            dom_type = 'CWB'
                        elif 'DUF' in line[5]:
                            dom_type = 'DUF'
                            
                        domain = domain_obj(locus, line[6], line[7], dom_type, types[dom_type][1], types[dom_type][2])
                        if locus not in locus_list:
                            protein = protein_obj(locus, line[2], [domain])
                            proteins[locus] = protein
                            locus_list.append(locus)
                        if domain not in proteins[locus].domains and domain.locus == locus:
                            proteins[locus].domains.append(domain)
      
                            
    outdir = os.path.dirname(dom_out)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    with open(dom_out, 'w') as dom_file:
        dom_file.write('DATASET_DOMAINS\n')
        dom_file.write('SEPARATOR COMMA\n')
        dom_file.write('DATASET_LABEL,gtf\n')
        dom_file.write('COLOR,#ff0000\n')
        dom_file.write('DATA\n')
        for l in proteins.keys():
            node = l
            length = proteins[l].length
            line = f'{node},{length},'
            for domain in proteins[l].domains:
                line += f'{domain.shape}|{domain.start}|{domain.end}|'
                line += f'{domain.color.lower()}|{domain.type},'
            line = line[:-1] + '\n'
            dom_file.write(line)
            