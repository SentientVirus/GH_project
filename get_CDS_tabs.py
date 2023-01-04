# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 21:57:29 2021

@author: usuario
"""
import os
import re

homedir = '/home/marina'
outdir = 'plots/tabfiles'

if not os.path.exists(outdir):
    os.makedirs(outdir)

#Create a dictionary to plot strands according to the phylogeny
phylogeny = {}
n = 1
with open('phylogeny.txt') as f:
    lines = f.readlines()
    for line in lines:
        line = line.replace(' ', '')
        line = line.replace('\n', '')
        no = str(n)
        while len(no) < 3:
            no = '0' + no
        n += 1
        phylogeny[line] = no

'''#Create a dictionary to plot strands according to the phylogeny
phylogeny =  {'G0102': '001', 'G0401': '002', 'G0803': '003', 'G0804': '004',
              'G0404': '005', 'G0801': '006', 'G0802': '007', 'G0402': '008',
              'G0602': '009', 'G0408': '009', 'G0'}'''

#Change to the directory where the gbk files are
os.chdir('../Akunkeei_files/gbff') #Replace with a parameter from Snakemake

class CDS:
    def __init__(self, name, start, end, strand, locus):
        self.name = name #3-letter code for the gene
        self.start = start #Start position of the protein gene
        self.end = end #End position of the gene
        self.strand = strand #Strand where the gene is (+ or -)
        self.locus = locus #Name of the file (genome)
    def __str__(self):
        return f'Protein: {self.name} in genome {self.locus}, {self.start}-{self.end}, strand {self.strand}'

omit = ['Fhon13', 'H3B2-03M-01', 'H3B2-03M-04'] #Or use shadow: minimal in Snakemake

for filename in os.listdir(os.getcwd()):
    if filename.split('_')[0] not in omit and filename.endswith('.gbff'):
        f = (row for row in open(filename))
        locus = next(f)[12:23] #Take file name (i.e. HXBX-XXX)
            
        for line in f:
            if 'strain=' in line:
                strain = ''.join(e for e in line if (e.isalnum() or e == '-'))[6:]
                with open(f'../../GH_project/{outdir}/{phylogeny[strain]}_{strain}_gpr.tab', 'w') as gn: #Start writing tab file
                    gn.write('name\tstart\tend\tstrand\tregion\n')
            if line[5:8] == 'CDS':
                check = False
                if 'complement' in line: #If the gene is in the complementary strand
                    line = re.sub(r'[(complement)]', '', line) #Remove the word "complement"
                    strand = '-' #Save strand
                else:
                    strand = '+'
                position = line[21:-1].split('..') #Save position of CDS
                for i in range(8):
                    gene = next(f) #Go to the next line
                    if '/gene' in gene: #If the next line includes gene name
                        gene = gene[:-1].replace(' ', '') #Remove spaces
                        gene = gene.replace('/gene=', '') #Remove other text
                        gene = re.sub(r'[("")]', '', gene) #Remove ""
                        check = True
                        break
                        #print(gene)
                    if '/locus_tag' in gene: #If the next line includes gene name
                        locus_tag = gene[:-1].replace(' ', '') #Remove spaces
                        locus_tag = locus_tag.replace('/locus_tag=', '') #Remove other text
                        locus_tag = re.sub(r'[("")]', '', locus_tag) #Remove ""
                        print(locus_tag)
                if check == False:
                    gene = locus_tag
                p = CDS(gene, position[0], position[1], strand, locus) #Save info in object

                with open(f'{phylogeny[strain]}_{filename[:-4]}_gpr.tab', 'a') as gn: #Add object info to file
                    gn.write(f'{p.name}\t{p.start}\t{p.end}\t{p.strand}\t{p.locus}\n')
