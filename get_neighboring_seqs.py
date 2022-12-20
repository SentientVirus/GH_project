# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 00:34:10 2021

@author: usuario
"""
import os
import re
import csv
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import as_fasta

class gbk_entry:
    def __init__(self, name, subset, locus_tag, translation, gene_seq, description = ''):
        self.name = name
        self.subset = subset
        self.locus_tag = locus_tag
        self.translation = translation
        self.gene_seq = gene_seq
        self.length = len(self.translation)
        self.description = description
    def __str__(self):
        return f'<{self.name} at {self.locus_tag} with length {self.length}, {self.subset} subset>'

same_sequences = {} #Keep track of sequences that are the same
added = [] #Keep track of what has already been added
count = 0
initial_path = '/home/marina/GH_project'
outdir = '/home/marina/GH_project/data/fasta/other_genes'
gbk_dir = '/home/marina/Akunkeei_files/gbff'

gene_dict = {}

if not os.path.exists(outdir):
    os.makedirs(outdir)

# min_length = 0 #667 #Set minimum length for the proteins

# gene_names = {} #Genes to be retrieved and their positions

genes = ['ohrR', 'yifK', 'yhdG', 'GH39', 'nox', 'hpcG', 
         'oppA', 'sbnD', 'epsE', 'epsL', 'ywqE', 'ywqD', 'ywqC', 'tagU'] #Run without GH39, plnV (before and after nox) and ytgP (before sbnD)

gtypes = ['GH70', 'GH32']
for gtype in gtypes:
    if gtype == 'GH70':
        GH70_subset = ['G0403', 'H1B3-02M', 'H3B2-03J', 'H3B1-04J', 'H3B1-04X', 
                      'H4B2-02J', 'H4B2-04J', 'H4B5-04J', 'H4B5-05J']
    else:
        GH70_subset = ['A1001', 'A1805', 'H1B1-04J', 'H3B2-03M',
                      'H4B4-02J', 'H4B4-12M', 'H4B5-03X']
        
    check = False
    for gene_n in genes:
            
        
        #Change everything from here to read the GenBank directly
        
        # strain_list = []
        # gene_list = []
        
        # with open(file_out, 'w') as gn: #Reset file
        #     gn.write('')
        
        def sort_key(text:str):
            return(text[0].upper(), text[1:])
        
        os.chdir(gbk_dir)
        file_list = sorted(os.listdir(os.getcwd()), key=sort_key)
        for filename in file_list:
            sname = filename.split('_')[0]
            
                    
            if filename.endswith('.gbff') and sname in GH70_subset:
                # if ('H3B1-04' in filename or 'H3B2-03J' in filename) and gene_n == 'oppA_2':
                #     gene_n = 'oppA_3'
                # elif 'H1B3-02M' in filename and gene_n == 'oppA_2':
                #     gene_n = 'oppA'
                with open(filename) as handle:
                    for record in SeqIO.parse(handle, 'genbank'):
                        for feature in record.features:
                            if 'gene' in feature.qualifiers.keys() and 'translation' in feature.qualifiers.keys():
                                gene_str = feature.qualifiers['gene'][0] 
                                locus_tag = feature.qualifiers['locus_tag'][0][3:]
                                #Redo this part and use dictionaries to take the highest or lowest locus tag depending on the gene
                                if (gene_str == gene_n and (int(locus_tag.split('_')[1]) > 12000)) or ('mhpD' in gene_str and gene_n == 'hpcG'):
                                    if not (gtype == 'GH70' and gene_n == 'yhdG' and (int(locus_tag.split('_')[1]) < 12200)):
                                        if (gene_n, gtype) not in gene_dict.keys():
                                            gene_dict[(gene_n, gtype)] = []
                                            file_out = f'{outdir}/a_kunkeei_{gene_n}_{gtype}.faa'
                                            file_fna = file_out.replace('.faa', '.fna')
                                            if gene_n != 'yhdG':
                                                with open(file_out, 'w') as gn:
                                                    gn.write('')
                                        elif (sname.replace('-', '') in [gbk.locus_tag.split('_')[0] for gbk in gene_dict[(gene_n, gtype)]]):
                                            gene_n = f'{gene_n}_2'
                                            if (gene_n, gtype) not in gene_dict.keys() and gene_n != 'sbnD_2' and gene_n != 'GH39_2':
                                                gene_dict[(gene_n, gtype)] = []
                                                file_out = f'{outdir}/a_kunkeei_{gene_n}_{gtype}.faa' #Add .fna files to code
                                                file_fna = file_out.replace('.faa', '.fna')
                                                with open(file_out, 'w') as gn, open(file_fna, 'w') as gn2:
                                                    gn.write('')
                                                    gn2.write('')
                                        
                                        gbkobj = gbk_entry(gene_n, gtype, locus_tag, 
                                                           feature.qualifiers['translation'][0], 
                                                           feature.location.extract(record).seq,
                                                           feature.qualifiers['product'][0])
                                        
                                        if gbkobj.name != 'sbnD_2' and gbkobj.name != 'GH39_2':
                                            gene_dict[(gene_n, gtype)].append(gbkobj)
                                        gene_n = gene_n.split('_')[0]
                                        
                                        #Filter out yhdG
                                        file_out = f'{outdir}/a_kunkeei_{gbkobj.name}_{gtype}.faa'
                                        file_fna = file_out.replace('.faa', '.fna')
                                        if (gbkobj.name != 'sbnD_2') and (gbkobj.name != 'yhdG') and ('GH39' not in gbkobj.name or (gbkobj.name == 'GH39' and gbkobj.length > 500)):
                                            with open(file_out, 'a') as outfile:
                                                sequence = Seq(gbkobj.translation)
                                                entry = SeqRecord(sequence, gbkobj.locus_tag, 
                                                                  gbkobj.name, description = gbkobj.description)
                                                outfile.write(as_fasta(entry))
                                            with open(file_fna, 'a') as outfile2:
                                                gene_entry = SeqRecord(gbkobj.gene_seq, gbkobj.locus_tag, 
                                                                  gbkobj.name, description = gbkobj.description)
                                                outfile2.write(as_fasta(gene_entry))
                                            print(gbkobj)
                                        
                                        
                            # with open(file_out, 'a') as gn:
                            #     if record.id == gene_n:
                            #         print('It works!')
                            #         file_out.write(record)
                            #         print(f'{record.id} with length {len(record.seq)} added.\n')
        #         if ('H3B1-04' in filename or 'H3B2-03J' in filename) and gene_n == 'oppA_2':
        #             gene_n = 'oppA_3'
        #         elif 'H1B3-02M' in filename and gene_n == 'oppA_2':
        #             gene_n = 'oppA'
        #         f = (row for row in open(filename))
        #         locus = next(f)[12:22] #Take file name (i.e. HXBX-XXX)
        #         locus = locus.replace(' ', '')
        
        #         for line in f: #Loop through lines in file
        #             #print(line[22:27])
        #             if line[22:27] == 'locus': #Get gene names
        #                 gene = line[:-1].replace(' ', '') #Remove spaces
        #                 gene = gene.replace('/locus_tag=', '') #Remove other text
        #                 gene = re.sub(r'[("")]', '', gene) #Remove ""
        #                 line2 = next(f)
        #                 line2 = next(f)
        #                 # if 'gene' in line2:
        #                 #     gene_name = line2.replace('/gene=', '')
        #                 #     print(gene_name)
        #                 if (gene_n in line2) or ('A1001' in filename and 'yhdG' in gene_n and 'yhdG_2' in line2) or ('mhpD' in line2 and gene_n == 'hpcG'): #If it is the gene(s) of interest
        #                     while 'translation' not in line2: #Jump to the start of the protein sequence
        #                         line2 = next(f)
        #                     prot_seq = '' #Define string
        #                     while line2[21:24] !='/EC' and line2[5:9] != 'gene': #Get full protein sequence
        #                         prot_seq += line2
        #                         line2 = next(f)
        #                     prot_seq = prot_seq.replace(' ', '')
        #                     prot_seq = prot_seq.replace('/translation=', '')
        #                     prot_seq = re.sub(r'[("\n")]', '', prot_seq)
        #                     if 'oppA' in gene_n:
        #                         gene_n = 'oppA_2'
        
        
        #                     with open(file_out, 'r') as gn: #Add sequence info to file
        #                         prot = prot_seq
        #                         n = 1
        #                         i = 0
        #                         while i < n:
        #                             if prot not in added:
        #                                 added.append(prot)
        #                                 same_sequences[prot] = [gene]
        #                             else:
        #                                 same_sequences[prot].append(gene)
        #                             with open(file_out, 'a+') as gn2:
        #                                 record = SeqRecord(Seq(prot), gene, gene, '')
        #                                 fasta_record = as_fasta(record)
        #                                 gn2.write(fasta_record)
        #                                 print(f'{gene} added, {len(prot)}')
        #                                 count += 1
        #                             gene_list.append(gene)
        #                             i += 1
        
        # x = f'{gene_n}_{gtype}_protein'
        # print(f'{outdir}/same_{x}_sequences.txt')
        # with open(f'{outdir}/same_{x}_sequences.txt', 'w') as gn2:
        #     x = x.replace('_', ' ')
        #     gn2.write(f'There are a total of {count} ({len(same_sequences.keys())}) different sequences at the {x} level.\n')
        #     gn2.write('The following genomes contain the same sequence:\n\n')
        #     for seq in same_sequences.keys():
        #         [gn2.write(f'{sequence} ') for sequence in same_sequences[seq]]
        #         gn2.write('\n')
        #         #[gn2.write(f'{g}\t{same_sequences[seq][0]}\n') for g in same_sequences[seq]]
