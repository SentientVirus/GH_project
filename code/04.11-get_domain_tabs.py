# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 16:04:16 2021

Script to retrieve the plotting information for the CDS of the genes
in a region of interest so that they can be later plotted next to
the tree in ete3.

Note: Unzip the .zip files with SignalP6 slow predictions before running this
script.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import Python libraries
# =============================================================================
import os
import csv
import pandas as pd

# =============================================================================
# 1. Define classes to store domain and protein information.
# =============================================================================
class domain_obj:
    '''Object with domain information.
        locus: Locus tag of the CDS where the domain is present.
        start: Start position of the domain in the gene.
        end: End position of the domain in the gene.
        dom_type: The type of domain (glycosyl hydrolase, signal peptide, etc.)
        length: Length of the domain.
        shape: Shape with which the domain should be plotted.
        color: Color of the shape that represents the domain.
        '''
    def __init__(self, locus, start, end, dom_type, shape, color): #Initialize domain object
        #Assign properties to object    
        self.locus = locus 
        self.start = int(start)
        self.end = int(end)
        self.length = max(self.start, self.end) - min(self.start, self.end)
        self.type = dom_type
        self.shape = shape
        self.color = color
    def __str__(self): #Object information that is printed
        return f'<{self.type} domain of length {self.length} at locus {self.locus}>'
    
class protein_obj:
    '''Object with information about a full protein.
        locus: Locus tag of the protein.
        length: Total length of the protein.
        domains: List with presence/absence of signal peptide, no. of GH domains
        and no. of glucan-binding and cell wall-binding domains.'''
    def __init__(self, locus, length, domains = []):
        self.locus = locus
        self.length = int(length)
        self.domains = domains
    def __str__(self): #Object information that is printed
        return f'<{len(self.domains)}-domain protein of length {self.length} at locus {self.locus}>'

# =============================================================================
# 2. Define inputs and formatting
# =============================================================================
types = {'SP': ['Signal peptide', 'DI', '#00197F'], #Dictionary with formatting information to plot the domains 
         'GH70': ['Glycosyl hydrolase family 70', 'EL', '#FF006F'], 
         'DUF': ['DUF5776', 'TR', '#7D8292'],
         'GH32': ['Glycosyl hydrolase family 32', 'EL', '#EE3900'], 
         'CB': ['Choline-binding', 'RE', '#FFAB00'], 
         'CWB': ['Cell-wall binding', 'RE', '#FF8900']}

short_genes = ['A0901_14250', 'A1001_13210', 'A1003_13400', 'A1202_14500',
               'A1401_13710', 'A1404_14290', 'A1805_13670', 'K2W83_RS06625',
               'FHON2_14510', 'G0101_13740', 'G0403_13980', 'H1B104J_13940',
               'H1B105A_13110', 'H1B302M_13810', 'H3B101A_14180', 
               'H3B104J_13910', 'H3B104X_14110', 'H3B202X_13780', 
               'H3B203J_14310', 'H3B203M_13420', 'H3B206M_13730', 
               'H3B209X_14340', 'H4B111J_14460', 'H4B202J_13700', 
               'H4B204J_14250', 'H4B205J_13920', 'H4B206J_14400', 
               'H4B211M_13940', 'H4B402J_13530', 'H4B405J_14150',
               'H4B406M_14380', 'H4B412M_14170', 'H4B501J_13790 ',
               'H4B503X_13610', 'H4B504J_14270', 'H4B505J_13820', 
               'LDX55_06780', 'APS55_RS03400']

rep_strains = ['A0901', 'A1001', 'A1003', 'A1202', 'A1401', 'A1404', 'A1805', 
               'DSMZ12361', 'fhon2', 'G0101', 'G0403', 'H1B1-04J', 'H1B1-05A', 
               'H1B3-02M', 'H3B1-01A', 'H3B1-04J', 'H3B1-04X', 'H3B2-02X', 
               'H3B2-03J', 'H3B2-03M', 'H3B2-06M', 'H3B2-09X', 'H4B1-11J', 
               'H4B2-02J', 'H4B2-04J', 'H4B2-05J', 'H4B2-06J', 'H4B2-11M', 
               'H4B4-02J', 'H4B4-05J', 'H4B4-06M', 'H4B4-12M', 'H4B5-01J', 
               'H4B5-03X', 'H4B5-04J', 'H4B5-05J', 'IBH001', 'MP2']

workdir = os.path.expanduser('~') + '/GH_project' #Working directory
indir = f'{workdir}/interproscan' #Directory with inputs
outdir = f'{workdir}/data/tabs'

GS3 = ['H3B206M_12830', 'H4B111J_13560', 'H4B405J_13350', 'H4B406M_13450'] #List of GS3 genes
GH_types = ['GH70', 'GH32'] #List of glycosyl hydrolase domains

short_tree = f'{workdir}/data/fasta/GH70/trees/complete_short_repset.mafft.faa.treefile'

signalp_dir = f'{workdir}/signalp_pred/slow'

if not os.path.exists(outdir): #If the output directory doesn't exist
    os.makedirs(outdir) #Create it

# =============================================================================
# 3. Save information to file
# =============================================================================
signalp_dict = {}
signalp = {}
ref_tree1 = f'{workdir}/data/fasta/GH70/trees/GH70_functional_repset.mafft.faa.treefile' #Load reference tree of GH70 genes
ref_tree2 = f'{workdir}/data/fasta/GH32/trees/GH32_repset.mafft.faa.treefile' #Set the right name of the input treefile

dom_out1 = f'{outdir}/GH70_domain_file.txt' #Set the path to the output file

for GH_type in GH_types:
    signalp_inputs = f'{signalp_dir}/{GH_type}'
    for filename in os.listdir(signalp_inputs):
        if filename.endswith('plot.txt'):
            with open(f'{signalp_inputs}/{filename}') as sp:
                locus = filename.replace('output_', '').replace('_plot.txt', '')
                df = pd.read_csv(sp, sep = '\t', skiprows = 1)
                for index, row in df.iterrows():
                    if row['Other'] > 0.6:
                        signalp_dict[locus] = row['# pos']
                        break
                
    ref_tree =  f'{workdir}/data/fasta/{GH_type}/trees/{GH_type}_functional_repset.mafft.faa.treefile' #Load reference tree of GH70 genes
    dom_out = f'{outdir}/{GH_type}_domain_file.txt' #Set the path to the output file
    
    if GH_type == 'GH32':
        ref_tree = ref_tree.replace('_functional', '')
        
    proteins = {} #Create an empty dictionary to store protein information
    locus_list = [] #Create empty list of locus tags
    
    with open(ref_tree) as ref: #Open the tree with all domains as a reference
        ref_text = ref.read() #Read tree
        print(ref_text) #Print text in the treefile

    infiles = [filename for filename in os.listdir(indir) if filename.split('.')[0] in rep_strains]
    for filename in infiles: #Loop through files in input directory that end with .tsv
        with open(f'{indir}/{filename}') as file: #Open file
            domain_tab = csv.reader(file, delimiter = '\t') #Read file
            for line in domain_tab: #Loop through lines in file
                locus = line[0].replace('-', '').upper() #Remove - character from locus tags in the files
                
                if locus == 'MP2_13360': #If the locus is the GS2 from MP2
                    locus = 'APS55_RS03845' #Manually assign the right locus tag to it
                elif locus == 'MP2_13350': #If the locus is the GS1 from MP2
                    locus = 'APS55_RS03850' #Do the same
                elif locus == 'MP2_14250':
                    locus = 'APS55_RS03400'
                
                signalp[locus] = False
                    
                add = False #Boolean to check if the domain should be retrieved
                if locus not in locus_list and locus in ref_text:
                    proteins[locus] = protein_obj(locus, line[2], []) #Create a protein object for the locus
                    locus_list.append(locus) #Append locus to the list of loci
                    if locus in signalp_dict and signalp_dict[locus] > 1 and signalp[locus] == False:
                        signalp[locus] = True
                        sp = domain_obj(locus, 1, signalp_dict[locus], 'SP', types['SP'][1], types['SP'][2]) #Create a domain object with domain information
                        if sp not in proteins[locus].domains and sp.locus == locus: #If the domain is not in the list of domains for the protein
                            proteins[locus].domains.append(sp) #Add domain to the list
                    
                if locus in ref_text:
                    #Condition to get any of the domains of interest
                    if GH_type == 'GH32' and ((line[3] == 'Pfam' and 'Choline' in line[5]) or (line[3] == 'SMART' and 'glyco_32' in line[5])):
                        if line[3] == 'SMART': #If the method to retrieve the domain was SMART
                            dom_type = 'GH32' #The domain type is set to GH32
                            add = True 
                        elif line[3] == 'Pfam': #If the method to retrieve the domain was Pfam
                            if 'Choline' in line[5]: #If it is a choline-binding domain
                                dom_type = 'CB' #Set the gene type to CB
                                add = True
                            elif 'cell wall' in line[5]:
                                dom_type = 'CWB'
                                add = True
                            
                    #Condition to get the domains of interest
                    elif GH_type == 'GH70' and (locus in short_genes) or (line[3] == 'Pfam'):
                        #Condition to retrieve signal peptides (prioritize TIGRFAM)
                        if 'family 70' in line[5]: #If the domain is a GH70 domain
                            dom_type = 'GH70' #Save it as GH70
                            add = True
                        elif 'Choline' in line[5]:
                            dom_type = 'CB'
                            add = True
                        elif 'cell wall' in line[5]: #If the domain is annotated as a cell wall-binding domain
                            dom_type = 'CWB' #Save it as CWB
                            add = True
                        elif 'DUF' in line[5]: #If it is a domain of unknown function
                            dom_type = 'DUF' #Save it as DUF
                            add = True
                        
                    if add: #If the domain should be added
                        print(locus, dom_type)
                        domain = domain_obj(locus, line[6], line[7], dom_type, types[dom_type][1], types[dom_type][2]) #Create a domain object with domain information
                        if domain not in proteins[locus].domains and domain.locus == locus: #If the domain is not in the list of domains for the protein
                            proteins[locus].domains.append(domain) #Add domain to the list
    
    with open(dom_out, 'w') as dom_file: #Open output file in write mode
        #Write file headers    
        dom_file.write('DATASET_DOMAINS\n')
        dom_file.write('SEPARATOR COMMA\n')
        dom_file.write(f'DATASET_LABEL,{GH_type}\n')
        dom_file.write('COLOR,#ff0000\n')
        dom_file.write('DATA\n')
        for locus in locus_list:
            line = f'{proteins[locus].locus},{proteins[locus].length},'

            for domain in proteins[locus].domains: #Loop through domains in protein
                line += f'{domain.shape}|{domain.start}|{domain.end}|' #Add domain formatting information to file
                line += f'{domain.color.lower()}|{domain.type},'
            line = line[:-1] + '\n' #Add a line break at the end of the line
            dom_file.write(line) #Write line to file
                