# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 16:04:16 2021

Script to retrieve the plotting information for the CDS of the genes
in a region of interest so that they can be later plotted next to
the tree in ete3.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import Python libraries
# =============================================================================
import os
import csv
import pathlib

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

workdir = os.path.expanduser('~') + '/GH_project' #Working directory
indir = f'{workdir}/interproscan' #Directory with inputs
outdir = f'{workdir}/data/tabs'

GS3 = ['H3B206M_12830', 'H4B111J_13560', 'H4B405J_13350', 'H4B406M_13450'] #List of GS3 genes
GH_types = ['GH70', 'GH32'] #List of glycosyl hydrolase domains

short_tree = f'{workdir}/data/fasta/GH70/trees/complete_short_repset.mafft.faa.treefile'

if not os.path.exists(outdir): #If the output directory doesn't exist
    os.makedirs(outdir) #Create it

# =============================================================================
# 3. Save information to file
# =============================================================================
for GH_type in GH_types: #Loop through GH types
    ref_tree = f'{workdir}/data/fasta/{GH_type}/trees/{GH_type}_functional_repset.mafft.faa.treefile' #Load reference tree of GH70 genes
    if GH_type == 'GH32': #If the genes have GH32 domains
        ref_tree = ref_tree.replace('functional_', '') #Set the right name of the input treefile
    dom_out = f'{outdir}/{GH_type}_domain_file.txt' #Set the path to the output file
    
    proteins = {} #Create an empty dictionary to store protein information
    locus_list = [] #Create empty list of locus tags
        
    with open(ref_tree) as ref: #Open the tree with all domains as a reference
        ref_text = ref.read() #Read tree
        print(ref_text) #Print text in the treefile
    
    for filename in pathlib.Path(indir).glob('*.tsv'): #Loop through files in input directory that end with .tsv
        with open(filename) as file: #Open file
            domain_tab = csv.reader(file, delimiter = '\t') #Read file
            for line in domain_tab: #Loop through lines in file
                sig_pep = False #Boolean to check if a signal peptide is predicted by any method
                locus = line[0].replace('-', '').upper() #Remove - character from locus tags in the files
                if locus == 'MP2_13360': #If the locus is the GS2 from MP2
                    locus = 'APS55_RS03845' #Manually assign the right locus tag to it
                elif locus == 'MP2_13350': #If the locus is the GS1 from MP2
                    locus = 'APS55_RS03850' #Do the same
                elif locus == 'MP2_14250':
                    locus = 'APS55_RS03400'
                
                if locus in ref_text or locus in short_genes: #If the locus tag is in the reference tree
                    print(locus) #Print the locus
                    add = False #Boolean to check if the domain should be retrieved
                    #Condition to get any of the domains of interest
                    if ('GH32' in ref_tree) and ((line[3] == 'Pfam' and 'Choline' in line[5]) or ('A1404' in locus and line[3] == 'Phobius' and line[4] == 'SIGNAL_PEPTIDE') or (line[3] == 'Pfam' and 'signal peptide' in line[5] and sig_pep == False) or (line[3] == 'TIGRFAM' and 'signal peptide' in line[5] and sig_pep == False) or (line[3] == 'SMART' and 'glyco_32' in line[5])):
                        if line[3] == 'TIGRFAM': #If the method to detect the domain is TIGRFAM
                            dom_type = 'SP' #Assign it as a signal peptide domain
                            sig_pep = True #Set the boolean for SP presence to true
                            add = True #Set the boolean to retrieve the domain to true
                        elif line[3] == 'SMART': #If the method to retrieve the domain was SMART
                            dom_type = 'GH32' #The domain type is set to GH32
                            add = True 
                        elif line[3] == 'Pfam': #If the method to retrieve the domain was Pfam
                            if 'Choline' in line[5]: #If it is a choline-binding domain
                                dom_type = 'CB' #Set the gene type to CB
                                add = True
                            elif 'cell wall' in line[5]:
                                dom_type = 'CWB'
                                add = True
                            else: #If no signal peptide has been retrieved
                                dom_type = 'SP' #Retrieve the signal peptide
                                sig_pep = True 
                                add = True
                        elif line[3] == 'Phobius': #If the method was Phobius
                            dom_type = 'SP' #Set the gene type to signal peptide
                            sig_pep = True
                            add = True
                            
                    #Condition to get the domains of interest
                    elif ('GH70' in ref_tree or locus in short_genes) and (line[3] == 'Pfam' or (line[3] == 'TIGRFAM' and 'signal peptide' in line[5]) or (locus in GS3 and line[3] == 'Phobius' and line[4] == 'SIGNAL_PEPTIDE')):
                        #Condition to retrieve signal peptides (prioritize TIGRFAM)
                        if (line[3] == 'TIGRFAM' or ('LDX55' not in locus and 'K2W83' not in locus and line[3] == 'Pfam' and 'signal peptide' in line[5]) or (locus in GS3 and line[3] == 'Phobius' and line[4] == 'SIGNAL_PEPTIDE')) and sig_pep == False:
                            dom_type = 'SP'
                            sig_pep = True
                            add = True
                        elif 'family 70' in line[5]: #If the domain is a GH70 domain
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
                        domain = domain_obj(locus, line[6], line[7], dom_type, types[dom_type][1], types[dom_type][2]) #Create a domain object with domain information
                        if locus not in locus_list: #If the locus has not been recorded yet
                            protein = protein_obj(locus, line[2], [domain]) #Create a protein object for the locus
                            proteins[locus] = protein #Save the object to a dictionary where the locus is the key and the object is the value
                            locus_list.append(locus) #Append locus to the list of loci
                        if domain not in proteins[locus].domains and domain.locus == locus: #If the domain is not in the list of domains for the protein
                            proteins[locus].domains.append(domain) #Add domain to the list
        
    with open(dom_out, 'w') as dom_file: #Open output file in write mode
        #Write file headers    
        dom_file.write('DATASET_DOMAINS\n')
        dom_file.write('SEPARATOR COMMA\n')
        dom_file.write('DATASET_LABEL,gtf\n')
        dom_file.write('COLOR,#ff0000\n')
        dom_file.write('DATA\n')
        
        for l in proteins.keys(): #Loop through proteins in the dictionary
            length = proteins[l].length #Get the length of the protein
            line = f'{l},{length},' #Add locus and length information to file
            for domain in proteins[l].domains: #Loop through domains in protein
                line += f'{domain.shape}|{domain.start}|{domain.end}|' #Add domain formatting information to file
                line += f'{domain.color.lower()}|{domain.type},'
            line = line[:-1] + '\n' #Add a line break at the end of the line
            dom_file.write(line) #Write line to file
            