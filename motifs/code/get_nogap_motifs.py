# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:39:34 2021

@author: usuario
"""
import re
from Bio import SeqIO

class motif:
    def __init__(self, num, seq, start, end): #Define motif class
        self.num = num #Motif number (1-7)
        self.seq = seq #Motif sequence
        self.start = start #Start position of the motif in a protein sequence
        self.end = end #End position of the motif in the sequence
    def __str__(self): #What print(class) returns
        return f'<Motif {self.num} from position {self.start} to {self.end}>'

class prot_sequence: #All motifs associated to protein (1 sequence/position per motif)
    def __init__(self, locus, I, II, III, IV, V, VI, VII):
        self.locus = locus #Locus tag
        self.I = I #Domain 1
        self.II = II #Domain 2
        self.III = III
        self.IV = IV
        self.V = V
        self.VI = VI
        self.VII = VII #Domain 7
        self.n_motifs = 0 #Number of motifs
        self.bool_list = [0, 0, 0, 0, 0, 0, 0, 0]
    def calculate_n_motifs(self): #Function to calculate number of motifs in the protein
        motif_list = [self.I, self.II, self.III, self.IV, self.V, self.VI, self.VII]
        bool_list = [0 + len(motif) for motif in motif_list]
        self.bool_list = [0] + bool_list
        self.n_motifs = sum(bool_list)
    def __str__(self): #What print(class) returns
        return f'<{self.locus} with {self.n_motifs} motifs>'
    
def search_motif(seq_str, num): #Function to search for any motif
    motif_list = ['', r'[MAVT]D(.{1,1})V(.{1,1})[ND]Q', r'(.{1,1})R(.{1,1})DA(.{2,2})[NFS][MVI][DHNS]',
                  r'H[LI][SHNA][LI][LV]E(.*?)(([WGP](.{3,3})[DGTS])|(SLEAAT))',
                  r'[FIY][VIL][RHN][AV]HD(.{1,2})[AVIS]Q(.{2,2})[IV]',
                  r'[EDA][FL]LL[AG][NSD][DQ][IVE]DNSN[PV]', 
                  r'[LW]G[IVF](.{3,3})[EQ][FLM][AP][PG][HQ]Y',
                  r'[VTS][VTI]PR[VMI][YF]YGD']
    motifs = []

    motif_info = re.finditer(motif_list[num], seq_str)
    for motif_match in motif_info:
        if motif_match is not None: #If there is a match, store new values
            motif_seq = motif_match.group(0)
            start = motif_match.start() + 1
            end = motif_match.end() + 1
        else: #Otherwise use generic values
            motif_seq = ''
            start = 0
            end = 0
        motifs.append([num, motif_seq, start, end])
    return motifs

def get_group(locus):
    GS1 = ['A1001_12310', 'A1202_13520', 'A1401_12750', 'A1805_12820', 'fhon2_13540', 
	       'G0101_12800', 'G0102_12710','G0102_12720', 'G0403_13100', 'H1B1-04J_13010',
               'H1B1-05A_12300', 'H1B3-02M_12880', 'H2B1-05J_12840', 'H2B1-05J_12850',
               'H3B1-01A_13240', 'H3B1-01J_12980', 'H3B1-02X_13350', 'H3B1-03M_13200', 
               'H3B1-11A_13170', 'H3B2-02X_12850', 'H3B2-03M_12480', 'H3B2-09X_13340', 
               'H4B1-01A_12830', 'H4B2-02J_12880', 'H4B2-04J_13330', 'H4B4-02J_12600', 
               'H4B4-12M_13240', 'H4B5-03X_12670', 'H4B5-04J_13460', 'H4B5-05J_12880', 
               'MP2_13350']
    
    GS2_BRS = ['H3B1-01J_13010','H4B2-02J_12890', 'H4B2-04J_13350']
    
    GS = ['A1003_12540']

    GS2 = ['fhon2_13560', 'fhon2_13570', 'G0403_13120', 'H1B3-02M_12900', 'H3B1-01A_13260',
               'H3B1-03M_13220', 'H3B1-11A_13190', 'H3B2-02X_12860', 'H3B2-09X_13360',
               'H3B2-09X_13370', 'H4B5-04J_13480', 'H4B5-05J_12900', 'MP2_13360']

    NCB = ['A1401_04720', 'fhon2_04830', 'H3B1-03M_04730', 'A1202_05530',
               'H3B1-01A_04720', 'H3B1-02A_04730', 'A1805_04920', 'H4B5-04J_05510',
               'H3B2-02X_04770', 'H4B2-04J_04710', 'A1003_04750', 'H1B3-02M_04960',
               'G0101_04800', 'H1B1-05A_04750', 'A0901_05380', 'H3B2-03M_04710',
               'H3B1-02X_04940', 'H4B5-05J_04980', 'H4B1-01A_04950']

    BRS = ['A1401_12760', 'fhon2_13550', 'G0403_13110', 'H1B3-02M_12890', 'H3B1-01A_13250',
           'H3B1-01J_12990', 'H3B1-03M_13210', 'H3B1-11A_13180', 'H3B1-11M_12570', 
           'H3B1-11M_12580', 'H3B2-09X_13350', 'H4B2-04J_13340', 'H4B5-04J_13470',
           'H4B5-05J_12890']

    short = ['A1202_14500', 'A1404_14290', 'A1805_13670', 'H4B2-04J_14250', 'H2B1-05J_13740', 'H4B1-01A_13720',
             'A0901_14250', 'A1001_13210', 'G0403_13980', 'fhon2_14510',
             'A1003_13400', 'H1B1-04J_13940', 'H1B3-02M_13810', 'H4B4-02J_13530', 
             'H4B2-02J_13700', 'H3B1-01A_14180', 'H3B1-02X_14140', 'H3B1-01J_13900',
             'H3B1-03M_14110', 'H3B1-03M_14110']
    if locus in GS1:
        return 'GS1'
    elif locus in GS2:
        return 'GS2'
    elif locus in NCB:
        return 'NCB'
    elif locus in BRS:
        return 'BRS'
    elif locus in short:
        return 'short'
    elif locus in GS2_BRS:
        return 'GS2_BRS'
    elif locus in GS:
        return 'GS'

prot_list = []
GH_type = 'BRS'
input_file = f'../PoSE/{GH_type}_codon_nogap.fasta' #'a_kunkeei_GH70_codon_nogap.fasta' #'a_kunkeei_GH70_prot.fasta' #'a_kunkeei_gtf_prot.fasta'
output_file = f'{GH_type}_nogap_motif_presence.tsv' #'GH70_nogap_motif_presence.tsv' #'GH70_motif_presence.tsv' #'motif_presence.tsv'
output_file2 = f'{GH_type}_nogap_motifs.tsv' #'GH70_nogap_motifs.tsv' #'GH70_motifs.tsv' #'gtf_motifs.tsv'
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
with open(output_file, 'w') as out_file, open(output_file2, 'w') as out_file2:
    out_file.write('locus\ttype\tmotif_I\tmotif_II\tmotif_III\tmotif_IV\tmotif_V\tmotif_VI\tmotif_VII\n')
    out_file2.write('locus\ttype\tmotif_type\tsequence\tstart\tend\ttotal_length\n')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        category = get_group(name)
        if category == None:
            category = get_group(name[:-2])
        #print(sequence)
        motif_list = ['']
        for i in range(1, 8):
            motif_matches = search_motif(sequence, i)
            motif_list.append(motif_matches)
            for mot in motif_matches:
                mot = motif(*mot)
                if len(mot.seq) > 0:
                   out_file2.write(f'{name}\t{category}\t{mot.num}\t{mot.seq}\t{mot.start}\t{mot.end}\t{len(sequence)}\n')
        prot = prot_sequence(name, motif_list[1], motif_list[2], motif_list[3], motif_list[4], motif_list[5], motif_list[6], motif_list[7])
        prot.calculate_n_motifs()
        prot_list.append(prot)
        bool_list = list(map(int, prot.bool_list))
        out_file.write(f'{prot.locus}\t{category}\t{bool_list[1]}\t{bool_list[2]}\t'+
                       f'{bool_list[3]}\t{bool_list[4]}\t{bool_list[5]}\t'+
                       f'{bool_list[6]}\t{bool_list[7]}\n')

'''VI_motif = re.search(r'[LW]G[IVF][SDTY]((NI)|S)[FLV][EQ][LM][AP]PQY', string)
start = VI_motif.start()+1
end = VI_motif.end()+1
view = VI_motif.group(0)

I_motif = re.search(r'[MAV]D[FIY]VPDQ', 'LHAQGIKVLADIVPDQVYSLPNQ')
I_view = I_motif.group(0)
V_motif = re.search(r'[EDA][FL]LL[AG][NSD][DQ][IVE]DNSN[PV]', 'YKNYELLLANDIDNSNPTVQAEYKNYELLLANDIDNSNPTVQAE')
V_view = V_motif.group(0)
II_motif = re.search(r'[MIV]R[VIL]DA[AVILP][DS][NF][MVI][DHNS]', 'DGVRLDALDNVDADVVNI')
II_view = II_motif.group(0)
III_motif = re.search(r'H[LI][SHNA][LI][LV]E(.*?)(([WGP](.{3,3})[DGTS])|(SLEAAT))', 'ANANKHLSILEDWSKNDAYYQ')
III_view = III_motif.group(0)
IV_motif = re.search(r'[FIY][VI][RH]AHD(.{2,2})[AVIS]Q(.{2,2})[IV]', 'GGMPNYSYVHAHDKGIQERVGQAIVDT')
IV_view = IV_motif.group(0)
VII_motif = re.search(r'[VTS][VTI]PR[VMI][YF]YGD', 'A')
VII_view = VII_motif.group(0)'''
