#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 13:35:39 2025

This script converts the protein alignments and DNA data into codon alignments
and then runs CodeML on said alignments.

It requires Pal2Nal in the working directory (or to change the path to Pal2nal)
and CodeML in Biopython.

@author: Marina Mota-Merlo
"""

# =============================================================================
# 0. Import required packages
# =============================================================================

import subprocess
import os
import Bio.Phylo.PAML.codeml as codeml

# =============================================================================
# 1. Define inputs and paths
# =============================================================================

workdir = os.path.expanduser('~') + '/GH_project' #working directory
inpath = f'{workdir}/all_core/38_strains/fasta' #Path to input files
in_tree = f'{workdir}/trees/empty.tree'
outpath = inpath.replace('fasta', 'codons') #Path to Pal2Nal output
codeml_dir = inpath.replace('fasta', 'results') #Path to CodeML outputs

if not os.path.exists(outpath): #Create output directory if it does not exist
    os.makedirs(outpath)

infiles = [file for file in os.listdir(inpath) if file.endswith('.mafft.faa')] #Get list of input protein alignments

# =============================================================================
# 2. Define CodeML functions
# =============================================================================

def codeml_object(workdir, inpath, seqfile, outfile, treefile = in_tree):
    '''Function to create CodeML objects.
    Inputs: workdir  - working directory
            seqfile  - input codon alignment
            outfile  - output txt file with CodeML log
            treefile - phylogeny on which the calculations will be based (by
            default, an empty tree)'''
    cml = codeml.Codeml()                  #Initialize object
    cml.alignment = f'{inpath}/{seqfile}' #Assign alignmentfile to object
    cml.tree = str(treefile)               #Assign tree (optional)
    cml.out_file = str(outfile)            #Set output name
    cml.working_dir = str(workdir)         #CodeML working directory
    
    return cml

def set_codeml_options(cml, rmd, seq, mdl = 1, noise = 0, verb = False, cf = 2, icd = 0):
    '''Function to set the CodeML running options for a CodeML object.
    Inputs: cml   - CodeML object
            rmd   - Runmode
            seq   - Sequence type, codons (1) or amino acids (2)
            mdl   - Number of omega ratios per branch, by default 1
            noise - How much of the information is printed on-screen
            verb  - How detailed the output is, by default, less detailed
            cf    - Codon frequency model
            icd   - Translation code, universal by default'''
    cml.set_options(noisy = noise)    #How much info is printed on screen
    cml.set_options(verbose = verb)   #False: Do not show run output on-screen
    cml.set_options(runmode = rmd)    #-2: Pairwise comparisons
    cml.set_options(seqtype = seq)    #1: Codons, 2: AA
    cml.set_options(CodonFreq = cf)   #2: Average nt sequences at the 3 codon positions
    cml.set_options(model = mdl)      #1: One omega ration per branch
    cml.set_options(ncatG = 8)        #Number of gamma categories
    cml.set_options(Small_Diff = 1e-6)#Change per iteration
    cml.set_options(icode = icd)      #0: Universal translation code (GeneBank)
    cml.set_options(fix_kappa = 0)    #0: Kappa not fixed
    cml.set_options(kappa = 2)        #Value of initial/fixed kappa
    cml.set_options(fix_omega = 0)    #Estimate omega
    cml.set_options(omega = 0.7)      #Initial omega (check)
    
    return cml

# =============================================================================
# 3. Run Pal2nal and CodeML
# =============================================================================

for aln in infiles: #Loop through fasta protein alignments
    fna = aln.replace('.mafft.faa', '.fna') #Name of the fna file
    outfile = aln.replace('.mafft.faa', '.pal2nal') #Name of the Pal2nal output
    subprocess.run(f'{workdir}/pal2nal.v14/pal2nal.pl {inpath}/{aln} {inpath}/{fna} -output paml -nogap > {outpath}/{outfile};',
                   shell = True) #Run Pal2nal
    
    codeml_out = outfile.replace('.pal2nal', '.txt') #Name of the CodeML output
    
    outdir = f'{codeml_dir}/{codeml_out.replace(".txt", "")}' #Path to outputs
    
    if not os.path.exists(outdir): #If the output path doesn't exist
        os.makedirs(outdir) #Create it
    
    cml = codeml_object(outdir, outpath, outfile, codeml_out) #Create CodeML object
    cml = set_codeml_options(cml, -2, 1) #Set run options
    cml.run(verbose = True, parse = False) #Run CodeML
