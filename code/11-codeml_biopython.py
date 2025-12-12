# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 10:30:36 2021

This is a script to run CodeML and calculate pairwise substitution rates of 
GH70 and GH32 genes. 

@author: Marina Mota Merlo
"""

import os, logging, traceback

##CodeML in BioPython

import Bio
import Bio.Phylo.PAML.codeml as codeml

# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(filename = snakemake.log[0], level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                          *traceback.format_exception(exc_type, exc_value, exc_traceback)
                          ]))

sys.excepthook = handle_exception

sys.stdout = open(snakemake.log[0], 'a')

# =============================================================================
# Define CodeML objects
# =============================================================================

def codeml_object(workdir, seqfile, outfile, treefile = 'trees/empty.tree'):
    cml = codeml.Codeml()            #Initialize object
    cml.alignment = '/home/marina/GH_project/' + str(seqfile)     #Assign alignmentfile to object
    cml.tree = str(treefile)         #Assign tree (optional)
    cml.out_file = '/home/marina/GH_project/' + str(outfile)      #Set output name
    cml.working_dir = str(workdir)   #Revise this after setting up snakemake
    
    return cml

# =============================================================================
# Set CodeML options (same role as control file codeml.ctl)
# =============================================================================

def set_codeml_options(cml, rmd, seq, mdl = 1, noise = 0, verb = False, cf = 2, icd = 0):
    cml.set_options(noisy = noise)    #Revise this option
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
    cml.set_options(omega = 0.7)      #Initial omega
    
    return cml

# =============================================================================
# Apply these functions to snakemake pipeline
# =============================================================================

tree_inputs = list(snakemake.input['trees_GH70']) + list(snakemake.input['trees_GH32'])

for i in range(len(tree_inputs)):
    cml = codeml_object(snakemake.params.outdir + '/' + os.path.basename(tree_inputs[i]).split('_')[0],
                        snakemake.input['codons'][i], snakemake.output['summary'][i])
    cml = set_codeml_options(cml, -2, 1)
    
    ##Run CodeML
    cml.run(verbose = True, parse = False)
