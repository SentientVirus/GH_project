# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 10:30:36 2021

This is a script to use CodeML through BioPython to calculate the
pairwise substitution rates between the cores genes in the different strains.

@author: Marina Mota Merlo
"""

# =============================================================================
# 0. Import required modules
# =============================================================================

import os, logging, traceback
import Bio.Phylo.PAML.codeml as codeml

# =============================================================================
# 0. Logging
# =============================================================================
print('The script is running...')
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

print('Start logging...')

# =============================================================================
# 1. Define CodeML objects
# =============================================================================

def codeml_object(workdir, seqfile, outfile, treefile = 'trees/empty.tree'):
    cml = codeml.Codeml()            #Initialize object
    cml.alignment = str(seqfile)     #Assign alignmentfile to object
    cml.tree = str(treefile)         #Assign tree (optional)
    cml.out_file = str(outfile)      #Set output name
    cml.working_dir = str(workdir)   #Revise this after setting up snakemake
    
    return cml

# =============================================================================
# 2. Set CodeML options (same role as control file codeml.ctl)
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
    cml.set_options(omega = 0.7)      #Initial omega (check)
    
    return cml

# =============================================================================
# 3. Apply these functions to the Snakemake pipeline
# =============================================================================

#Create output directory if it doesn't exist
if not os.path.exists(snakemake.output.outdir):
    os.makedirs(snakemake.output.outdir)

for i in range(len(snakemake.output.outfiles)): #Loop through output files
    sub_outdir = os.path.dirname(snakemake.output.outfiles[i]) #Retrieve the path to the output file
    print(sub_outdir)

    #Create the path to the output file if it doesn't exist
    if not os.path.exists(sub_outdir):
        os.makedirs(sub_outdir)

    #Create a CodeML object and set the CodeML parameters
    cml = codeml_object(sub_outdir, snakemake.input[i],
                        snakemake.output.outfiles[i], snakemake.input[-1])
    cml = set_codeml_options(cml, -2, 1)
    print(cml.alignment)
    
    #Run CodeML
    cml.run(verbose = True, parse = False)
