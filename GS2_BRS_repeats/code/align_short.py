#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:31:16 2024

@author: Marina Mota-Merlo

This script retrieves the sequences of GS2_BRS genes and those of G0403 GS2 and
BRS and A1401 BRS and performs an alignment of all the sequences. Requires
Mafft and PyMSAViz
"""

import logging, traceback, sys
import os
import subprocess
from pymsaviz import MsaViz

# =============================================================================
# 0. Logging
# =============================================================================
log = os.path.expanduser('~') + '/GH_project/GS2_BRS_repeats/logs/align_short.log'

logger = logging.getLogger()

logging.basicConfig(filename = log, level = logging.INFO,
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

sys.stdout = open(log, 'a')

# =============================================================================
# 1. Defining inputs and functions
# =============================================================================
workdir = os.path.expanduser('~') + '/GH_project/GS2_BRS_repeats/files'
in_faa = f'{workdir}/DSM_seqs.faa'
out_mafft = in_faa.replace('.faa', '.mafft.faa')
out_aln_fig = os.path.expanduser('~') + '/GH_project/GS2_BRS_repeats/plots/DSM_PTL_aln.png'

indir = os.path.dirname(in_faa)
log_path = os.path.dirname(log)
fig_path = os.path.dirname(out_aln_fig)
paths = [indir, log_path, fig_path]

threads = 12

[os.makedirs(file_path) for file_path in paths if not os.path.exists(file_path)]
    
with open(log, 'w+') as flog:
    flog.write('')
            
# =============================================================================
# 2. Run Mafft-linsi and save results
# =============================================================================
with open(in_faa) as fasta_seqs:
    subprocess.run(f'mafft-linsi --thread {threads} {in_faa} > {out_mafft} 2>> {log};', shell = True)
    
# =============================================================================
# 3. Generate plot
# =============================================================================
mv = MsaViz(out_mafft, color_scheme = 'Flower', start = 1, end = 400, wrap_length=100)#, 
            # show_count = True, show_consensus = True)

mv.add_text_annotation((6, 46), 'Signal peptide', text_size = 12,
                        text_color = 'black', range_color = 'black')
mv.add_text_annotation((287, 400), 'GH70 domain', text_size = 12,
                        text_color = 'black', range_color = 'black')
mv.savefig(out_aln_fig)
mv.savefig(out_aln_fig.replace('png', 'svg'))
mv.savefig(out_aln_fig.replace('png', 'tiff'))
