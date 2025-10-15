# GH project
Family 70 and 32 glycoside hydrolases in _Apilactobacillus kunkeei_.

## Folder organization
Global folder organization:
The `code` folder contains most of the code for this project, except for:
- ``RDP5_analysis/code`` -> Code to generate alignments and plots for RDP5
- ``GS2_BRS_repeats/code`` -> Code to get an alignment of several GS2 and BRS and plot it
- ``motifs/code`` -> Code to retrieve motifs from GH70 genes in _A. kunkeei_
- ``add_species/code`` -> Code to generate a GH70 tree adding sequences from other species and to retrieve motifs from the added sequences

## How to use this code
Apart for the organization into subfolders, the scripts are numbered in the order in which they should be ran.
When there are two numbers in the scripts (for example, 4.4), it means that the script relies on any previous steps labelled as 1-3.4.
If no scripts before 4.4 are labelled as 1-3.4, it indicates that the script depends on the main pipeline scripts 1-3, and the following steps should be named 5.4-N.4.
Each script indicates a description of what it does in the beginning of the file.
The main pipeline is included in a Snakemake pipeline (see ``Snakefile``).
The scripts for core genes and core protein trees have their own pipeline (see ``code/core_genes`` and ``core_Snakefile``).

