threads=16
in_og="results/domains/ingroup_domains.fasta"
in_GH70="files/GH70_30kb_repset.faa"  #"../data/fasta/GH70/GH70_functional_repset.faa"
in_fasta="results/domains/GH70_species.faa"
alignment="results/alignment/GH70_species.mafft.faa"
annot_log="logs/01-interproscan.log"
aln_log="logs/03.1-alignment.log"
tree_log="logs/03.2-tree.log"

mkdir -p logs
mkdir -p $(dirname -- $alignment)

# Put outgroups and GH70s from A. kunkeei in the same file
cat $in_og $in_GH70 > $in_fasta

# The following part required: conda activate codeml-run

# Run alignment
mafft-linsi --thread $threads $in_fasta > $alignment 2> $aln_log

# Create tree
iqtree -nt $threads -s $alignment -st AA -m LG+G4+F -bb 1000 -bnni 2> $tree_log

