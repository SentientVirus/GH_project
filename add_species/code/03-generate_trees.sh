threads=16
basedir=~/GH_project/add_species
parentdir=~/GH_project
in_og=$basedir"/results/domains/ingroup_domains.faa"
in_GH70=$parentdir"/data/fasta/GH70/GH70_functional_repset.faa" #"files/GH70_30kb_small_subset.faa"  #"../data/fasta/GH70/GH70_functional_repset.faa"
in_fasta=$basedir"/results/domains/GH70_species.faa"
alignment=$basedir"/results/alignment/GH70_species.mafft.faa"
aln_log=$basedir"/logs/03.1-alignment.log"
tree_log=$basedir"/logs/03.2-tree.log"

mkdir -p $basedir/logs
mkdir -p $(dirname -- $alignment)

# Put outgroups and GH70s from A. kunkeei in the same file
cat $in_og $in_GH70 > $in_fasta

# The following part required: conda activate codeml-run

# Run alignment
mafft-linsi --thread $threads $in_fasta > $alignment 2> $aln_log

# Create tree
iqtree --redo -nt $threads -s $alignment -st AA -m TEST -bb 1000 -bnni > $tree_log 2> $tree_log

# Former model -m LG+G4+F
