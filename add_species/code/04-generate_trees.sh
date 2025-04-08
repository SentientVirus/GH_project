#Environment: alignment_tree.yml
threads=48 #No. of threads
basedir=~/GH_project/add_species #Working directory
parentdir=~/GH_project #Parent directory of the working directory
in_og=$basedir"/results/domains/blastp_domains.faa" #Input fasta file with GH70 domains from other species
in_GH70=$parentdir"/data/fasta/GH70/GH70_functional_repset.faa" #Input fasta file with GH70 domains from the A. kunkeei representative dataset
in_fasta=$basedir"/results/domains/blastp_species.faa" #Fasta file with both A. kunkeei domains and domains from other species
alignment=$basedir"/results/alignment/blastp_species.mafft.faa" #Output alignment
aln_log=$basedir"/logs/04.1-alignment_blastp.log" #Mafft log file
tree_log=$basedir"/logs/04.2-tree_blastp.log" #IQtree log file
outgroups=$basedir"/files/gtfB-like_outgroups.faa" #File with outgroups

mkdir -p $basedir/logs #Create log directory if it doesn't exist
mkdir -p $(dirname -- $alignment) #Create alignment directory if it doesn't exist

# Put outgroups and GH70s from A. kunkeei in the same file
cat $in_og $in_GH70 $outgroups > $in_fasta

# Run alignment
mafft-linsi --thread $threads $in_fasta > $alignment 2> $aln_log

# Create tree
iqtree --redo -nt $threads -s $alignment -st AA -msub nuclear -bb 1000 -bnni > $tree_log 2> $tree_log #Former model: -m LG+G4+F
