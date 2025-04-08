#Environment: alignment_tree.yml
threads=16 #No. of threads to be used
workdir="/home/marina/GH_project/add_species"
in_og=$workdir"/files/blastp_formatted.faa" #Input fasta file with complete protein sequences
out_og=$workdir"/results/interproscan/blastp_species.tsv" #Output file
annot_log=$workdir"/logs/02-interproscan_blastp.log" #Log file

mkdir -p $(dirname -- $annot_log) #Create the log directory if it doesn't exist
mkdir -p $(dirname -- $out_og) #Create the output directory if it doesn't exist

# Get domain annotations of outgroups
~/interproscan-5.59-91.0/interproscan.sh -i $in_og -f tsv -o $out_og -dp -cpu $threads 2> $annot_log;
