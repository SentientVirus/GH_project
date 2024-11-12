threads=16 #No. of threads to be used
in_og="files/ingroups.faa" #Input fasta file with complete protein sequences
out_og="results/interproscan/GH70_species.tsv" #Output file
annot_log="logs/01-interproscan.log" #Log file

mkdir -p $(dirname -- $annot_log) #Create the log directory if it doesn't exist
mkdir -p $(dirname -- $out_og) #Create the output directory if it doesn't exist

# Get domain annotations of outgroups
~/interproscan-5.59-91.0/interproscan.sh -i $in_og -f tsv -o $out_og -dp -cpu $threads 2> $annot_log;
