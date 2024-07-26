threads=16
in_og="files/ingroups.faa"
out_og="results/interproscan/GH70_species.tsv"
annot_log="logs/01-interproscan.log"

mkdir -p $(dirname -- $annot_log)
mkdir -p $(dirname -- $out_og)

# Get domain annotations of outgroups
~/interproscan-5.59-91.0/interproscan.sh -i $in_og -f tsv -o $out_og -dp -cpu $threads 2> $annot_log;
