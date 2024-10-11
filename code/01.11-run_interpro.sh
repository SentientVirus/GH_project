threads=$4
in_og=$1
out_og=$2
annot_log=$3

mkdir -p $(dirname -- $annot_log)
mkdir -p $(dirname -- $out_og)

# Get domain annotations of outgroups
~/interproscan-5.59-91.0/interproscan.sh -i $in_og -f tsv -o $out_og -dp -cpu $threads 2> $annot_log > $annot_log;
