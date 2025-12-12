for strain in "IBH001" "DSMZ12361";
do
interproscan.sh -i ../Akunkeei_files/faa/${strain}_protein.faa -f tsv -b ${strain}_id -dp -cpu 64;
done 
