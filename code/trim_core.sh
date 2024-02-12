#!/bin/bash

input_dir="untrimmed_core"
output_dir="trimmed_core"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through files ending with genes.fasta in the input directory
for file in "$input_dir"/*genes.fasta; do
    filename=$(basename "$file")
    output_file="$output_dir/$filename"
    
    # Trim the file using Trimal with desired options
    trimal -in "$file" -out "$output_file" -gappyout

    echo "Trimmed file: ${output_file}"
done
