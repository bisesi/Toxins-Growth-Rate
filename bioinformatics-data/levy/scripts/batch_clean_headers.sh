#!/bin/zsh

# Directory containing the files (update if needed)
dir="/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/rocha-test/img-rocha/"

# Loop through all .genes.faa files
for faa_file in "$dir"/*.genes.faa; do
    # Extract the numeric prefix (e.g., 637000004 from 637000004.genes.faa)
    prefix=$(basename "$faa_file" .genes.faa)

    # Define corresponding .genes.fna file
    fna_file="$dir/${prefix}.genes.fna"

    # Check if the corresponding .genes.fna file exists
    if [[ -f "$fna_file" ]]; then
        echo "Processing: $faa_file & $fna_file"

        temp_file="${fna_file}.tmp"

        # Extract headers from genes.faa and store them in a temporary file
        awk '/^>/ {key=$1" "$2; sub(/^>/, "", key); map[key]=$0} END {for (k in map) print k "|" map[k]}' "$faa_file" > faa_headers.tmp

        # Process genes.fna and replace headers if a match is found
        awk '
            BEGIN {
                while (getline < "faa_headers.tmp") {
                    split($0, parts, "|");
                    map[parts[1]] = parts[2];
                }
                close("faa_headers.tmp");
            }
            /^>/ {
                key=$1" "$2;
                sub(/^>/, "", key);
                if (key in map) {
                    print map[key];  # Replace with the corrected header
                } else {
                    print $0;  # Keep original if no match
                }
                next;
            }
            { print $0 }  # Print sequence lines unchanged
        ' "$fna_file" > "$temp_file"

        # Replace the original genes.fna with the modified one
        mv "$temp_file" "$fna_file"

        echo "Updated: $fna_file"
    else
        echo "Skipping: No matching .genes.fna for $faa_file"
    fi
done

# Cleanup temporary file
rm -f faa_headers.tmp