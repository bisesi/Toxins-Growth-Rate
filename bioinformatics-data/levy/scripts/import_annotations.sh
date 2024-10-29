#!/bin/zsh

# make  sure we're in the right working directory
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy

# make annotations file
mkdir /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/zip_files
mkdir /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/annotations

# activate conda environment
conda activate ncbi_datasets

# set zip files directory
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/zip_files

# download all zip files into annotations folder
while read line; do
    datasets download genome taxon "$line" --reference --include cds --filename "$line".zip
done < /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/taxon_names.txt

# get into folder, unzip, and rename
for f in *\ *; do mv "$f" "${f// /_}"; done

for file in *.zip; do
    unzip $file
    [ -f ncbi_dataset/data/GCF*/*.fna ] && mv ncbi_dataset/data/GCF*/*.fna /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/annotations/"$file".fna || echo "No file"
    rm -r ncbi_dataset
    rm -f README.md
    rm -f $file
done

# deactivate environment
conda deactivate

#clean final file names
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/annotations
for cdsfile in *; do mv "${cdsfile}" "${cdsfile/.zip/}"; done

# delete zip files folder
rm -f /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/levy/zip_files