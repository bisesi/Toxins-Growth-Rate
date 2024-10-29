#!/bin/zsh

# get into folder
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/strep_files_ncbi/ncbi_dataset/data

# rename files according to folder name
for f in */*.fna ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp/$fp"_"$fl"; done #cds files
for f in */*.gbff ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp/$fp"_"$fl"; done #gbff files

# move files to correct directory
mv ./*/*.fna /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/strep_files_ncbi
mv ./*/*.gbff /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/strep_files_ncbi

# remove ncbi_dataset file
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/strep_files_ncbi
rm -r ncbi_dataset
rm -f README.md
