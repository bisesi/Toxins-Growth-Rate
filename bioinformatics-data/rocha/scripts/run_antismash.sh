#!/bin/zsh

# set path
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha/annotations

# create antimash_tables folder
mkdir /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha/antismash_tables

# activate conda
conda activate antismash

# run antismash on everything
for file in *.gbff; do
	antismash $file
done

#close conda
conda deactivate

# go into folders, rename summary html, and move into annotations folder
for f in */*.html ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp/$fp"_"$fl"; done
mv */*index.html /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha/antismash_tables

# remove all extraneous folders
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha/annotations 
rm -R -- *(/)
