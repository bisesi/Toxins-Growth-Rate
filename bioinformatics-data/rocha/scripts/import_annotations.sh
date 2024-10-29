#!/bin/zsh

# make  sure we're in the right working directory
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha

# activate conda environment
conda activate ncbi_datasets

# loop through each line of the file taxon_names.txt
while read line; do
    datasets download genome taxon "$line" --reference --include gbff,cds --filename "$line".zip
    unzip "$line".zip
    cd ncbi_dataset/data/GCF*
    mv *.fna "$line".fna
    mv *.gbff "$line".gbff
    for f in *\ *; do mv "$f" "${f// /_}"; done
    mv * /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha
    cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha
    rm -r ncbi_dataset
    rm -f README.md
    rm -f "$line".zip	
done < taxon_names.txt

# deactivate environment
conda deactivate

# make and move all files
mkdir annotations
mv *.fna /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha/annotations
mv *.gbff /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/rocha/annotations
