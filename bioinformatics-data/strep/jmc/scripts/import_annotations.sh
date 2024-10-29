#!/bin/zsh

# make  sure we're in the right working directory
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/jmc

# activate conda environment
conda activate ncbi_datasets

# loop through each line of the file taxon_names.txt
while read line; do
    datasets download genome accession "$line" --include gbff,cds --filename "$line".zip
    unzip "$line".zip
    cd ncbi_dataset/data/GCF*
    mv *.fna "$line".fna
    mv *.gbff "$line".gbff
    mv * /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/jmc
    cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/jmc
    rm -r ncbi_dataset
    rm -f README.md
    rm -f "$line".zip	
    rm md5sum.txt
done < jmc_accessions.txt

# deactivate environment
conda deactivate

# make and move all files
mkdir annotations
mv *.fna /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/jmc/annotations
mv *.gbff /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep/jmc/annotations
