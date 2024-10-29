#!/bin/zsh

# make sure we are in the working directory that we care about
cd /Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/bioinformatics-data/strep

# activate ncbi conda environment
conda activate ncbi_datasets

# import all streptomyces datasets using dehydrate with cds, gene and gbff formats
datasets download genome taxon streptomyces --assembly-level complete --include cds,gbff --dehydrated --filename complete_strep_dataset.zip

unzip complete_strep_dataset.zip -d strep_files_ncbi

datasets rehydrate --directory strep_files_ncbi

conda deactivate

# remove the zip file

rm -f complete_strep_dataset.zip


