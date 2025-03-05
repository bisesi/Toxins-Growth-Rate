# ATB
# get accession numbers from levy fna files
# 5 Feb 2025

#load packages
library("gRodon")
library("Biostrings")
library("tidyverse")

#create dataset
files <- list.files(here::here("bioinformatics-data", "levy", "annotations"), pattern = "\\.fna$")

accessions <- c()
for (i in files){
  path <- here::here("bioinformatics-data", "levy", "annotations", i)
  first_line <- readLines(path, n = 1)
  accession <- str_extract(first_line, "(?<=\\|)[^_]+_[^_]+")
  print(accession)
  accessions <- c(accessions, accession)
}

writeLines(accessions, "accession_numbers.txt")
