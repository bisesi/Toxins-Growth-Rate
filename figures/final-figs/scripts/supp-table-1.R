# supplemental table generation
# ATB

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")
library("cowplot")

# import datasets
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  pull(species_id) %>% unique()

strep <- data.frame(read.delim(here::here("bioinformatics-data", "strep", "all_strep_accession_numbers.txt"))) %>%
  pull(accessions)

levy <- read_csv(here::here("bioinformatics-data", "levy", "full_growth_toxin_dataset_levy.csv")) %>%
  pull(species_id) %>% unique()

supptable1 = rbind(data.frame(dataset = c("Viera-Silva and Rocha"), ncbi_taxon_name = rocha),
                   data.frame(dataset = c("Streptomyces"), ncbi_taxon_name = strep),
                   data.frame(dataset = c("Toxinome"), ncbi_taxon_name = levy))

write.csv(supptable1, file = here::here("figures", "final-figs", "tables", "supplemental-table-1.csv"))
