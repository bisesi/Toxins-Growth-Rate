# create integrated datasets

#load packages
library("tidyverse")

# rocha
growth_rates_rocha <- read_csv(here::here("bioinformatics-data", "rocha", "grodon_growthrate_data.csv")) %>%
  dplyr::select(-c(`...1`)) %>% unique()
toxins_rocha <- read_csv(here::here("bioinformatics-data", "rocha", "antismash_data.csv")) %>%
  dplyr::select(-c(`...1`)) %>%
  separate_longer_delim(type, delim = ",") %>% unique()
rocha_data <- growth_rates_rocha %>% inner_join(., toxins_rocha, by = "species") %>% mutate(dataset = "rocha") %>% rename(species_id = species)
rocha_experimental <- read_csv(here::here("bioinformatics-data", "rocha", "supp-table-1.csv")) %>% janitor::clean_names()
all_rocha <- rocha_experimental %>% select(d_h, ncbi_name, altered_species_name) %>% rename(species_id = ncbi_name) %>% 
  inner_join(., rocha_data, by = "species_id", relationship = "many-to-many") %>% filter(is.na(altered_species_name)) %>% unique()

#write rocha csv
write.csv(all_rocha, file = here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv"))

# strep
#atb
growth_rates_strep_atb <- read_csv(here::here("bioinformatics-data", "strep", "atb", "grodon_growthrate_data.csv")) %>% 
  dplyr::select(-c(`...1`)) %>% unique()
toxins_strep_atb <- read_csv(here::here("bioinformatics-data", "strep", "atb", "antismash_data.csv")) %>% 
  dplyr::select(-c(`...1`)) %>%
  separate_longer_delim(type, delim = ",") %>% unique()
species_ids <- growth_rates_strep_atb %>% inner_join(., toxins_strep_atb, by = "species_id") %>% mutate(dataset = "strep") %>%
  mutate(transformed_strings = gsub("[^0-9\\.]", "", species_id)) %>%
  mutate(type = ifelse(grepl("gcf", species_id) == TRUE, "GCF", "GCA")) %>%
  arrange(transformed_strings, desc(type == "GCF")) %>%
  distinct(transformed_strings, .keep_all = TRUE) %>% pull(species_id) %>% unique()
strep_data_atb <- growth_rates_strep_atb %>% inner_join(., toxins_strep_atb, by = "species_id") %>% mutate(dataset = "strep") %>%
  filter(species_id %in% species_ids)

#jmc
growth_rates_strep_jmc <- read_csv(here::here("bioinformatics-data", "strep", "jmc", "grodon_growthrate_data.csv")) %>% 
  dplyr::select(-c(`...1`)) %>% unique()
toxins_strep_jmc <- read_csv(here::here("bioinformatics-data", "strep", "jmc", "antismash_data.csv")) %>% 
  dplyr::select(-c(`...1`)) %>%
  separate_longer_delim(type, delim = ",") %>% unique()
strep_data_jmc <- growth_rates_strep_jmc %>% inner_join(., toxins_strep_jmc, by = "species_id") %>% mutate(dataset = "strep")

#write strep csv
write.csv(rbind(strep_data_atb %>% select(-c(species)), strep_data_jmc), file = here::here("bioinformatics-data", "strep", "full_growth_toxin_dataset_strep.csv"))


#levy 
species <- read_csv(here::here("bioinformatics-data", "levy", "genome_ids.csv")) %>%
  dplyr::select(-c(`...1`))
genes <- read_csv(here::here("bioinformatics-data", "levy", "all_genes.csv"))
pfams <- read_csv(here::here("bioinformatics-data", "levy", "all_pfams.csv"))
genes_pfams <- read_csv(here::here("bioinformatics-data", "levy", "all_genes_pfams.csv"))
unique_proteinID_pfamID_pairs <- genes_pfams %>% select(proteinID, pfamID) %>% unique()
levy_ids <- species %>% 
  mutate(species_id = word(organism, 1,2, sep=" ")) %>%
  filter(!str_detect(species_id, "\\b\\w+\\b sp\\b")) %>%
  filter(!str_detect(species_id, "\\d")) %>%
  filter(str_detect(species_id, "^[A-Z][a-z]+\\s[a-z]+$")) %>%
  mutate(species_id = tolower(species_id)) %>%
  select(genome, species_id)

unique_ids_per_genome <- genes %>% filter(protein_toxicity == "tox") %>%
  select(genome, proteinID, product_name) %>%
  inner_join(., unique_proteinID_pfamID_pairs, by = "proteinID") %>%
  inner_join(., levy_ids, by = "genome")

growth_rates_levy <- read_csv(here::here("bioinformatics-data", "levy", "grodon_growthrate_data.csv")) %>% 
  dplyr::select(-c(`...1`)) %>% unique()

#concatenate levy data and dump all plasmid maintenance toxins
complete_levy <- unique_ids_per_genome %>% inner_join(., growth_rates_levy, by = "species_id") %>%
  select(genome, product_name, species_id, predicted_d, pfamID, proteinID) %>%
  mutate(product_name = ifelse(is.na(product_name), "unknown", product_name)) %>%
  mutate(product_name = tolower(product_name)) %>% filter(grepl("plasmid maintenance", product_name) == FALSE & grepl("plasmid", product_name) == FALSE &
                                                            grepl("toxin-antitoxin", product_name) == FALSE & grepl("addiction", product_name) == FALSE) 

#write levy csv
write.csv(complete_levy, file = here::here("bioinformatics-data", "levy", "full_growth_toxin_dataset_levy.csv"))




