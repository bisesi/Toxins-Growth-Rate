# generate levy dataset

# load levy data sets
growth_rates_levy <- read_csv(here::here("bioinformatics-data", "levy", "grodon_growthrate_data.csv")) %>% 
  dplyr::select(-c(`...1`)) %>% unique()
toxins_levy <- read_csv(here::here("bioinformatics-data", "levy", "genome_ids_species_ids_levy.csv")) %>%
  dplyr::select(-c(`...1`)) %>% unique()
levy_data <- growth_rates_levy %>% inner_join(., toxins_levy, by = "species_id") %>%
  inner_join(., read_csv(here::here("bioinformatics-data", "levy", "group_names.csv")) %>% mutate(product_name = tolower(product_name)), by = "product_name")

#write levy csv
write.csv(levy_data, file = here::here("bioinformatics-data", "levy", "full_growth_toxin_dataset_levy.csv"))


