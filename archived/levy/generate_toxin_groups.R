# generate new toxin product group names for levy data

#load packages
library("tidyverse")
library("stringdist")
library("data.table")
library("stringi")

#load data
levy_ids <- read_csv(here::here("bioinformatics-data", "levy", "genome_ids.csv")) %>% 
  mutate(species_id = word(organism, 1,2, sep=" ")) %>%
  filter(!str_detect(species_id, "\\b\\w+\\b sp\\b")) %>%
  filter(!str_detect(species_id, "\\d")) %>%
  filter(str_detect(species_id, "^[A-Z][a-z]+\\s[a-z]+$")) %>%
  mutate(species_id = tolower(species_id)) %>%
  select(genome, species_id)
toxins_levy <- read_csv(here::here("bioinformatics-data", "levy", "all_genes.csv")) %>% inner_join(., levy_ids, by = "genome") %>% mutate(dataset = "levy")

write.csv(toxins_levy %>% select(genome, species_id, product_name) %>% mutate(product_name = tolower(product_name)), file = "genome_ids_species_ids_levy.csv")

aggregated_genomes <- toxins_levy %>% filter(protein_toxicity == "tox") %>% select(product_name, species_id) %>% mutate(product_name = tolower(product_name)) %>%
  filter(grepl("hypothetical", product_name) == FALSE & grepl("uncharacterized", product_name) == FALSE & grepl("plasmid stabili", product_name) == FALSE & 
           grepl("unknown", product_name) == FALSE & grepl("possible", product_name) == FALSE) %>% 
  filter(grepl("ase", product_name) == TRUE | grepl("addiction", product_name) == TRUE | grepl("acrd", product_name) == TRUE | grepl("cidin", product_name) == TRUE | grepl("ribo", product_name) == TRUE |
        grepl("multidrug", product_name) == TRUE | grepl("virulence", product_name) == TRUE | grepl("crispr", product_name) == TRUE | grepl("atp", product_name) == TRUE |
          grepl("abc", product_name) == TRUE | grepl("phage", product_name) == TRUE | grepl("efflux", product_name) == TRUE | grepl("toxin", product_name) == TRUE |
          grepl("death", product_name) == TRUE | grepl("lyso", product_name) == TRUE | grepl("cin", product_name) == TRUE | grepl("sin", product_name) == TRUE |
          grepl("damage", product_name) == TRUE | grepl("killer", product_name) == TRUE | grepl("repeat", product_name) == TRUE | grepl("abortive", product_name) == TRUE |
          grepl("cell wall", product_name) == TRUE | grepl("patho", product_name) == TRUE | grepl("glutinin", product_name) == TRUE | grepl("mate", product_name) == TRUE |
          grepl("regulat", product_name) == TRUE | grepl("transcrip", product_name) == TRUE) %>% unique()
 
#main groups 
aggregated_genomes <- aggregated_genomes %>% mutate(group = case_when(grepl("bacteriocin", product_name) == TRUE ~ "bacteriocin",
                                                grepl("colicin", product_name) == TRUE ~ "colicin",
                                                grepl("abc", product_name) == TRUE ~ "abc",
                                                grepl("atpase", product_name) == TRUE ~ "atpase",
                                                grepl("abortive", product_name) == TRUE ~ "abortive infection",
                                                grepl("addiction", product_name) == TRUE ~ "addiction module",
                                                grepl("adhesin", product_name) == TRUE ~ "adhesin",
                                                grepl("aerolysin", product_name) == TRUE ~ "aerolysin",
                                                grepl("atp-binding casette", product_name) == TRUE ~ "atp-binding casette",
                                                grepl("entericidin", product_name) == TRUE ~ "entericidin",
                                                grepl("endolysin", product_name) == TRUE ~ "endolysin",
                                                grepl("lysozyme", product_name) == TRUE ~ "lysozyme",
                                                grepl("cytolysin", product_name) == TRUE ~ "cytolysin",
                                                grepl("botulinum", product_name) == TRUE ~ "botulinum",
                                                grepl("multidrug", product_name) == TRUE | grepl("efflux", product_name) == TRUE ~ "multidrug",
                                                grepl("cell wall", product_name) == TRUE ~ "cell wall",
                                                grepl("enterotoxin", product_name) == TRUE ~ "enterotoxin",
                                                grepl("endotoxin", product_name, product_name) == TRUE ~ "endotoxin",
                                                grepl("exotoxin", product_name) == TRUE ~ "exotoxin",
                                                grepl("filamentous", product_name) == TRUE ~ "filamentous hemaglutinin",
                                                grepl("lysin", product_name) == TRUE ~ "lysin",
                                                grepl("leukocidin", product_name) == TRUE ~ "leukocidin",
                                                grepl("leukotoxin", product_name) == TRUE ~ "leukotoxin",
                                                grepl("pertussis", product_name) == TRUE ~ "pertussis",
                                                grepl("pyocin", product_name) == TRUE ~ "pyocin",
                                                grepl("ricin", product_name) == TRUE ~ "ricin",
                                                grepl("shiga", product_name) == TRUE ~ "shiga",
                                                grepl("ankyrin", product_name) == TRUE ~ "ankyrin",
                                                grepl("antitoxin", product_name) == TRUE ~ "toxin-antitoxin",
                                                grepl("nuclease", product_name) == TRUE ~ "nuclease",
                                                grepl("ribosom", product_name) == TRUE ~ "ribosomal",
                                                grepl("atp-binding", product_name) == TRUE ~ "atp-binding",
                                                grepl("insecticidal", product_name) == TRUE ~ "insecticidal",
                                                grepl("hydrolase", product_name) == TRUE ~ "hydrolase",
                                                grepl("cytotoxin", product_name) == TRUE ~ "cytotoxin",
                                                grepl("enterocin", product_name) == TRUE ~ "enterocin",
                                                grepl("protease", product_name) == TRUE ~ "protease",
                                                grepl("polymerase", product_name) == TRUE ~ "polymerase",
                                                grepl("ribonuclease", product_name) == TRUE ~ "ribonuclease",
                                                grepl("reductase", product_name) == TRUE ~ "reductase",
                                                grepl("dehydrogenase", product_name) == TRUE ~ "dehydrogenase",
                                                grepl("rnase", product_name) == TRUE ~ "rnase",
                                                grepl("kinase", product_name) == TRUE ~ "kinase",
                                                grepl("permease", product_name) == TRUE ~ "permease",
                                                grepl("interferase", product_name) == TRUE ~ "interferase",
                                                grepl("lipase", product_name) == TRUE ~ "lipase",
                                                grepl("recombinase", product_name) == TRUE ~ "recombinase",
                                                grepl("oxidase", product_name) == TRUE ~ "oxidase",
                                                grepl("galactosidase", product_name) == TRUE ~ "galactosidase",
                                                grepl("glucosidase", product_name) == TRUE ~ "glucosidase",
                                                grepl("arabinofuranosidase", product_name) == TRUE ~ "arabinofuranosidase",
                                                grepl("catalase", product_name) == TRUE ~ "catalase",
                                                grepl("cyclase", product_name) == TRUE ~ "cyclase",
                                                grepl("aromatase", product_name) == TRUE ~ "aromatase",
                                                grepl("mannanase", product_name) == TRUE ~ "mannanase",
                                                grepl("lyase", product_name) == TRUE ~ "lyase",
                                                grepl("peptidase", product_name) == TRUE ~ "peptidase",
                                                grepl("plasmid maintenance", product_name) == TRUE ~ "plasmid maintenance",
                                                grepl("synthetase", product_name) == TRUE ~ "synthetase",
                                                grepl("primase", product_name) == TRUE ~ "primase",
                                                grepl("helicase", product_name) == TRUE ~ "helicase",
                                                grepl("peptidase", product_name) == TRUE ~ "peptidase",
                                                grepl("integrase", product_name) == TRUE ~ "integrase",
                                                grepl("aminidase", product_name) == TRUE ~ "aminidase",
                                                grepl("xylanase", product_name) == TRUE ~ "xylanase",
                                                grepl("xylosidase", product_name) == TRUE ~ "xylosidase",
                                                grepl("glucanase", product_name) == TRUE ~ "glucanase",
                                                grepl("ligase", product_name) == TRUE ~ "ligase",
                                                grepl("rtx", product_name) == TRUE ~ "rtx",
                                                grepl("phosphodiesterase", product_name) == TRUE ~ "phosphodiesterase",
                                                grepl("fucosidase", product_name) == TRUE ~ "fucosidase",
                                                grepl("carboxylase", product_name) == TRUE ~ "carboxylase",
                                                grepl("pyrophosphatase", product_name) == TRUE ~ "pyrophosphatase",
                                                grepl("mannosidase", product_name) == TRUE ~ "mannosidase",
                                                grepl("amylase", product_name) == TRUE ~ "amylase",
                                                grepl("zeta", product_name) == TRUE ~ "zeta",
                                                grepl("chitinase", product_name) == TRUE ~ "chitinase",
                                                grepl("glutinin", product_name) == TRUE | grepl("gluttinin", product_name) == TRUE ~ "hemaglutinin",
                                                grepl("death", product_name) == TRUE ~ "death on curing",
                                                grepl("mutase", product_name) == TRUE ~ "mutase",
                                                grepl("deaminase", product_name) == TRUE ~ "deaminase",
                                                grepl("regulatory", product_name) == TRUE ~ "regulatory",
                                                TRUE ~ "none"))


#create product groups based on high jw similarity 
products_of_interest <- aggregated_genomes %>% filter(!is.na(product_name)) %>% filter(group == "none") %>%
  group_by(product_name) %>% summarize(n = n())

grouped_products <- map_dfr(products_of_interest$product_name, ~ {
  i <- which(stringdist(., products_of_interest$product_name, "jw") < 0.2)
  tibble(index = i, title = products_of_interest$product_name[i])
}, .id = "group") %>%
  distinct(index, .keep_all = T) %>% 
  mutate(group = as.integer(group))

numbered_groups <- aggregated_genomes %>% filter(group == "none") %>% inner_join(., grouped_products %>% rename(product_name = title), by = "product_name") %>%
  select(-c(group.x, index)) %>% rename(group = group.y)

#combine with orginal levy data 
with_group_number <- aggregated_genomes %>% filter(group != "none") %>%
  rbind(., numbered_groups) %>%
  select(product_name, species_id, group) %>% unique()

#write group names csv
write.csv(with_group_number, file = here::here("bioinformatics-data", "levy", "levy_toxin_data_with_groups.csv"))

# final dataset for evaluation
levy <- read_csv(here::here("bioinformatics-data", "levy", "full_growth_toxin_dataset_levy.csv")) %>% dplyr::select(-c(`...1`)) %>% 
  inner_join(., read_csv(here::here("bioinformatics-data", "levy", "genome_ids_species_ids_levy.csv")), by = c("species_id", "product_name")) %>% 
  dplyr::select(-c(`...1`)) %>% unique() %>%
  filter(species_id != "chloracidobacterium thermophilum") %>% select(-c(phon)) %>% 
  inner_join(., read_csv(here::here("bioinformatics-data", "levy", "levy_toxin_data_with_groups.csv")) %>% dplyr::select(-c(`...1`)), by = c("species_id", "product_name")) %>%
  mutate(group = case_when(grepl("yhav", product_name) == TRUE ~ "yhav",
                           grepl("yoeb", product_name) == TRUE ~ "yoeb",
                           grepl("chpk", product_name) == TRUE ~ "chpk",
                           grepl("ccdb", product_name) == TRUE ~ "ccdb",
                           grepl("cpta", product_name) == TRUE ~ "cpta",
                           grepl("yafq", product_name) == TRUE ~ "yafq",
                           grepl("pare", product_name) == TRUE ~ "pare",
                           grepl("ybaj", product_name) == TRUE ~ "ybaj",
                           grepl("ypjf", product_name) == TRUE ~ "ypjf",
                           grepl("yafo", product_name) == TRUE ~ "yafo",
                           grepl("exfoliative", product_name) == TRUE ~ "exfoliative",
                           grepl("esterase", product_name) == TRUE ~ "esterase",
                           grepl("pemk", product_name) == TRUE ~ "pemk",
                           grepl("vapc", product_name) == TRUE ~ "vapc",
                           TRUE ~ group)) %>% full_join(., read_csv(here::here("bioinformatics-data", "levy", "named_toxin_groups.csv")) %>% 
                                                          select(-c(`...1`)) %>% rename(group = group_number) %>% mutate(group = as.character(group)), by = "group") %>%
  mutate(group = ifelse(!is.na(group_name), group_name, group)) %>% select(-c(group_name)) %>% mutate(nums = as.numeric(group)) %>% filter(is.na(nums)) %>% select(-c(nums)) %>%
  filter(!is.na(species_id))

write.csv(levy, file = here::here("bioinformatics-data", "levy", "complete_levy_dataset.csv"))

