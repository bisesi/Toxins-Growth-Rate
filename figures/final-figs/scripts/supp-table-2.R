# ATB
# supp table 2
# list of marine taxa

library("tidyverse")
library("cowplot")
library("tidytext")
library("ggtext")
library("boot")
library("mosaic")

# load data
levy <- read_csv(here::here("bioinformatics-data", "levy", "full_growth_toxin_dataset_levy.csv")) %>%
  dplyr::select(-c(`...1`))

# marine snow or marine sediment associated
alphaproteo <- c("roseobacter", "temperatibacter", "yoonia",
                 "amphiplicatus", "aquisalinus", "hyphococcus", "marinibacterium",
                 "marinicaulis", "parvularcula", "rhodobacter", "roseovarius",
                 "oceanibacterium", "snaethiella", "sphingomonas", "sphingobium", "sphingorhabdus",
                 "sphingopyxis", "blastomonas", "lutibacterium", "sandarakinorhabdus",
                 "sandaracinobacter", "citromicrobium", "novosphingobium", "erythrobacter",
                 "rhodocista", "azospirillum", "rhodospirillum", "rhodospira", "pararhodospirillum",
                 "algihabitans", "defluviicoccus", "roseospira", "dongia", "caenispirillum",
                 "reyranella", "magnetospirillum", "roseospirillum", "inquilinus", "rhodovibrio",
                 "sulfitobacter", "planktomarina", "marinibacterium", "halovulum", "tropicimonas",
                 "rhodothalassium") # 194 species
flavobacteriales <- c("flavobacterium", "kriegella", "gelidibacter", "leeuwenhoekiella", "aequorivita",
                      "olleya", "salinimicrobium", "xanthomarina", "bergeyella", "gramella", 
                      "psychoroserpens", "sabulilitoribacter", "lutibacter", "aurantivirga", "nonlabens",
                      "mesoflavibacter", "imtechella", "psychroflexus", "joostella", "zunongwangia",
                      "lacinutrix", "formosa", "facecalibacter", "sinomicrobium") # 124
oceanospirallales <- c("oceanospirillum", "oleispira", "oceaniserpentilla", "halomonas", "halotalea", "carnimonas",
                       "marinospirillum", "nitrincola", "profundimonas", "alcanivorax", "oceanobacter", "balneatrix",
                       "kushneria", "bermanella", "neptuniibacter", "hahella", "thalassolituus", "litoricola",
                       "marinomonas", "zymobacter", "halovibrio", "neptunomonas", "amphritea") # 68 species
campylobacterales <- c("arcobacter", "nitratiruptor", "helicobacter", 
                       "nitratifractor", "campylobacter", "hydrogenimonas") # 58
actinobacteridae <- c("actinotalea", "brachybacterium", "brevibacterium",
                      "dietzia", "leucobacter", "microbacterium", "micrococcus",
                      "rhodococcus") # 82
vibrionales <- c("vibrio", "salinivibrio", "photobacterium",
                 "allivibrio", "enterovibrio", "thaumasiovibrio",
                 "listonella", "photococcus", "catenococcus", 
                 "echinomonas", "allomonas", "beneckea", "enhydrobacter") # 87 species
alteromonadales <- c("shewanella", "motilimonas", "pseudoalteromonas",
                     "alteromonas", "aestuariicella", "alishewanella",
                     "agarivorans", "haliea", "salinimonas", "alkalimarinus",
                     "lacimicrobium", "paraglaciecola", "marinobacter",
                     "microbulbifer", "catenovulum", "planctobacterium",
                     "glaciecola", "aestuariibacter", "marinobacterium",
                     "aliiglaciecola", "alginatibacterium", "colwellia") # 106 species

outliers <- c(2772190761)

all_marine <- unique(c(alphaproteo, flavobacteriales, oceanospirallales,
                       campylobacterales, actinobacteridae, vibrionales, alteromonadales))

subset <- levy %>% filter(str_extract(species_id, "^[^ ]+") %in% all_marine)

# export csv file
supptable2 <- subset %>% select(species_id) %>% unique() 

write.csv(supptable2, file = here::here("figures", "final-figs", "tables", "supplemental-table-2.csv"))

