# fix levy dataset for ncbi

dat <- read_csv("genome_ids.csv")

new_names <- dat %>% mutate(new_name = word(organism, 1,2, sep=" ")) %>%
  filter(!str_detect(new_name, "\\b\\w+\\b sp\\b")) %>%
  filter(!str_detect(new_name, "\\d")) %>%
  filter(str_detect(new_name, "^[A-Z][a-z]+\\s[a-z]+$"))

unique_names <- unique(new_names$new_name)

writeLines(tolower(unique_names), "taxon_names.txt")
