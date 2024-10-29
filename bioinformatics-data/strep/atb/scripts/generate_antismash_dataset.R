#generate antismash bcg data for strep data set

#load packages
library("rvest")
library("tidyverse")

#use a helper function to clean filenames
convert_filename <- function(filename) {
  base_name <- tools::file_path_sans_ext(filename)
  result <- gsub("_", " ", base_name)
  cleaned_string <- gsub("\\[|\\]", "", result)
  cleaned_string <- str_replace_all(cleaned_string, "\\bindex\\b", "")
  cleaned_string <- str_squish(cleaned_string)
  return(cleaned_string)
}

#use another helper function to grab the organism name from each genbank file
readGBFF <- function (file, text = readLines(file), verbose = FALSE){
  line_index <- grep("ORGANISM", text)
  organism_line <- text[line_index[1]]
  cleaned_string <- trimws(organism_line)
  cleaned_string <- sub("^\\s*ORGANISM\\s+", "", cleaned_string)
  return(cleaned_string)
}

# load all html tables
files <- list.files(here::here("bioinformatics-data", "strep", "antismash_tables"), pattern = "\\.html$")

predicted_secondary_metabolites <- data.frame()
no_bcgs <- vector()
for (i in files){
  temp <- read_html(here::here("bioinformatics-data", "strep", "antismash_tables", i)) %>% 
    html_table()
  if (length(temp) == 0) {
    no_bcgs <- append(no_bcgs, tolower(convert_filename(i)))
  } else {
    data <- last(temp)
    data$species_id <- tolower(convert_filename(i))
    predicted_secondary_metabolites <- rbind(predicted_secondary_metabolites, data %>% janitor::clean_names()) 
  }
  print(paste(which(files == i), "of", length(files)))
}

#grab organism name
files <- list.files(here::here("bioinformatics-data", "strep", "strep_files_ncbi"), pattern = "\\.gbff$")
species_names <- data.frame()
for (i in files){
  name <- readGBFF(here::here("bioinformatics-data", "strep", "strep_files_ncbi", i))
  numerical_string <- tolower(convert_filename(i))
  species_names <- rbind(species_names, data.frame(species = name, species_id = numerical_string))
  print(paste(which(files == i), "of", length(files)))
}

# combine dataset of species with and without bcgs
full_dataset <- predicted_secondary_metabolites %>% dplyr::select(c(region, type, from, to, species_id)) %>%
  inner_join(., species_names, by = "species_id", relationship = "many-to-many")

# export final data set
write.csv(full_dataset, file = here::here("bioinformatics-data", "strep", "antismash_data.csv"))


