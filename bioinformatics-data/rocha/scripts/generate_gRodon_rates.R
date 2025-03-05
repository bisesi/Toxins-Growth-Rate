# generate gRodon data for rocha data set

#load packages
library("gRodon")
library("Biostrings")
library("tidyverse")

#use a helper function to clean filenames
convert_filename <- function(filename) {
  base_name <- tools::file_path_sans_ext(filename)
  result <- gsub("_", " ", base_name)
  cleaned_string <- gsub("\\[|\\]", "", result)
  return(cleaned_string)
}

#create dataset
files <- list.files(here::here("bioinformatics-data", "rocha"), pattern = "\\.fna$")

predicted_growth_rate_data <- data.frame()
for (i in files){
  path <- here::here("bioinformatics-data", "rocha", "annotations", i)
  genes <- readDNAStringSet(path)
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
  output <- predictGrowth(genes, highly_expressed)
  growthrate <- output$d
  upperci <- output$UpperCI
  lowerci <- output$LowerCI
  temp <- data.frame(species = tolower(convert_filename(i)), predicted_d = growthrate, upperci = upperci, lowerci = lowerci)
  predicted_growth_rate_data <- rbind(predicted_growth_rate_data, temp)
  print(i)
}

#export dataset
write.csv(predicted_growth_rate_data, file = here::here("bioinformatics-data", "rocha", "grodon_growthrate_data.csv"))


