#generate gRodon data for levy data set

#load packages
library("gRodon")
library("Biostrings")
library("tidyverse")

#use a helper function to clean filenames
convert_filename <- function(filename) {
  base_name <- tools::file_path_sans_ext(filename)
  result <- gsub(".genes", "", base_name)
  return(result)
}

#create dataset
files <- list.files(here::here("bioinformatics-data", "levy", "rocha-test", "img-rocha"), pattern = "\\.fna$")

predicted_growth_rate_data <- data.frame()
for (i in files){
  skip_to_next <- FALSE
  tryCatch({
    path <- here::here("bioinformatics-data", "levy", "rocha-test", "img-rocha", i)
    genes <- readDNAStringSet(path)
    highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
    output <- predictGrowth(genes, highly_expressed)
    temp <- data.frame(genome = tolower(convert_filename(i)), predicted_d = output$d, upperci = output$UpperCI, lowerci = output$LowerCI)
    predicted_growth_rate_data <- rbind(predicted_growth_rate_data, temp)
    print(paste(which(files == i), "of", length(files)))
  }, error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next){next}
}

#export dataset
write.csv(predicted_growth_rate_data, file = here::here("bioinformatics-data", "levy", "rocha-test", "img_grodon_growthrates.csv"))


