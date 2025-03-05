#load packages
library("seqinr")
library("Biostrings")
library("tidyverse")
library("ape")
library("msa")

# convert filenames
convert_filename <- function(filename) {
  base_name <- tools::file_path_sans_ext(filename)
  result <- gsub("_", " ", base_name)
  cleaned_string <- gsub("\\[|\\]", "", result)
  cleaned_string <- gsub("cds from genomic", "", cleaned_string)
  cleaned_string <- str_squish(cleaned_string)
  return(cleaned_string)
}

convert_strepname <- function(name) {
  result <- gsub("genomic", "", name)
  return(str_squish(result))
}

# get used IDs
strep <- read_csv(here::here("bioinformatics-data", "strep", "full_growth_toxin_dataset_strep.csv")) %>%
  dplyr::select(-c(`...1`)) %>% pull(species_id) %>% unique()

strep_names <- c()
for (i in strep){
  index <- which(strep == i)
  strep_names[index] <- convert_strepname(i)
}

#create dataset
files <- list.files(here::here("bioinformatics-data", "strep", "annotations"), pattern = "\\.fna$")

correct_list <- c()
for (i in files){
  if (tolower(convert_filename(i)) %in% strep_names){
    len <- length(correct_list)
    correct_list[len + 1] <- i
  }
}

annotation_list = list()
new_names = c()
no_16S = c()
for (i in correct_list){
  index <- which(correct_list == i)
  path <- here::here("bioinformatics-data", "strep", "annotations", i)
  genes <- readDNAStringSet(path)
  ribo <- grepl("16S",names(genes),ignore.case = T)
  if (length(which(ribo == TRUE)) < 1){
    len <- length(no_16S)
    no_16S[len + 1] <- index
    next
  }
  concat <- DNAStringSet(paste(as.character(genes[ribo]), collapse = ""))
  annotation_list[[index]] <- concat
  new_names[index] <- convert_filename(i)
}

annotation_list <- annotation_list[-c(no_16S)]
new_names <- new_names[-c(no_16S)]
full_16S <- DNAStringSet(unlist(lapply(annotation_list, as.character)))
names(full_16S) <- basename(tolower(new_names))

align <- msa(full_16S, method = "ClustalW")
tweak_align <- msaConvert(align, type = "seqinr::alignment")
dist <- seqinr::dist.alignment(tweak_align, matrix = "identity")
tree_subset <- njs(dist)

png(here::here("figures", "final-figs", "imgs", "strep-phylo.png"), res = 300, width = 4000, height = 16000)
plot.phylo(tree_subset, main="Phylogenetic Tree", 
           use.edge.length = T)
dev.off()



