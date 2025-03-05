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
  return(cleaned_string)
}

#create dataset
files <- list.files(here::here("bioinformatics-data", "rocha", "annotations"), pattern = "\\.fna$")

annotation_list = list()
new_names = c()
for (i in files){
  index <- which(files == i)
  path <- here::here("bioinformatics-data", "rocha", "annotations", i)
  genes <- readDNAStringSet(path)
  ribo <- grepl("16S",names(genes),ignore.case = T)
  if (length(which(ribo == TRUE)) < 1){
    next
  }
  concat <- DNAStringSet(paste(as.character(genes[ribo]), collapse = ""))
  annotation_list[[index]] <- concat
  new_names[index] <- convert_filename(i)
}

annotation_list <- annotation_list[-c(50)]
new_names <- new_names[-c(50)]
full_16S <- DNAStringSet(unlist(lapply(annotation_list, as.character)))
names(full_16S) <- basename(tolower(new_names))

align <- msa(full_16S, method = "ClustalW")
tweak_align <- msaConvert(align, type = "seqinr::alignment")
dist <- seqinr::dist.alignment(tweak_align, matrix = "identity")
tree_subset <- njs(dist)

png(here::here("figures", "final-figs", "imgs", "rocha-phylo.png"), res = 300, width = 4000, height = 16000)
plot.phylo(tree_subset, main="Phylogenetic Tree", 
           use.edge.length = T)
dev.off()




