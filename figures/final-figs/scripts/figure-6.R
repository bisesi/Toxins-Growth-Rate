# ATB
# figure 6
# propagating error

# load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")
library("reshape2")
library("scales")
library("deldir")
library("ggforce")

# functions
generate_uniform_samples <- function(df, bootstrap_number) {
  # Ensure the data frame has the required columns
  if (!all(c("species_id", "lowerci", "upperci") %in% names(df))) {
    stop("The data frame must contain 'species_id', 'lowerci', and 'upperci' columns.")
  }
  
  # Generate the uniform samples for each row
  result <- lapply(1:nrow(df), function(i) {
    species_id <- df$species_id[i]
    lower <- df$lowerci[i]
    upper <- df$upperci[i]
    
    if (lower >= upper) {
      stop(paste("Invalid range for species_id:", species_id))
    }
    
    samples <- runif(bootstrap_number, min = lower, max = upper)
    data.frame(species_id = species_id, sample = samples)
  })
  
  # Combine all results into a single data frame
  result_df <- do.call(rbind, result)
  return(result_df)
}

generate_normal_samples <- function(df, bootstrap_number, sd_type, sd_provided = NA) {
  # Ensure the data frame has the required columns
  if (!all(c("species_id", "lowerci", "upperci", "predicted_d") %in% names(df))) {
    stop("The data frame must contain 'species_id', 'lowerci', 'predicted_d', and 'upperci' columns.")
  }
  
  # Generate the uniform samples for each row
  result <- lapply(1:nrow(df), function(i) {
    species_id <- df$species_id[i]
    lower <- df$lowerci[i]
    upper <- df$upperci[i]
    mean <- df$predicted_d[i]
    
    if (lower >= upper) {
      stop(paste("Invalid range for species_id:", species_id))
    }
    
    if (sd_type == "from_interval") {
      samples <- rnorm(bootstrap_number, mean = mean, sd = (upper - lower)/(2 * 1.96))
    } else if (sd_type == "preset") {
      samples <- rnorm(bootstrap_number, mean = mean, sd = sd_provided)
    }
    
    data.frame(species_id = species_id, sample = samples)
  })
  
  # Combine all results into a single data frame
  result_df <- do.call(rbind, result)
  return(result_df)
}

# load data and set seed
strep <- read_csv(here::here("bioinformatics-data", "strep", "full_growth_toxin_dataset_strep.csv")) %>%
  dplyr::select(-c(`...1`))
set.seed(123)

# get bcg presence and absence
all_bcgs <- strep %>%
  filter(!is.na(type)) %>% pull(type) %>% unique()
bcg_presence_absence <- data.frame()
for (i in all_bcgs){
  if (i == "terpene" | i == "siderophore") {
    all <- data.frame(species_id = strep %>% pull(species_id) %>% unique(), present = 1, type = i)
    bcg_presence_absence <- rbind(bcg_presence_absence, all)
    next
  }
  true_set <- strep %>% filter(type == i) %>% pull(species_id) %>% unique()
  false_set <- setdiff(strep %>% pull(species_id) %>% unique(), true_set)
  bool_0 <- data.frame(species_id = false_set, present = as.integer(0), type = i)
  bool_1 <- data.frame(species_id = true_set, present = as.integer(1), type = i)
  bcg_presence_absence <- rbind(bcg_presence_absence, bool_0, bool_1)
}

strep_bcgs <- strep %>% dplyr::select(species_id, predicted_d) %>% unique() %>%
  inner_join(., bcg_presence_absence, by = "species_id") %>% 
  filter(present == 1) %>%
  pull(type) %>% unique()

bootstrap_number <- 1000

# significant models in original data
sig <- c("hserlactone", "aminocoumarin", "RRE-containing",
         "prodigiosin", "arylpolyene", "indole", "betalactone",
         "lanthipeptide-class-ii", "lanthipeptide-class-iii",
         "furan", "thiopeptide", "phosphonate", "other", "NRPS-like",
         "blactam", "LAP", "T2PKS", "NAPAA", "transAT-PKS-like",
         "amglyccycl", "lanthipeptide-class-iv", "cyanobactin", "T1PKS",
         "butyrolactone", "NRPS", "melanin", "ranthipeptide", "phosphoglycolipid", "RiPP-like")

# do all sampling
#uniform
uniform_samples <- generate_uniform_samples(df = strep %>% dplyr::select(species_id, upperci, lowerci) %>% unique(), 
                                            bootstrap_number = bootstrap_number)

sample_data <- uniform_samples %>% group_by(species_id) %>% mutate(replicate_number = 1:n()) %>%
  ungroup() %>% mutate(rate = log(2) / sample)

all_models <- data.frame()
for (i in strep_bcgs){
  presence <- bcg_presence_absence %>% filter(type == i)
  full_data <- sample_data %>% inner_join(., presence, by = "species_id")
  models <- full_data %>% dplyr::select(species_id, present, type, rate, replicate_number) %>%
    group_by(replicate_number) %>%
    nest() %>% 
    mutate(log_model = map(data, ~glm(cbind(present, 1 - present) ~ rate, data = ., family = "binomial"))) %>% 
    mutate(tidy = map(log_model, broom::tidy)) %>% 
    unnest(tidy) %>%
    dplyr::select(replicate_number, term, estimate, p.value) %>%
    filter(term != "(Intercept)") %>%
    ungroup() %>% mutate(type = i)
  all_models <- rbind(all_models, models)
  print(i)
}

adjusted_models_uniform <- data.frame()
for (i in c(1:bootstrap_number)){
  data <- all_models %>% filter(replicate_number == i)
  data$p_adjusted <- p.adjust(data$p.value, method = "BH")
  adjusted_models_uniform <- rbind(adjusted_models_uniform, data)
}

#calculated sd
calculated_sd_samples <- generate_normal_samples(df = strep %>% dplyr::select(species_id, upperci, lowerci, predicted_d) %>% unique(), 
                                                 bootstrap_number = bootstrap_number, sd_type = "from_interval")

sample_data <- calculated_sd_samples %>% group_by(species_id) %>% mutate(replicate_number = 1:n()) %>%
  ungroup() %>% mutate(rate = log(2) / sample)

all_models <- data.frame()
for (i in strep_bcgs){
  presence <- bcg_presence_absence %>% filter(type == i)
  full_data <- sample_data %>% inner_join(., presence, by = "species_id")
  models <- full_data %>% dplyr::select(species_id, present, type, rate, replicate_number) %>%
    group_by(replicate_number) %>%
    nest() %>% 
    mutate(log_model = map(data, ~glm(cbind(present, 1 - present) ~ rate, data = ., family = "binomial"))) %>% 
    mutate(tidy = map(log_model, broom::tidy)) %>% 
    unnest(tidy) %>%
    dplyr::select(replicate_number, term, estimate, p.value) %>%
    filter(term != "(Intercept)") %>%
    ungroup() %>% mutate(type = i)
  all_models <- rbind(all_models, models)
  print(i)
}

adjusted_models_calcnorm <- data.frame()
for (i in c(1:bootstrap_number)){
  data <- all_models %>% filter(replicate_number == i)
  data$p_adjusted <- p.adjust(data$p.value, method = "BH")
  adjusted_models_calcnorm <- rbind(adjusted_models_calcnorm, data)
}

# set sd
deviation <- 0.25
preset_sd_samples <- generate_normal_samples(df = strep %>% dplyr::select(species_id, upperci, lowerci, predicted_d) %>% unique(), 
                                             bootstrap_number = bootstrap_number, sd_type = "preset", sd_provided = deviation)

sample_data <- preset_sd_samples %>% mutate(sample = abs(sample)) %>% group_by(species_id) %>% mutate(replicate_number = 1:n()) %>%
  ungroup() %>% mutate(rate = log(2) / sample)

all_models <- data.frame()
for (i in strep_bcgs){
  presence <- bcg_presence_absence %>% filter(type == i)
  full_data <- sample_data %>% inner_join(., presence, by = "species_id")
  models <- full_data %>% dplyr::select(species_id, present, type, rate, replicate_number) %>%
    group_by(replicate_number) %>%
    nest() %>% 
    mutate(log_model = map(data, ~glm(cbind(present, 1 - present) ~ rate, data = ., family = "binomial"))) %>% 
    mutate(tidy = map(log_model, broom::tidy)) %>% 
    unnest(tidy) %>%
    dplyr::select(replicate_number, term, estimate, p.value) %>%
    filter(term != "(Intercept)") %>%
    ungroup() %>% mutate(type = i)
  all_models <- rbind(all_models, models)
  print(i)
}

adjusted_models_setnorm <- data.frame()
for (i in c(1:bootstrap_number)){
  data <- all_models %>% filter(replicate_number == i)
  data$p_adjusted <- p.adjust(data$p.value, method = "BH")
  adjusted_models_setnorm <- rbind(adjusted_models_setnorm, data)
}

# part A - distribution of significant models
partA <- rbind(adjusted_models_calcnorm %>% mutate(distribution = "calculated normal"),
      adjusted_models_uniform %>% mutate(distribution = "uniform"),
      adjusted_models_setnorm %>% mutate(distribution = "preset normal")) %>%
  filter(p_adjusted < 0.05) %>%
  mutate(direction = ifelse(estimate > 0, "beta > 0", "beta < 0")) %>%
  group_by(replicate_number, direction, distribution) %>% summarize(n = n()) %>% ungroup() %>%
  ggplot(aes(n)) + geom_histogram(binwidth = 1) + facet_grid(distribution~direction, scales = "free") + theme_bw(base_size = 16) +
  xlab("# of significant models") + ylab("# of bootstrap replicates")

# positive calc norm - mean 20.9, median 21, max 26, min 16 
# negative calc norm - mean 5.6, median 5, max 9, min 4
# positive set norm - mean 10.7, median 11, max 21, min 1
# negative set norm - mean 3.6, median 4, max 8, min 1
# positive uni - mean 21.0, median 21, max 27, min 16
# negative uni - mean 5.6, median 5, max 9, min 4

# part B - 
partB <- rbind(adjusted_models_calcnorm %>% mutate(distribution = "calculated normal"),
               adjusted_models_uniform %>% mutate(distribution = "uniform"),
               adjusted_models_setnorm %>% mutate(distribution = "preset normal")) %>%
  filter(p_adjusted < 0.05) %>% mutate(direction = ifelse(estimate > 0, "beta > 0", "beta < 0")) %>%
  group_by(type, direction, distribution) %>% summarize(n = n()) %>%
  ungroup() %>% mutate(sig = ifelse(type %in% sig, "yes", "no")) %>%
  filter(((n / bootstrap_number) > 0.05) | sig == "yes") %>%
  mutate(type = case_when(type == "T1PKS" ~ "type I PKS",
                          type == "T3PKS" ~ "type III PKS",
                          type == "hserlactone" ~ "homoserine lactone",
                          type == "LAP" ~ "linear azol(in)e",
                          type == "hglE-KS" ~ "heterocyst glycolipid synthase-like PKS",
                          type == "NAPAA" ~ "non-alpha poly-amino acids",
                          type == "T2PKS" ~ "type II PKS",
                          type == "amglyccycl" ~ "aminoglycoside",
                          TRUE ~ type)) %>%
  ggplot(aes(x = fct_reorder(type, n), y = (n / bootstrap_number)*100, fill = sig)) + geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  coord_flip() + facet_grid(direction~distribution, scales = "free") + ylab("% of bootstrap replicates") + theme_bw(base_size = 16) + 
  scale_fill_manual(values = c("grey", "black")) +
  theme(axis.title.y = element_blank(), legend.position = "none")

#final figure
figure6 <- plot_grid(partA, partB, ncol = 2, label_size = 26, labels = c("A", "B"), rel_widths = c(0.5,1))

png(here::here("figures", "final-figs", "imgs", "figure-6.png"), res = 300, width = 4250, height = 2700)
figure6
dev.off()
