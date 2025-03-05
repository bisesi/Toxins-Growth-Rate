# ATB
# supp figure 5
# levy correlations

library("tidyverse")
library("cowplot")
library("tidytext")
library("ggtext")
library("boot")
library("mosaic")

#load growth rate data
levy <- read_csv(here::here("bioinformatics-data", "levy", "full_growth_toxin_dataset_levy.csv")) %>%
  dplyr::select(-c(`...1`))

outliers <- c(637000000, 2593338849)

# clusters/toxins per genome histogram
partA <- levy %>% dplyr::select(genome, proteinID) %>% group_by(genome) %>% summarize(n = n()) %>%
  filter(!genome %in% outliers) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  ylab("# of genomes") +
  labs(fill = "")+
  xlab("# of toxin genes per genome") +
  theme_bw(base_size = 16)

# clusters/toxins per genome histogram
partB <- levy %>% dplyr::select(genome, pfamID) %>% unique() %>% group_by(genome) %>% summarize(n = n()) %>%
  filter(!genome %in% outliers) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  ylab("# of genomes") +
  labs(fill = "")+
  xlab("# of unique toxin types per genome") +
  theme_bw(base_size = 16)

# growth rate histogram
partC <- levy %>% filter(!genome %in% outliers) %>% 
  dplyr::select(species_id, predicted_d) %>% unique() %>% 
  drop_na() %>%
  mutate(rate = log(2) / predicted_d) %>%
  ggplot(aes(rate)) +
  geom_histogram() +
  ylab("# of genomes") +
  xlab("growth rate") +
  theme_bw(base_size = 16)

# correlation with number of toxin genes 
d_vs_numbertoxgenes <- levy %>% filter(!genome %in% outliers) %>% 
  select(genome, proteinID, predicted_d) %>% unique() %>% group_by(genome, predicted_d) %>% 
  summarize(n = n()) %>%
  ungroup() %>% mutate(rate = log(2) / predicted_d)

bootstrap_total = do(10000)*glm(n ~ rate, family = "poisson", data=mosaic::resample(d_vs_numbertoxgenes))
lower_total <- confint(bootstrap_total, level = 0.95)[2,2]
upper_total <- confint(bootstrap_total, level = 0.95)[2,3]
est_total <- coef(glm(n ~ rate, family = "poisson", data = d_vs_numbertoxgenes))[2]
sig = "***"

total <- data.frame(estimates = c(lower_total, upper_total, est_total), type = c("lower", "upper", "estimate"))

partD <- d_vs_numbertoxgenes %>%
  ggplot(aes(x = rate, y = n)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  ylab("# of toxin genes") +
  xlab("growth rate") +
  theme_bw(base_size = 16) +
  annotate("text", label = paste0("beta: ", round(total$estimates[3], 3), sig, " [", round(total$estimates[1], 3), " , ", round(total$estimates[2], 3), "]"), x = 4, y = 60)

# correlation with unique types of toxin genes 
d_vs_numbertoxgenes_perproduct <- levy %>% filter(!genome %in% outliers) %>% 
  select(genome, pfamID, predicted_d) %>% unique() %>% group_by(genome, predicted_d) %>% 
  summarize(n = n()) %>%
  ungroup() %>% mutate(rate = log(2) / predicted_d)

bootstrap_unique = do(1000)*glm(n ~ rate, family = "poisson", data=mosaic::resample(d_vs_numbertoxgenes_perproduct))
lower_unique <- confint(bootstrap_unique, level = 0.95)[2,2]
upper_unique <- confint(bootstrap_unique, level = 0.95)[2,3]
est_unique <- coef(glm(n ~ rate, family = "poisson", data = d_vs_numbertoxgenes_perproduct))[2]
sig = "***"

unique <- data.frame(estimates = c(lower_unique, upper_unique, est_unique), type = c("lower", "upper", "estimate"))

partE <- d_vs_numbertoxgenes_perproduct %>%
  ggplot(aes(x = rate, y = n)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  ylab("# of unique toxin types") +
  xlab("growth rate") +
  theme_bw(base_size = 16)+
  annotate("text", label = paste0("beta: ", round(unique$estimates[3], 3), sig, " [", round(unique$estimates[1], 3), " , ", round(unique$estimates[2], 3), "]"), x = 4, y = 60)

# pfamIDs
all_bcgs <- levy %>%
  filter(!is.na(pfamID)) %>% pull(pfamID) %>% unique()
bcg_presence_absence <- data.frame()
for (i in all_bcgs){
  true_set <- levy %>% filter(pfamID == i) %>% pull(genome) %>% unique()
  false_set <- setdiff(levy %>% pull(genome) %>% unique(), true_set)
  bool_0 <- data.frame(genome = false_set, present = as.integer(0), type = i)
  bool_1 <- data.frame(genome = true_set, present = as.integer(1), type = i)
  bcg_presence_absence <- rbind(bcg_presence_absence, bool_0, bool_1)
}

# binomial models
binomial_models_levy <- data.frame()
for (i in all_bcgs){
  print(paste(which(all_bcgs == i), "of", length(all_bcgs)))
  data <- levy %>% 
    dplyr::select(genome, predicted_d) %>% unique() %>% mutate(rate = log(2) / predicted_d) %>% 
    inner_join(., bcg_presence_absence, by = "genome") %>% 
    filter(type == i)
  model <- glm(cbind(present, 1 - present) ~ rate, data = data, family = "binomial")
  sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], pfamID = i, beta = summary(model)$coefficients[2,1])
  binomial_models_levy <- rbind(binomial_models_levy, sum_table)
}

binomial_models_levy$p_adjusted <- p.adjust(binomial_models_levy$p_value, method = "BH")

# partF, total models
partF <- binomial_models_levy %>%
  mutate(significant = ifelse(p_adjusted <= 0.05, "sig", "ns"),
         sign = ifelse(beta >= 0, "positive\nbeta", "negative\nbeta")) %>%
  group_by(significant, sign) %>%
  summarize(n = n()) %>%
  rbind(., data.frame(significant = "sig", sign = "negative\nbeta", n = 0)) %>%
  ggplot(aes(x = sign, y = n, fill = significant)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  theme_bw(base_size = 16) +
  xlab("beta sign") +
  ylab("# of models") +
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "none", axis.title.x = element_blank())

#final figure
top <- plot_grid(partA, partB, partC, ncol = 3, labels = c("A", "B", "C"), label_size = 26)
bottom <- plot_grid(partD, partE, partF, ncol = 3, labels = c("D", "E", "F"), label_size = 26, rel_widths = c(1,1,1))
suppfigure5 <- plot_grid(top, bottom, ncol = 1)

png(here::here("figures", "final-figs", "imgs", "supp-figure-5.png"), res = 300, width = 4000, height = 2500)
suppfigure5
dev.off()


