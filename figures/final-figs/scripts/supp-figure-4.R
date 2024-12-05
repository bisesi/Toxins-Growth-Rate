# ATB
# supp figure 4
# rocha bioinformatics analysis

#load packages
library("tidyverse")
library("cowplot")
library("tidytext")
library("ggtext")

# load data
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  dplyr::select(-c(`...1`)) %>% filter(!is.na(type)) %>% mutate(d_h = ifelse(species_id == "vibrio vulnificus", 20/60, d_h))

#part A
partA <- rocha %>% dplyr::select(species_id, type) %>% group_by(species_id) %>% summarize(n = n()) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  ylab("# of genomes") +
  labs(fill = "")+
  xlab("# of clusters per genome") +
  theme_bw(base_size = 16)

# growth rate histogram
partB <- rocha %>% dplyr::select(species_id, d_h) %>% unique() %>% mutate(dataset = "rocha") %>%
  drop_na() %>%
  ggplot(aes(log(2) / d_h)) +
  geom_histogram() +
  ylab("# of genomes") +
  xlab("growth rate") +
  theme_bw(base_size = 16)

# get bgc presence absence
all_bcgs <- rocha %>%
  filter(!is.na(type)) %>% pull(type) %>% unique()
bcg_presence_absence <- data.frame()
for (i in all_bcgs){
  true_set <- rocha %>% filter(type == i) %>% pull(species_id) %>% unique()
  false_set <- setdiff(rocha %>% pull(species_id) %>% unique(), true_set)
  bool_0 <- data.frame(species_id = false_set, present = as.integer(0), type = i)
  bool_1 <- data.frame(species_id = true_set, present = as.integer(1), type = i)
  bcg_presence_absence <- rbind(bcg_presence_absence, bool_0, bool_1)
}

# binomial models
binomial_models_rocha <- data.frame()
rocha_bcgs <- rocha %>% dplyr::select(species_id, d_h) %>% unique() %>% inner_join(., bcg_presence_absence, by = "species_id") %>% 
  filter(present == 1) %>%
  pull(type) %>% unique()

for (i in rocha_bcgs){
  print(i)
  data <- rocha %>% dplyr::select(species_id, d_h) %>% unique() %>% mutate(rate = log(2) / d_h) %>% inner_join(., bcg_presence_absence, by = "species_id") %>% filter(type == i)
  model <- glm(cbind(present, 1 - present) ~ rate, data = data, family = "binomial")
  sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1],
                          upperci = confint(model, level = 0.95, parm = "rate")[2], lowerci = confint(model, level = 0.95, parm = "rate")[1])
  binomial_models_rocha <- rbind(binomial_models_rocha, sum_table)
}

binomial_models_rocha$p_adjusted <- p.adjust(binomial_models_rocha$p_value, method = "BH")

# part C
partC <- binomial_models_rocha %>%
  mutate(significant = ifelse(p_adjusted <= 0.05, "sig", "ns"),
         sign = ifelse(beta >= 0, "positive\nbeta", "negative\nbeta")) %>%
  group_by(significant, sign) %>%
  summarize(n = n()) %>%
  rbind(., data.frame(significant = "sig", sign = "negative\nbeta", n = 0)) %>%
  rbind(., data.frame(significant = "sig", sign = "positive\nbeta", n = 0)) %>%
  ggplot(aes(x = sign, y = n, fill = significant)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  theme_bw(base_size = 16) +
  xlab("beta sign") +
  ylab("# of models") +
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "none", axis.title.x = element_blank())

#final figure
suppfigure4 <- plot_grid(partA, partB, partC, ncol = 3, labels = c("A", "B", "C"), label_size = 26)

png(here::here("figures", "final-figs", "imgs", "supp-figure-4.png"), res = 300, width = 3000, height = 1000)
suppfigure4
dev.off()

