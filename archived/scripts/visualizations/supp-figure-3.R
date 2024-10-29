# ATB
# supplemental figure 3
# rocha bioinformatics analysis, truncated d

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")

# load data
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  dplyr::select(-c(`...1`))

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
rocha_bcgs <- rocha %>% select(species_id, predicted_d) %>% unique() %>% inner_join(., bcg_presence_absence, by = "species_id") %>% 
  filter(present == 1) %>% group_by(type) %>% summarize(n = n()) %>% filter(n > 3) %>%
  pull(type) %>% unique()

for (i in rocha_bcgs){
  print(i)
  data <- rocha %>% select(species_id, predicted_d) %>% unique() %>% inner_join(., bcg_presence_absence, by = "species_id") %>% filter(type == i) %>% 
    mutate(predicted_d = ifelse(predicted_d > 5, 5, predicted_d))
  model <- glm(data$present ~ log2(data$predicted_d), family = "binomial")
  sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1])
  binomial_models_rocha <- rbind(binomial_models_rocha, sum_table)
}

binomial_models_rocha$p_adjusted <- p.adjust(binomial_models_rocha$p_value, method = "BH")

# part A
partA <- binomial_models_rocha %>%
  mutate(significant = ifelse(p_adjusted < 0.05, TRUE, FALSE)) %>%
  ggplot(aes(x = fct_reorder(type, p_adjusted), y = p_adjusted, fill = significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 16) +
  xlab("BGC") +
  ylab("adjusted p-value") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed")+
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "none")

# part B
partB <- binomial_models_rocha %>%
  mutate(significant = ifelse(p_adjusted < 0.05, TRUE, FALSE)) %>%
  ggplot(aes(x = fct_reorder(type, beta), y = beta, fill = significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 16) +
  xlab("BGC") +
  ylab("beta") +
  geom_hline(yintercept = 0.0, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "none")

# part C
sig_bcgs_rocha <- binomial_models_rocha %>%
  filter(p_adjusted < 0.05) %>% arrange(beta) %>% pull(type)

partC <- rocha %>% select(species_id, predicted_d) %>% unique() %>%
  inner_join(., bcg_presence_absence, by = "species_id") %>%
  filter(type %in% sig_bcgs_rocha) %>%
  mutate(predicted_d = ifelse(predicted_d > 5, 5, predicted_d)) %>%
  ggplot(aes(x = log2(predicted_d), y = present)) +
  scale_y_discrete(limits = c(0, 1))+
  geom_point(shape = 1) +
  facet_wrap(~factor(type, levels = sig_bcgs_rocha), ncol = 2) +
  theme_bw(base_size = 16) +
  xlab("log2(predicted doubling time (hours))") +
  ylab("BGC presence") +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

#final figure
left <- plot_grid(partA, partB, ncol = 1, labels = c("A", "B"), label_size = 26)
suppfigure3 <- plot_grid(left, partC, labels = c("", "C"), label_size = 26, ncol = 2, rel_widths = c(1, 0.7))

png(here::here("figures", "final-figs", "imgs", "supp-figure-3.png"), res = 300, width = 3500, height = 3500)
suppfigure3
dev.off()


