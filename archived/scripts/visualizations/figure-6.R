# ATB
# figure 6
# levy bioinformatics analysis

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")
library("cowplot")

# load data
levy <- read_csv(here::here("bioinformatics-data", "levy", "full_growth_toxin_dataset_levy.csv")) %>%
  dplyr::select(-c(`...1`)) %>% rename(group = subgroup) %>% mutate(group = ifelse(is.na(group), "lysozyme", group)) %>%
  mutate(group = ifelse(group == "cytoxoxin", "cytotoxin", group)) %>% 
  filter(!group %in% c("muramidase", "small molecules", "trypsin", "subtisilin", "metallopeptidase", "siderophore", "ankyrin",
                       "transporter"))
  
# get bcg presence and absence
groups <- levy %>% pull(group) %>% unique()
tox_presence_absence <- data.frame()
for (i in groups){
  true_set <- levy %>% filter(group == i) %>% pull(genome) %>% unique()
  false_set <- setdiff(levy %>% pull(genome) %>% unique(), true_set)
  bool_0 <- data.frame(genome = false_set, present = as.integer(0), toxin = i)
  bool_1 <- data.frame(genome = true_set, present = as.integer(1), toxin = i)
  tox_presence_absence <- rbind(tox_presence_absence, bool_0, bool_1)
  print(paste(which(groups == i), "of", length(groups)))
}

# levy binomial
binomial_models_levy <- data.frame()
levy_bcgs <- levy %>% select(genome, predicted_d) %>% unique() %>%
  inner_join(., tox_presence_absence, by = "genome") %>% 
  filter(present == 1) %>% group_by(toxin) %>% summarize(n = n()) %>% filter(n > 3) %>%
  pull(toxin) %>% unique()

for (i in levy_bcgs){
  print(i)
  data <- levy %>% select(genome, predicted_d) %>% unique() %>% 
    inner_join(., tox_presence_absence, by = "genome") %>% filter(toxin == i)
  model <- glm(data$present ~ log2(data$predicted_d), family = "binomial")
  sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1])
  binomial_models_levy <- rbind(binomial_models_levy, sum_table)
}

binomial_models_levy$p_adjusted <- p.adjust(binomial_models_levy$p_value, method = "BH")

# part A
partA <- binomial_models_levy %>% mutate(sig = ifelse(p_adjusted < 0.05, "significant", "non-significant")) %>% 
  mutate(corre = ifelse(beta < 0, "negative", "positive")) %>% group_by(sig, corre) %>% summarize(n = n()) %>%
  ggplot(aes(x = corre, y = n, fill = sig)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  theme_bw(base_size = 16) +
  xlab("sign of beta") +
  ylab("number of toxin groups") +
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "none")

# part B
partB <- binomial_models_levy %>%
  filter(p_adjusted < 0.05 & beta < 0 & type %in% well_represented_groups) %>%
  mutate(significant = ifelse(p_adjusted < 0.05, TRUE, FALSE)) %>%
  ggplot(aes(x = fct_reorder(type, beta), y = beta, fill = significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 16) +
  xlab("BGC") +
  ylab("beta") +
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "none")

# part C
sig_bcgs_levy <- binomial_models_levy  %>%
  filter(p_adjusted < 0.05 & beta < -0.5 & type %in% well_represented_groups) %>% arrange(beta) %>% pull(type)

partC <- levy %>% select(genome, predicted_d) %>% unique() %>%
  inner_join(., tox_presence_absence, by = "genome") %>%
  filter(toxin %in% sig_bcgs_levy) %>%
  ggplot(aes(x = log2(predicted_d), y = present)) +
  scale_y_discrete(limits = c(0, 1))+
  geom_point(shape = 1) +
  facet_wrap(~factor(toxin, levels = sig_bcgs_levy), ncol = 3) +
  theme_bw(base_size = 16) +
  xlab("log2(predicted doubling time (hours))") +
  ylab("BGC presence") +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

#final figure
left <- plot_grid(partA, partB, ncol = 1, labels = c("A", "B"), label_size = 26)
figure6 <- plot_grid(left, partC, labels = c("", "C"), label_size = 26, ncol = 2, rel_widths = c(1, 0.7))

png(here::here("figures", "final-figs", "imgs", "figure-6.png"), res = 300, width = 4500, height = 4500)
figure6
dev.off()

