# ATB
# supplemental figure 2
# rocha bioinformatics analysis

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")

# load data
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  dplyr::select(-c(`...1`)) %>% mutate(d_h = ifelse(species_id == "vibrio vulnificus", 20/60, d_h))

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
rocha_bcgs <- rocha %>% dplyr::select(species_id, predicted_d) %>% unique() %>% inner_join(., bcg_presence_absence, by = "species_id") %>% 
  filter(present == 1) %>% 
  pull(type) %>% unique()

for (i in rocha_bcgs){
  print(i)
  data <- rocha %>% dplyr::select(species_id, predicted_d) %>% unique() %>% mutate(rate = log(2) / predicted_d) %>% inner_join(., bcg_presence_absence, by = "species_id") %>% filter(type == i)
  model <- glm(cbind(present, 1 - present) ~ rate, data = data, family = "binomial")
  sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1],
                          upperci = confint(model, level = 0.95, parm = "rate")[2], lowerci = confint(model, level = 0.95, parm = "rate")[1])
  binomial_models_rocha <- rbind(binomial_models_rocha, sum_table)
}

binomial_models_rocha$p_adjusted <- p.adjust(binomial_models_rocha$p_value, method = "BH")

# part A
partA <- binomial_models_rocha %>%
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

# part B
partB <- binomial_models_rocha %>%
  filter(p_adjusted < 0.05) %>%
  ggplot(aes(x = fct_reorder(type, beta, .desc = TRUE), y = beta)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymax = upperci, ymin = lowerci)) +
  coord_flip() +
  theme_bw(base_size = 16) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ylab("beta") +
  theme(legend.position = "none", axis.title.y = element_blank())

# part C
sig_bcgs_rocha <- binomial_models_rocha %>%
  filter(p_adjusted < 0.05) %>% arrange(beta) %>% pull(type)

partC <- rocha %>% dplyr::select(species_id, predicted_d) %>% unique() %>%
  inner_join(., bcg_presence_absence, by = "species_id") %>%
  filter(type %in% sig_bcgs_rocha) %>%
  mutate(rate = log(2) / predicted_d) %>%
  ggplot(aes(x = log10(rate), y = present)) +
  scale_y_discrete(limits = c(0, 1))+
  geom_point(shape = 1) +
  facet_wrap(~factor(type, levels = sig_bcgs_rocha), ncol = 4) +
  theme_bw(base_size = 16) +
  xlab("predicted log10-transformed growth rate") +
  ylab("BGC presence") +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

#final figure
suppfig2 <- plot_grid(plot_grid(partA, partB, ncol = 2, labels = c("A", "B"), label_size = 26, rel_widths = c(0.8,1)), 
                      partC, ncol = 1, labels = c("", "C"), label_size = 26)

png(here::here("figures", "final-figs", "imgs", "supp-figure-2.png"), res = 300, width = 3000, height = 1500)
suppfig2
dev.off()


