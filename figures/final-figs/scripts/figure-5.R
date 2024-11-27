# ATB
# figure 5
# strep bioinformatics analysis

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")
library("ggtext")

# load data
strep <- read_csv(here::here("bioinformatics-data", "strep", "full_growth_toxin_dataset_strep.csv")) %>%
  dplyr::select(-c(`...1`))

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

# strep binomial
binomial_models_strep <- data.frame()
no_ci <- data.frame()
strep_bcgs <- strep %>% dplyr::select(species_id, predicted_d) %>% unique() %>%
  inner_join(., bcg_presence_absence, by = "species_id") %>% 
  filter(present == 1) %>%
  pull(type) %>% unique()

for (i in strep_bcgs){
  print(i)
  skip_to_next <- FALSE
  data <- strep %>% dplyr::select(species_id, predicted_d) %>% unique() %>% mutate(rate = log(2) / predicted_d) %>% 
    inner_join(., bcg_presence_absence, by = "species_id") %>% filter(type == i)
  model <- glm(cbind(present, 1 - present) ~ rate, data = data, family = "binomial")
  tryCatch({
    sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1],
                            upperci = confint(model, level = 0.95, parm = "rate")[2], lowerci = confint(model, level = 0.95, parm = "rate")[1])
    binomial_models_strep <- rbind(binomial_models_strep, sum_table)
  }, error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next){
    sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1],
                            upperci = NA, lowerci = NA)
    no_ci <- rbind(no_ci, sum_table)
  }
}

binomial_models_strep <- rbind(binomial_models_strep, no_ci)
binomial_models_strep$p_adjusted <- p.adjust(binomial_models_strep$p_value, method = "BH")

# part A
partA <- binomial_models_strep %>%
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

toxins <- c("RiPP-like", "NRPS", "NRPS-like", "butyrolactone", "cyanobactin",
            "amglyccycl", "blactam", "thiopeptide", "betalactone",
            "transAT-PKS-like", "furan", "LAP", "T2PKS", "T1PKS", "indole", "RRE-containing", "aminocoumarin",
            "lanthipeptide-class-ii", "lanthipeptide-class-iii", "lanthipeptide-class-iv")

robust <- c("NRPS", "RRE-containing", "RiPP-like", "aminoglycoside", "arylpolyene",
            "butrylactone", "homoserine lactone", "indole", "lanthipeptide-class-iv", "linear azol(in)e",
            "melanin", "non-alpha poly-amino acids", "other", "phosphoglycolipid", "ranthipeptide",
            "thiopeptide", "transAT-PKS-like", "type I PKS", "type II PKS")
# part B
partB <- binomial_models_strep %>%
  filter(p_adjusted < 0.05) %>%
  mutate(toxin = ifelse(type %in% toxins, TRUE, FALSE)) %>%
  mutate(type = case_when(type == "T1PKS" ~ "type I PKS",
                          type == "T3PKS" ~ "type III PKS",
                          type == "hserlactone" ~ "homoserine lactone",
                          type == "LAP" ~ "linear azol(in)e",
                          type == "hglE-KS" ~ "heterocyst glycolipid synthase-like PKS",
                          type == "NAPAA" ~ "non-alpha poly-amino acids",
                          type == "T2PKS" ~ "type II PKS",
                          type == "amglyccycl" ~ "aminoglycoside",
                          TRUE ~ type)) %>%
  mutate(type = ifelse(type %in% robust, paste0("***", type, "***"), type)) %>%
  ggplot(aes(x = fct_reorder(type, beta, .desc = TRUE), y = beta, color = toxin)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymax = upperci, ymin = lowerci)) +
  coord_flip() +
  theme_bw(base_size = 16) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ylab("beta") +
  scale_color_manual(values = c("black", "#E69F00")) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text = element_markdown())

# part C
pos_bcgs_strep <- c("type II PKS", "NRPS", "melanin", "RiPP-like")

partC <- strep %>% dplyr::select(species_id, predicted_d) %>% unique() %>%
  inner_join(., bcg_presence_absence, by = "species_id") %>%
  mutate(type = ifelse(type == "T2PKS","type II PKS",type)) %>%
  filter(type %in% pos_bcgs_strep) %>%
  mutate(rate = log(2) / predicted_d) %>%
  ggplot(aes(x = rate, y = present)) +
  scale_y_discrete(limits = c(0, 1))+
  geom_point(shape = 1) +
  facet_wrap(~factor(type, levels = pos_bcgs_strep), ncol = 6) +
  theme_bw(base_size = 16) +
  xlab("growth rate") +
  ylab("BGC presence") +
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

# part D
neg_bcgs_strep <- c("homoserine lactone", "RRE-containing", "arylpolyene", "indole")

partD <- strep %>% dplyr::select(species_id, predicted_d) %>% unique() %>%
  inner_join(., bcg_presence_absence, by = "species_id") %>%
  mutate(type = ifelse(type == "hserlactone","homoserine lactone",type)) %>%
  filter(type %in% neg_bcgs_strep) %>%
  mutate(rate = log(2) / predicted_d) %>%
  ggplot(aes(x = rate, y = present)) +
  scale_y_discrete(limits = c(0, 1))+
  geom_point(shape = 1) +
  facet_wrap(~factor(type, levels = neg_bcgs_strep), ncol = 6) +
  theme_bw(base_size = 16) +
  xlab("growth rate") +
  ylab("BGC presence") +
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  geom_smooth(method = "glm", method.args = list(family = "binomial"))


#final figure
figure5 <- plot_grid(plot_grid(partA, partB, ncol = 2, labels = c("A", "B"), label_size = 26, rel_widths = c(0.6,1)), 
                     partC, partD, ncol = 1, labels = c("", "C", "D"), label_size = 26, rel_heights = c(1, 0.5, 0.5))

png(here::here("figures", "final-figs", "imgs", "figure-5.png"), res = 300, width = 3000, height = 3000)
figure5
dev.off()
