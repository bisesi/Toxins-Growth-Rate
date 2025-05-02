# ATB
# figure 3
# bioinformatics shape data

#load packages
library("tidyverse")
library("cowplot")
library("tidytext")
library("ggtext")

# load rocha data sets, note that 6 genomes do not have toxin data due to issues with antismash
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  dplyr::select(-c(`...1`)) %>% mutate(d_h = ifelse(species_id == "vibrio vulnificus", 20/60, d_h))
strep <- read_csv(here::here("bioinformatics-data", "strep", "full_growth_toxin_dataset_strep.csv")) %>%
  dplyr::select(-c(`...1`))

# clusters/toxins per genome histogram
partA1 <- rocha %>% dplyr::select(species_id, type) %>% group_by(species_id) %>% summarize(n = n()) %>%
  rbind(., rocha %>% filter(is.na(type)) %>% dplyr::select(species_id) %>% mutate(n = 0)) %>% mutate(dataset = "Vieira-Silva and Rocha") %>%
  rbind(., strep %>% dplyr::select(species_id, type) %>% group_by(species_id) %>% summarize(n = n()) %>% mutate(dataset = "*Streptomyces*")) %>%
  mutate(dataset = factor(dataset, levels = c("Vieira-Silva and Rocha", "*Streptomyces*"))) %>%
  filter(dataset == "Vieira-Silva and Rocha") %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  facet_wrap(~dataset, scales = "free") +
  ylab("Genomes") +
  labs(fill = "")+
  xlab("Clusters per genome") +
  theme_bw(base_size = 16) + theme(strip.text = element_markdown(), axis.text = element_text(colour="black"))

partA2 <- rocha %>% dplyr::select(species_id, type) %>% group_by(species_id) %>% summarize(n = n()) %>%
  rbind(., rocha %>% filter(is.na(type)) %>% dplyr::select(species_id) %>% mutate(n = 0)) %>% mutate(dataset = "Vieira-Silva and Rocha") %>%
  rbind(., strep %>% dplyr::select(species_id, type) %>% group_by(species_id) %>% summarize(n = n()) %>% mutate(dataset = "*Streptomyces*")) %>%
  mutate(dataset = factor(dataset, levels = c("Vieira-Silva and Rocha", "*Streptomyces*"))) %>%
  filter(dataset == "*Streptomyces*") %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  facet_wrap(~dataset, scales = "free") +
  ylab("Genomes") +
  labs(fill = "")+
  xlab("Clusters per genome") +
  theme_bw(base_size = 16) + theme(strip.text = element_markdown(), axis.title.y = element_blank(), axis.text = element_text(colour="black"))

png(here::here("figures", "final-figs", "imgs", "figure-3A1.png"), res = 300, width = 1000, height = 1000)
partA1
dev.off()

png(here::here("figures", "final-figs", "imgs", "figure-3A2.png"), res = 300, width = 1000, height = 1000)
partA2
dev.off()

# growth rate histogram
partB1 <- rocha %>% dplyr::select(species_id, d_h) %>% unique() %>% mutate(dataset = "Vieira-Silva and Rocha") %>%
  rbind(., strep %>% dplyr::select(species_id, predicted_d) %>% rename(d_h = predicted_d) %>% unique() %>% mutate(dataset = "*Streptomyces*")) %>%
  drop_na() %>%
  mutate(dataset = factor(dataset, levels = c("Vieira-Silva and Rocha", "*Streptomyces*"))) %>%
  filter(dataset == "Vieira-Silva and Rocha") %>%
  mutate(rate = log(2) / d_h) %>%
  ggplot(aes(rate)) +
  facet_wrap(~dataset, scales = "free") +
  geom_histogram() +
  ylab("Genomes") +
  xlab(expression(paste("Growth rate (", hr^{-1}, ")"))) +
  theme_bw(base_size = 16) + theme(strip.text = element_markdown(), axis.text = element_text(colour="black"))

partB2 <- rocha %>% dplyr::select(species_id, d_h) %>% unique() %>% mutate(dataset = "Vieira-Silva and Rocha") %>%
  rbind(., strep %>% dplyr::select(species_id, predicted_d) %>% rename(d_h = predicted_d) %>% unique() %>% mutate(dataset = "*Streptomyces*")) %>%
  drop_na() %>%
  mutate(dataset = factor(dataset, levels = c("Vieira-Silva and Rocha", "*Streptomyces*"))) %>%
  filter(dataset == "*Streptomyces*") %>%
  mutate(rate = log(2) / d_h) %>%
  ggplot(aes(rate)) +
  facet_wrap(~dataset, scales = "free") +
  geom_histogram() +
  ylab("Genomes") +
  xlab(expression(paste("Growth rate (", hr^{-1}, ")"))) +
  theme_bw(base_size = 16) + theme(strip.text = element_markdown(), axis.title.y = element_blank(),
                                   axis.text = element_text(colour="black"))

png(here::here("figures", "final-figs", "imgs", "figure-3B1.png"), res = 300, width = 1000, height = 1000)
partB1
dev.off()

png(here::here("figures", "final-figs", "imgs", "figure-3B2.png"), res = 300, width = 1000, height = 1000)
partB2
dev.off()

# number of genes per type
partC1 <- rocha %>% filter(!is.na(type)) %>% mutate(total_number_of_genomes = rocha %>% pull(species_id) %>% unique() %>% length()) %>%
  dplyr::select(species_id, type, total_number_of_genomes) %>% unique() %>% group_by(type) %>% mutate(n = n()) %>% dplyr::select(-c(species_id)) %>% unique() %>%
  mutate(dataset = "Vieira-Silva and Rocha") %>%
  rbind(., strep %>% mutate(total_number_of_genomes = strep %>% pull(species_id) %>% unique() %>% length()) %>%
          dplyr::select(species_id, type, total_number_of_genomes) %>% unique() %>% group_by(type) %>% mutate(n = n()) %>% dplyr::select(-c(species_id)) %>% unique() %>%
          mutate(dataset = "*Streptomyces*")) %>%
  mutate(percent = n / total_number_of_genomes) %>%
  mutate(dataset = factor(dataset, levels = c("Vieira-Silva and Rocha", "*Streptomyces*"))) %>%
  filter(dataset == "Vieira-Silva and Rocha") %>%
  group_by(dataset) %>%
  arrange(desc(percent)) %>%
  slice_head(n = 25) %>%
  mutate(type = case_when(type == "T1PKS" ~ "type I PKS",
                          type == "T3PKS" ~ "type III PKS",
                          type == "hserlactone" ~ "homoserine lactone",
                          type == "LAP" ~ "linear azol(in)e",
                          type == "hglE-KS" ~ "heterocyst glycolipid synthase-like PKS",
                          type == "NAPAA" ~ "non-alpha poly-amino acids",
                          type == "T2PKS" ~ "type II PKS",
                          TRUE ~ type)) %>%
  ggplot(aes(x = reorder_within(as.character(substr(type, 1, 40)), percent, within = dataset), y = percent * 100)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~dataset, scales = "free") +
  theme_bw(base_size = 16) +
  scale_x_reordered() +
  xlab("toxin") +
  ylab("Proportion of genomes (%)") +
  theme(axis.title.y = element_blank()) + theme(strip.text = element_markdown(), axis.text = element_text(colour="black"))

partC2 <- rocha %>% filter(!is.na(type)) %>% mutate(total_number_of_genomes = rocha %>% pull(species_id) %>% unique() %>% length()) %>%
  dplyr::select(species_id, type, total_number_of_genomes) %>% unique() %>% group_by(type) %>% mutate(n = n()) %>% dplyr::select(-c(species_id)) %>% unique() %>%
  mutate(dataset = "Vieira-Silva and Rocha") %>%
  rbind(., strep %>% mutate(total_number_of_genomes = strep %>% pull(species_id) %>% unique() %>% length()) %>%
          dplyr::select(species_id, type, total_number_of_genomes) %>% unique() %>% group_by(type) %>% mutate(n = n()) %>% dplyr::select(-c(species_id)) %>% unique() %>%
          mutate(dataset = "*Streptomyces*")) %>%
  mutate(percent = n / total_number_of_genomes) %>%
  mutate(dataset = factor(dataset, levels = c("Vieira-Silva and Rocha", "*Streptomyces*"))) %>%
  filter(dataset == "*Streptomyces*") %>%
  group_by(dataset) %>%
  arrange(desc(percent)) %>%
  slice_head(n = 25) %>%
  mutate(type = case_when(type == "T1PKS" ~ "type I PKS",
                          type == "T3PKS" ~ "type III PKS",
                          type == "hserlactone" ~ "homoserine lactone",
                          type == "LAP" ~ "linear azol(in)e",
                          type == "hglE-KS" ~ "heterocyst glycolipid synthase-like PKS",
                          type == "NAPAA" ~ "non-alpha poly-amino acids",
                          type == "T2PKS" ~ "type II PKS",
                          TRUE ~ type)) %>%
  ggplot(aes(x = reorder_within(as.character(substr(type, 1, 40)), percent, within = dataset), y = percent * 100)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~dataset, scales = "free") +
  theme_bw(base_size = 16) +
  scale_x_reordered() +
  xlab("toxin") +
  ylab("Proportion of genomes (%)") +
  theme(axis.title.y = element_blank()) + theme(strip.text = element_markdown(), axis.title.y = element_blank(),
                                                axis.text = element_text(colour="black"))

png(here::here("figures", "final-figs", "imgs", "figure-3C1.png"), res = 300, width = 1750, height = 1500)
partC1
dev.off()

png(here::here("figures", "final-figs", "imgs", "figure-3C2.png"), res = 300, width = 1750, height = 1500)
partC2
dev.off()

#final figure
top <- plot_grid(partA, partB, ncol = 1, labels = c("A", "B"), label_size = 26)
figure3 <- plot_grid(top, partC, labels = c("", "C"), label_size = 26, ncol = 2, rel_widths = c(0.5,1))

partA <- data.frame(x = 1, y = 1) %>% ggplot(aes(x = x, y = y)) +
  theme_void() + geom_blank()

place_holders <- 

png(here::here("figures", "final-figs", "imgs", "figure-3.png"), res = 300, width = 4750, height = 1500)
figure3
dev.off()

