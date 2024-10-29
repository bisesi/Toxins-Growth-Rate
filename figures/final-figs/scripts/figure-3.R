# ATB
# figure 3
# bioinformatics shape data

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")
library("ggtext")

# load rocha data sets, note that 6 genomes do not have toxin data due to issues with antismash
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  dplyr::select(-c(`...1`)) %>% mutate(d_h = ifelse(species_id == "vibrio vulnificus", 20/60, d_h))
strep <- read_csv(here::here("bioinformatics-data", "strep", "full_growth_toxin_dataset_strep.csv")) %>%
  dplyr::select(-c(`...1`))

# clusters/toxins per genome histogram
partA <- rocha %>% dplyr::select(species_id, type) %>% group_by(species_id) %>% summarize(n = n()) %>%
  rbind(., rocha %>% filter(is.na(type)) %>% dplyr::select(species_id) %>% mutate(n = 0)) %>% mutate(dataset = "Viera-Silva & Rocha") %>%
  rbind(., strep %>% dplyr::select(species_id, type) %>% group_by(species_id) %>% summarize(n = n()) %>% mutate(dataset = "*Streptomyces*")) %>%
  mutate(dataset = factor(dataset, levels = c("Viera-Silva & Rocha", "*Streptomyces*"))) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  facet_wrap(~dataset, scales = "free") +
  ylab("# of genomes") +
  labs(fill = "")+
  xlab("# of clusters per genome") +
  theme_bw(base_size = 16) + theme(strip.text = element_markdown())

# growth rate histogram
partB <- rocha %>% dplyr::select(species_id, d_h) %>% unique() %>% mutate(dataset = "Viera-Silva & Rocha") %>%
  rbind(., strep %>% dplyr::select(species_id, predicted_d) %>% rename(d_h = predicted_d) %>% unique() %>% mutate(dataset = "*Streptomyces*")) %>%
  drop_na() %>%
  mutate(dataset = factor(dataset, levels = c("Viera-Silva & Rocha", "*Streptomyces*"))) %>%
  mutate(rate = log(2) / d_h) %>%
  ggplot(aes(rate)) +
  facet_wrap(~dataset, scales = "free") +
  geom_histogram() +
  ylab("# of genomes") +
  xlab("growth rate") +
  theme_bw(base_size = 16) + theme(strip.text = element_markdown())

# number of genes per type
partC <- rocha %>% filter(!is.na(type)) %>% mutate(total_number_of_genomes = rocha %>% pull(species_id) %>% unique() %>% length()) %>%
  dplyr::select(species_id, type, total_number_of_genomes) %>% unique() %>% group_by(type) %>% mutate(n = n()) %>% dplyr::select(-c(species_id)) %>% unique() %>%
  mutate(dataset = "Viera-Silva & Rocha") %>%
  rbind(., strep %>% mutate(total_number_of_genomes = strep %>% pull(species_id) %>% unique() %>% length()) %>%
          dplyr::select(species_id, type, total_number_of_genomes) %>% unique() %>% group_by(type) %>% mutate(n = n()) %>% dplyr::select(-c(species_id)) %>% unique() %>%
          mutate(dataset = "*Streptomyces*")) %>%
  mutate(percent = n / total_number_of_genomes) %>%
  mutate(dataset = factor(dataset, levels = c("Viera-Silva & Rocha", "*Streptomyces*"))) %>%
  group_by(dataset) %>%
  arrange(desc(percent)) %>%
  slice_head(n = 25) %>%
  ggplot(aes(x = reorder_within(as.character(substr(type, 1, 30)), percent, within = dataset), y = percent * 100)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~dataset, scales = "free") +
  theme_bw(base_size = 16) +
  scale_x_reordered() +
  xlab("toxin") +
  ylab("% of genomes") +
  theme(axis.title.y = element_blank()) + theme(strip.text = element_markdown())

#final figure
top <- plot_grid(partA, partB, ncol = 1, labels = c("A", "B"), label_size = 26)
figure3 <- plot_grid(top, partC, labels = c("", "C"), label_size = 26, ncol = 2, rel_widths = c(0.5,1))

png(here::here("figures", "final-figs", "imgs", "figure-3.png"), res = 300, width = 4250, height = 1500)
figure3
dev.off()

