# ATB
# supplemental figure 1
# grodon predictions

#load packages
library("tidyverse")
library("cowplot")
library("tidytext")
library("ggtext")

# load rocha data sets, note that 6 genomes do not have toxin data due to issues with antismash
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  dplyr::select(-c(`...1`)) %>% mutate(d_h = ifelse(species_id == "vibrio vulnificus", 20/60, d_h))

# rocha experimental scatterplot
partA <- rocha %>% 
  mutate(pred_rate = log(2) / predicted_d, exp_rate = log(2) / d_h) %>%
  ggplot(aes(x = exp_rate, y = pred_rate)) + 
  geom_point() + geom_abline(slope = 1, linetype = "dashed", color = "red") +
  ylab("predicted growth rate") +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(0, NA)) +
  xlab("observed growth rate") +
  theme_bw(base_size = 16)

# rocha experimental scatterplot under 5
partB <- rocha %>% 
  filter(predicted_d < 5 & d_h < 5) %>%
  unique() %>% 
  mutate(pred_rate = log(2) / predicted_d, exp_rate = log(2) / d_h) %>%
  ggplot(aes(x = exp_rate, y = pred_rate)) + 
  geom_point() + geom_abline(slope = 1, linetype = "dashed", color = "red") +
  ylab("predicted growth rate") +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(0, NA)) +
  xlab("observed growth rate") +
  theme_bw(base_size = 16)

# rocha experimentally determined vs predict
partC <- rocha %>%
  dplyr::select(species_id, predicted_d) %>% unique() %>%
  mutate(pred_rate = log(2) / predicted_d) %>%
  ggplot(aes(pred_rate)) + 
  geom_histogram() +
  ylab("# of genomes") +
  xlab("predicted growth rate") +
  theme_bw(base_size = 16)

#final figure
suppfig1 <- plot_grid(partA, partB, partC, ncol = 3, labels = c("A", "B", "C"), label_size = 26)

png(here::here("figures", "final-figs", "imgs", "supp-figure-1.png"), res = 300, width = 3000, height = 1250)
suppfig1
dev.off()





