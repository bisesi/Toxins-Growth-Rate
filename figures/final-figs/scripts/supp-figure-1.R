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
  ylab(expression(paste("Predicted growth rate (", hr^{-1}, ")"))) +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(0, NA)) +
  xlab(expression(paste("Observed growth rate (", hr^{-1}, ")"))) +
  theme_bw(base_size = 16)

# rocha experimental scatterplot
partB <- rocha %>% 
  mutate(pred_rate = log(2) / predicted_d, exp_rate = log(2) / d_h) %>%
  mutate(direction = case_when(round(pred_rate, 2) > round(exp_rate, 2) ~ "overestimate",
                               round(pred_rate, 2) < round(exp_rate, 2) ~ "underestimate", 
                               round(pred_rate, 2) == round(exp_rate, 2) ~ "accurate")) %>%
  select(pred_rate, exp_rate, species_id, direction) %>% unique() %>%
  group_by(direction) %>% summarize(n = n()) %>%
  ggplot(aes(x = direction, y = n)) + 
  geom_bar(stat = "identity") +
  ylab("Genomes") +
  theme_bw(base_size = 16) + theme(axis.title.x = element_blank())

# trend 
partC <- rocha %>% 
  mutate(pred_rate = log(2) / predicted_d, exp_rate = log(2) / d_h) %>% mutate(fold_diff = log2(pred_rate / exp_rate)) %>% 
  select(fold_diff, exp_rate, species_id, pred_rate) %>% unique() %>% 
  ggplot(aes(x = exp_rate, y = fold_diff)) + geom_point() +
  ylab(expression(paste("Log2(pred / obs growth rate(", hr^{-1}, "))"))) +
  xlab(expression(paste("Observed growth rate (", hr^{-1}, ")"))) +
  theme_bw(base_size = 16) + geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.14, color = "red", linetype = "dashed")

# rocha experimentally determined vs predict
partD <- rocha %>%
  dplyr::select(species_id, predicted_d) %>% unique() %>%
  mutate(pred_rate = log(2) / predicted_d) %>%
  ggplot(aes(pred_rate)) + 
  geom_histogram() +
  ylab("Genomes") +
  xlab(expression(paste("Predicted growth rate (", hr^{-1}, ")"))) +
  theme_bw(base_size = 16)

#final figure
suppfig1 <- plot_grid(partA, partB, partC, partD, ncol = 2, labels = c("A", "B", "C", "D"), label_size = 26)

png(here::here("figures", "final-figs", "imgs", "supp-figure-1.png"), res = 300, width = 3000, height = 3000)
suppfig1
dev.off()





