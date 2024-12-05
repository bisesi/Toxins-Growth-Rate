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

# clusters/toxins per genome histogram
partA <- levy %>% dplyr::select(genome, product_name) %>% group_by(genome) %>% summarize(n = n()) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  ylab("# of genomes") +
  labs(fill = "")+
  xlab("# of toxin genes per genome") +
  theme_bw(base_size = 16)

# clusters/toxins per genome histogram
partB <- levy %>% dplyr::select(genome, product_name) %>% unique() %>% group_by(genome) %>% summarize(n = n()) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_continuous(limits = c(0, NA)) +
  ylab("# of genomes") +
  labs(fill = "")+
  xlab("# of unique toxin types per genome") +
  theme_bw(base_size = 16)

# growth rate histogram
partC <- levy %>% dplyr::select(species_id, predicted_d) %>% unique() %>% 
  drop_na() %>%
  mutate(rate = log(2) / predicted_d) %>%
  ggplot(aes(rate)) +
  geom_histogram() +
  ylab("# of genomes") +
  xlab("growth rate") +
  theme_bw(base_size = 16)

# correlation with number of toxin genes 
d_vs_numbertoxgenes <- levy %>% group_by(genome) %>% mutate(n = n()) %>%
  ungroup() %>% dplyr::select(predicted_d, n) %>% unique() %>% filter(predicted_d < 200) %>%
  mutate(rate = log(2) / predicted_d) %>% dplyr::select(n, rate)

partD <- d_vs_numbertoxgenes %>%
  ggplot(aes(x = rate, y = n)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm") +
  ylab("# of toxin genes") +
  xlab("growth rate") +
  theme_bw(base_size = 16)

# correlation with unique types of toxin genes 
d_vs_numbertoxgenes_perproduct <- levy %>% unique() %>% group_by(genome) %>% mutate(n = n()) %>%
  ungroup() %>% dplyr::select(predicted_d, n) %>% unique() %>% filter(predicted_d < 200) %>%
  mutate(rate = log(2) / predicted_d)

partE <- d_vs_numbertoxgenes_perproduct %>%
  ggplot(aes(x = rate, y = n)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm") +
  ylab("# of unique toxin types") +
  xlab("growth rate") +
  theme_bw(base_size = 16)

# confidence intervals
bootstrap_total = do(10000)*lm(n ~ rate, data=mosaic::resample(d_vs_numbertoxgenes))
lower_total <- confint(bootstrap_total, level = 0.95)[2,2]
upper_total <- confint(bootstrap_total, level = 0.95)[2,3]
est_total <- coef(lm(n ~ rate, data = d_vs_numbertoxgenes))[2]

total <- data.frame(estimates = c(lower_total, upper_total, est_total), type = c("lower", "upper", "estimate"))

bootstrap_unique = do(10000)*lm(n ~ rate, data=mosaic::resample(d_vs_numbertoxgenes_perproduct))
lower_unique <- confint(bootstrap_unique, level = 0.95)[2,2]
upper_unique <- confint(bootstrap_unique, level = 0.95)[2,3]
est_unique <- coef(lm(n ~ rate, data = d_vs_numbertoxgenes_perproduct))[2]

unique <- data.frame(estimates = c(lower_unique, upper_unique, est_unique), type = c("lower", "upper", "estimate"))

partF <- rbind(total %>% mutate(comparison = "total"), unique %>% mutate(comparison = "unique")) %>%
  pivot_wider(names_from = type, values_from = estimates) %>%
  ggplot(aes(x = fct_reorder(comparison, estimate), y = estimate)) +
  geom_point(size = 3) +
  geom_pointrange(aes(ymax = upper, ymin = lower)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ylab("estimate") +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_blank())

#final figure
top <- plot_grid(partA, partB, partC, ncol = 3, labels = c("A", "B", "C"), label_size = 26)
bottom <- plot_grid(partD, partE, partF, ncol = 3, labels = c("D", "E", "F"), label_size = 26, rel_widths = c(1,1,0.5))
suppfigure5 <- plot_grid(top, bottom, ncol = 1)

png(here::here("figures", "final-figs", "imgs", "supp-figure-5.png"), res = 300, width = 4000, height = 2500)
suppfigure5
dev.off()


