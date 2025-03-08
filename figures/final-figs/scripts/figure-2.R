# ATB
# figure 2
# relationships between distance and biomass

#load packages
library("tidyverse")
library("cowplot")
library("tidytext")

#source functions
source(here::here("functions", "functions.R"))

#part A - coeff fitness
biomass <- get_colony_biomass(here::here("toxin-simulations", "simulations_spatial", "figure2", "total_biomass_figure2_coeff.csv"), toxin_coeff) 
stationary <- biomass %>% select(-c(`...1`)) %>% group_by(toxin_coeff, growth_rate, spatial_seed) %>% filter(cycle == max(cycle)) %>% 
  select(cycle, growth_rate, spatial_seed, toxin_coeff) %>% unique()

coeff <- biomass %>% inner_join(., stationary, by = c("cycle", "growth_rate", "spatial_seed", "toxin_coeff")) %>% 
  group_by(growth_rate, spatial_seed, species, toxin_coeff) %>% 
  summarize(total = sum(value)) %>% ungroup() %>% 
  pivot_wider(names_from = species, values_from = total) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         final_percent = producer / (producer + resistant + susceptible)) %>%
  group_by(growth_rate, toxin_coeff) %>%
  summarize(mean_fitness = mean(relative_fitness),
            mean_pop = mean(final_percent),
            se_fitness = sd(relative_fitness) / sqrt(3),
            se_pop = sd(final_percent) / sqrt(3)) %>% mutate(environment = "spatial")

partA1 <- coeff %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "vs resistant") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(toxin_coeff %in% c(1, 15, 30, 50, 75)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_coeff), fill = mean)) +
  xlab("growth rate") +
  ylab("mmol toxin / gram biomass") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())

partA2 <- coeff %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs cheater", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "total pop") %>%
  mutate(fitness_type = "vs all") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(toxin_coeff %in% c(1, 15, 30, 50, 75)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_coeff), fill = mean)) +
  xlab("growth rate") +
  ylab("mmol toxin / gram biomass") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 0.1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank())

partA <- plot_grid(partA2, partA1, ncol = 2, rel_widths = c(1, 0.75))
partA_labeled <- ggdraw(add_sub(partA, "growth rate", vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5, size = 16))

# part B - toxin diffusion rate
biomass <- get_colony_biomass(here::here("toxin-simulations", "simulations_spatial", "figure2", "total_biomass_figure2_diff.csv"), toxin_diff)
stationary <- biomass %>% group_by(toxin_diff, growth_rate, spatial_seed) %>% filter(cycle == max(cycle)) %>% 
  select(cycle, growth_rate, spatial_seed, toxin_diff) %>% unique()

diff <- biomass %>% inner_join(., stationary, by = c("cycle", "growth_rate", "spatial_seed", "toxin_diff")) %>% 
  group_by(growth_rate, spatial_seed, species, toxin_diff) %>% 
  summarize(total = sum(value)) %>% ungroup() %>% 
  pivot_wider(names_from = species, values_from = total) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         final_percent = producer / (producer + resistant + susceptible)) %>%
  group_by(growth_rate, toxin_diff) %>%
  summarize(mean_fitness = mean(relative_fitness),
            mean_pop = mean(final_percent),
            se_fitness = sd(relative_fitness) / sqrt(3),
            se_pop = sd(final_percent) / sqrt(3)) %>% mutate(environment = "spatial")

partB1 <- diff %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "vs resistant") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(toxin_diff %in% c(5e-8, 5e-7, 5e-6, 2.5e-5, 5e-5)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_diff / 5e-6), fill = mean)) +
  xlab("growth rate") +
  ylab("toxin / metabolite diffusion rate") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())

partB2 <- diff %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  filter(fitness_type == "total pop" & name == "mean_pop") %>%
  mutate(fitness_type = "vs all") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(toxin_diff %in% c(5e-8, 5e-7, 5e-6, 2.5e-5, 5e-5)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_diff / 5e-6), fill = value)) +
  xlab("growth rate") +
  ylab("toxin / metabolite diffusion rate") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 0.1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank())

partB <- plot_grid(partB2, partB1, ncol = 2, rel_widths = c(1, 0.75))
partB_labeled <- ggdraw(add_sub(partB, "growth rate", vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5, size = 16))

# part C - density
biomass <- read_csv(here::here("toxin-simulations", "simulations_spatial", "figure2", "total_biomass_figure2_density.csv"))
stationary <- biomass %>% group_by(density, growth_rate, spatial_seed) %>% filter(cycle == max(cycle)) %>% 
  select(cycle, growth_rate, spatial_seed, density) %>% unique()

density <- biomass %>% inner_join(., stationary, by = c("cycle", "growth_rate", "spatial_seed", "density")) %>% 
  group_by(growth_rate, spatial_seed, variable, density) %>% 
  summarize(total = sum(value)) %>% ungroup() %>% 
  pivot_wider(names_from = variable, values_from = total) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         final_percent = producer / (producer + resistant + susceptible)) %>%
  group_by(growth_rate, density) %>%
  summarize(mean_fitness = mean(relative_fitness),
            mean_pop = mean(final_percent),
            se_fitness = sd(relative_fitness) / sqrt(3),
            se_pop = sd(final_percent) / sqrt(3)) %>% mutate(environment = "spatial")

partC1 <- density %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "vs resistant") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(density %in% c(1, 3, 5, 7, 9)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(density * 10), fill = mean)) +
  xlab("growth rate") +
  ylab("number of initial colonies") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())

partC2 <- density %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  filter(fitness_type == "total pop" & name == "mean_pop") %>%
  mutate(fitness_type = "vs all") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(density %in% c(1, 3, 5, 7, 9)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(density * 10), fill = value)) +
  xlab("growth rate") +
  ylab("number of initial colonies") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 0.1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank())

partC <- plot_grid(partC2, partC1, ncol = 2, rel_widths = c(1, 0.75))
partC_labeled <- ggdraw(add_sub(partC, "growth rate", vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5, size = 16))

#part D - resistant frequency fitness
biomass <- read_csv(here::here("toxin-simulations", "for_reviewers", "total_biomass_low_ratio.csv"))

stationary <- biomass %>% select(-c(`...1`)) %>% group_by(growth_rate, spatial_seed, prod_ratio) %>% filter(cycle == max(cycle)) %>% 
  select(cycle, growth_rate, spatial_seed, prod_ratio) %>% unique()

low_ratio <- biomass %>% inner_join(., stationary, by = c("cycle", "growth_rate", "spatial_seed", "prod_ratio")) %>% 
  group_by(growth_rate, spatial_seed, variable, prod_ratio) %>% 
  summarize(total = sum(value)) %>% ungroup() %>% 
  pivot_wider(names_from = variable, values_from = total) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         final_percent = producer / (producer + resistant + susceptible)) %>% 
  group_by(growth_rate, prod_ratio) %>%
  summarize(mean_fitness = mean(relative_fitness),
            mean_pop = mean(final_percent),
            se_fitness = sd(relative_fitness) / sqrt(3),
            se_pop = sd(final_percent) / sqrt(3)) %>% mutate(environment = "spatial")

biomass <- read_csv(here::here("toxin-simulations", "simulations_spatial", "figure2", "total_biomass_figure2_ratio.csv"))
stationary <- biomass %>% group_by(cycle, growth_rate, spatial_seed, prod_ratio) %>% summarize(total = sum(value)) %>%
  group_by(growth_rate, spatial_seed, prod_ratio) %>% filter(cycle == max(cycle)) %>% slice_head() %>% ungroup()

ratio <- biomass %>% inner_join(., stationary, by = c("cycle", "growth_rate", "spatial_seed", "prod_ratio")) %>% 
  group_by(growth_rate, spatial_seed, variable, prod_ratio) %>% 
  summarize(total = sum(value)) %>% ungroup() %>% 
  pivot_wider(names_from = variable, values_from = total) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         final_percent = producer / (producer + resistant + susceptible)) %>%
  group_by(growth_rate, prod_ratio) %>%
  summarize(mean_fitness = mean(relative_fitness),
            mean_pop = mean(final_percent),
            se_fitness = sd(relative_fitness) / sqrt(3),
            se_pop = sd(final_percent) / sqrt(3)) %>% mutate(environment = "spatial")

partD1 <- rbind(low_ratio, ratio) %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", prod_ratio, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "vs resistant") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(prod_ratio %in% c(0.1, 0.2, 0.3, 0.4, 0.06, 0.02)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(1 - 2*prod_ratio), fill = mean)) +
  xlab("growth rate") +
  ylab("starting susceptible frequency") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())

partD2 <- rbind(low_ratio, ratio) %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", value / prod_ratio, 0.5)) %>% 
  filter(fitness_type == "total pop" & name == "mean_pop") %>%
  mutate(fitness_type = "vs all") %>%
  filter(growth_rate %in% c(0.125, 1)) %>%
  filter(prod_ratio %in% c(0.1, 0.2, 0.3, 0.4, 0.06, 0.02)) %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(1 - 2*prod_ratio), fill = intercept)) +
  xlab("growth rate") +
  ylab("starting susceptible frequency") +
  geom_tile() +
  scale_fill_gradient2(low = "#AA4499",
                       mid = "white",
                       high= "#117733",
                       midpoint = 1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank())

partD <- plot_grid(partD2, partD1, ncol = 2, rel_widths = c(1, 0.75))
partD_labeled <- ggdraw(add_sub(partD, "growth rate", vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5, size = 16))

# legend
data <- data.frame(
  x = rep(1:5, each = 5),
  y = rep(1:5, times = 5),
  value = runif(25, min = -10, max = 10) # Numerical values
)

legend <- get_plot_component(ggplot(data, aes(x = x, y = y, fill = value)) +
                               geom_tile() +
                               scale_fill_gradient2(
                                 low = "#AA4499", mid = "white", high = "#117733", midpoint = 0,
                                 name = "", # Removes the default legend title
                                 labels = c("low\nfitness", "high\nfitness"), # Custom labels
                                 breaks = range(data$value)) +
                               theme_bw(base_size = 16),
                             'guide-box-right', return_all = TRUE)

#final figure
figure2 <- plot_grid(plot_grid(partA_labeled, partB_labeled, partC_labeled, partD_labeled, ncol = 2, label_size = 26, 
                               labels = c("A", "B", "C", "D")), legend, ncol = 2, rel_widths = c(1,0.1))

png(here::here("figures", "final-figs", "imgs", "figure-2.png"), res = 300, width = 4500, height = 2500)
figure2
dev.off()





