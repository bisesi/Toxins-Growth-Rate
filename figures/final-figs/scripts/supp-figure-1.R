# ATB
# supplemental figure 1
# parameter sweeps

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")
library("reshape2")
library("scales")
library("deldir")
library("ggforce")

#source functions
source(here::here("functions", "functions.R"))

#part A - resistant frequency fitness
biomass <- read_csv(here::here("toxin-simulations", "simulations_spatial", "suppfigure1", "total_biomass_figure2_ratio.csv")) %>%
  dplyr::select(-c(`...1`))
stationary <- biomass %>% group_by(cycle, growth_rate, spatial_seed, prod_ratio) %>% summarize(total = sum(value)) %>%
  group_by(growth_rate, spatial_seed, prod_ratio) %>% filter(total == max(total)) %>% slice_head() %>% ungroup()

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

partA1 <- ratio %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "vs resistant") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(1 - 2*prod_ratio), fill = mean)) +
  xlab("growth rate") +
  ylab("starting susceptible frequency") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

partA2 <- ratio %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", value / prod_ratio, 0.5)) %>% 
  filter(fitness_type == "total pop" & name == "mean_pop") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(1 - 2*prod_ratio), fill = intercept)) +
  xlab("growth rate") +
  ylab("starting susceptible frequency") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank())

partA <- plot_grid(partA2, partA1, ncol = 2, rel_widths = c(1, 0.85))

# diff
biomass <- get_colony_biomass(here::here("toxin-simulations", "simulations_spatial", "suppfigure1", "total_biomass_figure2_diff.csv"), toxin_diff)
stationary <- get_stationary_timepoint(here::here("toxin-simulations", "simulations_spatial", "suppfigure1", "total_biomass_figure2_diff.csv"), toxin_diff)

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
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_diff / 5e-6), fill = mean)) +
  xlab("growth rate") +
  ylab("toxin:metabolite diffusion rate") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

partB2 <- diff %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  filter(fitness_type == "total pop" & name == "mean_pop") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_diff / 5e-6), fill = value)) +
  xlab("growth rate") +
  ylab("toxin:metabolite diffusion rate") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank())

partB <- plot_grid(partB2, partB1, ncol = 2, rel_widths = c(1, 0.85))

# density
biomass <- get_colony_biomass(here::here("toxin-simulations", "simulations_spatial", "suppfigure1", "total_biomass_figure2_volume.csv"), width)
stationary <- get_stationary_timepoint(here::here("toxin-simulations", "simulations_spatial", "suppfigure1", "total_biomass_figure2_volume.csv"), width)

volume <- biomass %>% inner_join(., stationary, by = c("cycle", "growth_rate", "spatial_seed", "width")) %>% 
  group_by(growth_rate, spatial_seed, species, width) %>% 
  summarize(total = sum(value)) %>% ungroup() %>% 
  pivot_wider(names_from = species, values_from = total) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         final_percent = producer / (producer + resistant + susceptible)) %>%
  group_by(growth_rate, width) %>%
  summarize(mean_fitness = mean(relative_fitness),
            mean_pop = mean(final_percent),
            se_fitness = sd(relative_fitness) / sqrt(3),
            se_pop = sd(final_percent) / sqrt(3)) %>% mutate(environment = "spatial")

partC1 <- volume %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "vs resistant") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(width), fill = mean)) +
  xlab("growth rate") +
  ylab("environment size") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

partC2 <- volume %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  filter(fitness_type == "total pop" & name == "mean_pop") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(width), fill = value)) +
  xlab("growth rate") +
  ylab("environment size") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank())

partC <- plot_grid(partC2, partC1, ncol = 2, rel_widths = c(1, 0.85))

# density
biomass <- read_csv(here::here("toxin-simulations", "simulations_spatial", "suppfigure1", "total_biomass_figure2_density.csv")) %>%
  dplyr::select(-c(`...1`))
stationary <- biomass %>% group_by(cycle, growth_rate, spatial_seed, density) %>% summarize(total = sum(value)) %>%
  group_by(growth_rate, spatial_seed, density) %>% filter(total == max(total)) %>% slice_head() %>% ungroup()

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

partD1 <- density %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "vs resistant") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(density * 10), fill = mean)) +
  xlab("growth rate") +
  ylab("number of initial colonies") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

partD2 <- density %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  filter(fitness_type == "total pop" & name == "mean_pop") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(density * 10), fill = value)) +
  xlab("growth rate") +
  ylab("number of initial colonies") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank())

partD <- plot_grid(partD2, partD1, ncol = 2, rel_widths = c(1, 0.85))

#final figure
suppfigure1 <- plot_grid(partA, partB, partC, partD, ncol = 2, label_size = 26, labels = c("A", "B", "C", "D"))

png(here::here("figures", "final-figs", "imgs", "supp-figure-1.png"), res = 300, width = 4250, height = 2700)
suppfigure1
dev.off()



