# ATB
# figure 2
# relationships between distance and biomass

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

#part A - coeff fitness
locations <- get_colony_locations(here::here("toxin-simulations", "simulations_spatial", "figure1", "starting_locations_figure1.csv"))
susceptibles <- data.frame()
for (i in unique(locations$spatial_seed)){
  data <- locations %>% filter(growth_rate == 0.25 & spatial_seed == i)
  sus <- get_susceptible_to_producer(data, "mean")
  susceptibles <- rbind(susceptibles, data.frame(sus, spatial_seed = i))
}

biomass <- get_colony_biomass(here::here("toxin-simulations", "simulations_spatial", "figure2", "total_biomass_figure2_coeff.csv"), toxin_coeff)
stationary <- get_stationary_timepoint(here::here("toxin-simulations", "simulations_spatial", "figure2", "total_biomass_figure2_coeff.csv"), toxin_coeff)

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
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_coeff), fill = mean)) +
  xlab("growth rate") +
  ylab("mmol toxin / gram biomass") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.5)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

partA2 <- coeff %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs cheater", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(fitness_type == "total pop") %>%
  ggplot(aes(x = as.factor(growth_rate), y = as.factor(toxin_coeff), fill = mean)) +
  xlab("growth rate") +
  ylab("mmol toxin / gram biomass") +
  geom_tile() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high= "darkgreen",
                       midpoint = 0.1)+
  facet_wrap(~fitness_type) +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none", strip.background = element_blank())

partA <- plot_grid(partA2, partA1, ncol = 2, rel_widths = c(1, 0.85))

# part B - toxin coeff susceptibles
distances_vs_biomass_susceptibles <- biomass %>% mutate(label = paste(species, colony_number, sep = "_")) %>% rename(biomass = value) %>%
  inner_join(., susceptibles %>% rename(distance = value), by = c("spatial_seed", "label")) %>%
  group_by(growth_rate, spatial_seed) %>% filter(cycle == max(cycle)) %>% ungroup()

slopes_coeff_susceptibles <- distances_vs_biomass_susceptibles %>%
  mutate(biomass = biomass * 1e7) %>%
  group_by(species, growth_rate, toxin_coeff, spatial_seed) %>%
  nest() %>% 
  mutate(linear_model = map(data, ~lm(biomass ~ distance, data = .))) %>% 
  mutate(tidy = map(linear_model, broom::tidy)) %>% 
  unnest(tidy) %>%
  dplyr::select(species, growth_rate, estimate, term, toxin_coeff, spatial_seed) %>%
  filter(term != "(Intercept)") %>%
  ungroup() 

partB <- slopes_coeff_susceptibles %>% group_by(growth_rate, toxin_coeff) %>% summarize(mean = mean(estimate), se = sd(estimate) / sqrt(3)) %>%
  mutate(type = "susceptibles") %>% filter(toxin_coeff %in% c(1,15,75)) %>%
  ggplot(aes(x = growth_rate, y = mean, shape = as.factor(toxin_coeff))) +
  geom_point(size = 3, stroke = 1.5) + theme_bw(base_size = 16) + geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_linerange(aes(ymin = mean - se, ymax = mean + se)) +
  scale_shape_manual(values = c(5,7,8))+
  ylab("effect of proximity") + xlab("growth rate") + theme(legend.position = "none", strip.background = element_blank()) + facet_wrap(~type)

legend2 <- get_plot_component(slopes_coeff_susceptibles %>% group_by(growth_rate, toxin_coeff) %>% summarize(mean = mean(estimate), se = sd(estimate) / sqrt(3)) %>%
                               mutate(type = "susceptibles") %>% filter(toxin_coeff %in% c(1,15,75)) %>%
                               ggplot(aes(x = growth_rate, y = mean, shape = as.factor(toxin_coeff))) +
                               geom_point(size = 3, stroke = 1.5) + theme_bw(base_size = 16) + geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
                               geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.01) +
                               scale_shape_manual(values = c(5,7,8))+ labs(shape = "mmol / gram")+
                               ylab("benefit of distance") + xlab("growth rate") + theme(legend.position = "bottom", strip.background = element_blank()) + facet_wrap(~type),
                             'guide-box-bottom', return_all = TRUE)

# part C toxin coeff resistants
resistant_biomass <- biomass %>% filter(species == "resistant") %>%
  group_by(growth_rate, spatial_seed,toxin_coeff) %>% filter(cycle == max(cycle)) %>% ungroup()

partC <- resistant_biomass %>% filter(toxin_coeff %in% c(1,15,75)) %>% group_by(toxin_coeff, growth_rate, spatial_seed) %>% 
  summarize(cv = sd(value) / mean(value)) %>% group_by(toxin_coeff, growth_rate) %>%
  mutate(mean_cv = mean(cv), se = sd(cv) / sqrt(3)) %>% 
  mutate(type = "resistant colony size") %>%
  ggplot(aes(x = growth_rate, y = mean_cv, shape = as.factor(toxin_coeff))) +
  geom_point(size = 3, stroke = 1.5) +
  scale_shape_manual(values = c(5,7,8))+
  theme_bw(base_size = 16) + xlab("growth rate") + ylab("coefficient of variation") +
  geom_linerange(aes(ymin = mean_cv - se, ymax = mean_cv + se)) +
  theme(legend.position = "none", strip.background = element_blank()) + facet_wrap(~type)

#final figure
figure2 <- plot_grid(partA, plot_grid(plot_grid(partB, partC, ncol = 2, label_size = 26, labels = c("B", "C")), legend2, rel_heights = c(1,0.1), ncol = 1),label_size = 26, labels = c("A",""), ncol = 2, rel_widths = c(1,1))

png(here::here("figures", "final-figs", "imgs", "figure-2.png"), res = 300, width = 4250, height = 1200)
figure2
dev.off()








