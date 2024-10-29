# ATB
# figure 1
# example/intro figure

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

# part A sub spatial
locations <- get_colony_locations(here::here("toxin-simulations", "simulations_spatial", "figure1", "starting_locations_figure1.csv"))
  
partA_spatial <- locations %>% 
  filter(growth_rate == 0.125 & spatial_seed == 1) %>%
  ggplot(aes(x = x, y = y, color = strain)) + geom_point(size = 2) + 
  theme_bw(base_size = 16) + 
  scale_color_manual(values = c("producer" = "#004488", "susceptible" = "#DDAA33", "resistant" = "#BB5566")) +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), 
        strip.background = element_blank(), axis.ticks = element_blank(), strip.text.y = element_blank())

# part B - toxin localization
media <- read_csv(here::here("toxin-simulations", "simulations_spatial", "figure1", "metabolites_figure1.csv")) %>%
  dplyr::select(-c(`...1`))

metabolite_variation <- media %>% group_by(metabolite, growth_rate, spatial_seed, cycle) %>%
  filter(cycle > 1) %>%
  summarize(cv = sd(conc_mmol) / mean(conc_mmol)) %>% group_by(metabolite, growth_rate, spatial_seed) %>%
  filter(cv == max(cv)) %>% ungroup() %>% mutate(metabolite = ifelse(metabolite == "toxin_e", "toxin", "carbon")) %>%
  dplyr::select(-c(cycle))

partB <- metabolite_variation %>%
  group_by(growth_rate, metabolite) %>%
  summarize(mean_cv = mean(cv), se_cv = sd(cv) / sqrt(3)) %>%
  mutate(metabolite = factor(metabolite, levels = c("colony size", "carbon", "toxin"))) %>%
  filter(metabolite == "toxin") %>%
  ggplot(aes(x = growth_rate, y = mean_cv)) +
  geom_point(size = 4) + 
  theme_bw(base_size = 16) + 
  ylab("coefficient of variation") +
  geom_linerange(aes(ymin = mean_cv - se_cv, ymax = se_cv + mean_cv)) +
  xlab("growth rate") +
  theme(legend.position = "none", strip.background = element_blank()) + facet_wrap(~metabolite, scales = "free")


# part C - localization metric with producer
biomass <- get_colony_biomass(here::here("toxin-simulations", "simulations_spatial", "figure1", "total_biomass_figure1.csv"))

susceptibles <- data.frame()
for (i in unique(locations$spatial_seed)){
  data <- locations %>% filter(growth_rate == 0.25 & spatial_seed == i)
  sus <- get_susceptible_to_producer(data, "mean")
  susceptibles <- rbind(susceptibles, data.frame(sus, spatial_seed = i))
}

distances_vs_biomass <- biomass %>% mutate(label = paste(species, colony_number, sep = "_")) %>% rename(biomass = value) %>%
  inner_join(., susceptibles %>% rename(distance = value), by = c("spatial_seed", "label")) %>%
  group_by(growth_rate, spatial_seed) %>% filter(cycle == max(cycle)) %>% ungroup()

fast_slope <- coef(lm(biomass ~ distance, data = distances_vs_biomass %>% filter(growth_rate == 1 & spatial_seed == 5) %>% mutate(biomass = biomass * 1e7)))[2]
slow_slope <- coef(lm(biomass ~ distance, data = distances_vs_biomass %>% filter(growth_rate == 0.125 & spatial_seed == 5) %>% mutate(biomass = biomass * 1e7)))[2]

partC <- distances_vs_biomass %>% filter(growth_rate %in% c(0.125, 1) & spatial_seed == 5) %>% 
  mutate(type = "susceptibles") %>%
  ggplot(aes(x = distance, y = biomass * 1e7, color = as.factor(growth_rate))) +
  geom_point() + theme_bw(base_size = 16) + geom_smooth(method = "lm", se = FALSE) + labs(color = "growth rate") +
  scale_color_manual(values = c("darkgrey", "black")) +
  annotate(geom = "text", x = 35, y = 275, color = "black", size = 5, label = paste("slope:", round(fast_slope, 2))) +
  annotate(geom = "text", x = 43, y = 160, color = "darkgrey", size = 5, label = paste("slope:", round(slow_slope, 2))) +
  ylab("final scaled biomass") + xlab("mean distance to producer") + theme(legend.position = "none", strip.background = element_blank()) + facet_wrap(~type)

legend <- get_plot_component(distances_vs_biomass %>% filter(growth_rate %in% c(0.125, 1) & spatial_seed == 5) %>% 
                               mutate(type = "susceptibles") %>%
                               mutate(growth_rate = ifelse(growth_rate == 0.125, "slow", "fast")) %>%
                               ggplot(aes(x = distance, y = biomass * 1e5, color = as.factor(growth_rate))) +
                               geom_point(size = 4) + theme_bw(base_size = 16) + 
                               scale_color_manual(values = c("black", "darkgrey")) + labs(color = "")+
                               ylab("final scaled biomass") + xlab("mean distance to producer") + theme(legend.position = "bottom", strip.background = element_blank()) + facet_wrap(~type),
                             'guide-box-bottom', return_all = TRUE)

# part D - penalty of proximity
slopes <- distances_vs_biomass %>%
  mutate(biomass = biomass * 1e7) %>%
  dplyr::select(growth_rate, spatial_seed, biomass, distance) %>%
  group_by(growth_rate, spatial_seed) %>%
  nest() %>% 
  mutate(linear_model = map(data, ~lm(biomass ~ distance, data = .))) %>% 
  mutate(tidy = map(linear_model, broom::tidy)) %>% 
  unnest(tidy) %>%
  dplyr::select(growth_rate, estimate, term, spatial_seed) %>%
  filter(term != "(Intercept)") %>%
  ungroup() 

partD <- slopes %>% group_by(growth_rate) %>% summarize(mean = mean(estimate), se = sd(estimate) / sqrt(3)) %>%
  mutate(type = "susceptibles") %>%
  ggplot(aes(x = growth_rate, y = mean)) +
  geom_point(size = 4) + theme_bw(base_size = 16) + geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_linerange(aes(ymin = mean - se, ymax = mean + se)) +
  annotate(geom = "text", color = "black", x = 0.75, y = 0.1, label = "better to be far") +
  annotate(geom = "text", color = "black", x = 0.75, y = -0.1, label = "better to be close") +
  ylab("effect of proximity") + xlab("growth rate") + theme(legend.position = "none", strip.background = element_blank()) + facet_wrap(~type)

# part E - liquid vs spatial
liquid <- read_csv(here::here("toxin-simulations", "simulations_liquid", "liquid_biomass_file_figure.csv")) %>% 
  dplyr::select(-c(`...1`)) %>% filter(cycle %in% c(300, 500)) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         starting_percent = 1e-9 / (1e-9 + 1e-9 + 8e-9),
         final_percent = producer / (producer + resistant + susceptible)) %>% dplyr::select(growth_rate, relative_fitness, final_percent) %>%
  mutate(environment = "liquid",
         se_fitness = 0,
         se_pop = 0) %>%
  mutate(mean_pop = min(final_percent),
         mean_fitness = relative_fitness) %>% dplyr::select(-c(relative_fitness, final_percent))

stationary <- get_stationary_timepoint(here::here("toxin-simulations", "simulations_spatial", "figure1", "total_biomass_figure1.csv"))

spatial <- biomass %>% inner_join(., stationary, by = c("cycle", "growth_rate", "spatial_seed")) %>% 
  group_by(growth_rate, spatial_seed, species) %>% 
  summarize(total = sum(value)) %>% ungroup() %>% 
  pivot_wider(names_from = species, values_from = total) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
       final_percent = producer / (producer + resistant + susceptible)) %>%
  group_by(growth_rate) %>%
  summarize(mean_fitness = mean(relative_fitness),
            mean_pop = mean(final_percent),
            se_fitness = sd(relative_fitness) / sqrt(3),
            se_pop = sd(final_percent) / sqrt(3)) %>% mutate(environment = "spatial")

partE <- rbind(spatial, liquid) %>%
  pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
  mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs resistant", "total pop")) %>%
  mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
  mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  ggplot(aes(x = growth_rate, y = mean, shape = environment)) +
  xlab("growth rate") +
  ylab("producer fraction") +
  facet_wrap(~fitness_type, ncol = 2, scales = "free") +
  geom_point(size = 4) + 
  scale_shape_manual(values = c(17, 16))+
  theme_bw(base_size = 16) + 
  geom_linerange(aes(ymin = mean - se, ymax = se + mean)) +
  geom_hline(aes(yintercept = intercept), color = "red", linetype = "dashed") +
  theme(legend.position = "none", strip.background = element_blank())

legend2 <- get_plot_component(rbind(spatial, liquid) %>%
                                pivot_longer(cols = c(mean_fitness:se_pop), names_to = "name") %>%
                                mutate(fitness_type = ifelse(name == "mean_fitness" | name == "se_fitness", "vs cheater", "total pop")) %>%
                                mutate(intercept = ifelse(fitness_type == "total pop", 0.1, 0.5)) %>% 
                                mutate(stat = ifelse(name == "se_pop" | name == "se_fitness", "se", "mean")) %>% dplyr::select(-c(name)) %>%
                                pivot_wider(names_from = stat, values_from = value) %>%
                                filter(fitness_type == "total pop") %>%
                                ggplot(aes(x = growth_rate, y = mean, shape = environment)) +
                                xlab("growth rate") +
                                ylab("producer fraction") +
                                facet_wrap(~fitness_type, ncol = 2, scales = "free") +
                                geom_point(size = 4) + 
                                scale_shape_manual(values = c(17, 16))+
                                theme_bw(base_size = 16) + labs(shape = "") +
                                geom_linerange(aes(ymin = mean - se, ymax = se + mean)) +
                                geom_hline(aes(yintercept = intercept), color = "red", linetype = "dashed") +
                                theme(legend.position = "bottom", strip.background = element_blank()), 'guide-box-bottom', return_all = TRUE)

#blank space for part A
partA <- data.frame(x = 1, y = 1) %>% ggplot(aes(x = x, y = y)) +
  theme_void() + geom_blank()

#final figure
#top <- plot_grid(partA, plot_grid(partA_spatial, partA, ncol = 1, rel_heights = c(1,0.75)), plot_grid(partB, legend, ncol = 1, rel_heights = c(1,0.1)), ncol = 3, rel_widths = c(1, 0.4, 0.8), label_size = 26, labels = c("A", "", "B"))
top <- plot_grid(partA, partB, ncol = 2, rel_widths = c(1, 0.6), label_size = 26, labels = c("A", "B"))
bottom <- plot_grid(plot_grid(partC, legend, ncol = 1, rel_heights = c(1,0.1)), partD, plot_grid(partE, legend2, ncol = 1, rel_heights = c(1,0.1)), ncol = 3, rel_widths = c(0.6,0.5,1), label_size = 26, labels = c("C", "D", "E"))
figure1 <- plot_grid(top, bottom, ncol = 1)

png(here::here("figures", "final-figs", "imgs", "figure-1.png"), res = 300, width = 3500, height = 2000)
figure1
dev.off()


