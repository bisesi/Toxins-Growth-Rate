# ATB
# generate cleaned simulation datasets

#load packages
library("tidyverse")
library("stringdist")
library("cowplot")
library("tidytext")
library("reshape2")

# all stationary simulations - liquid
stationary_phase_conditions <- read_csv(here::here("toxin-simulations", "simulations_liquid", "raw", "liquid_biomass_file.csv")) %>% mutate(dataset = "biomass1") %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_liquid", "raw", "liquid_biomass_file_to500.csv")) %>% mutate(dataset = "biomass2")) %>%
  dplyr::select(-c(`...1`)) %>% mutate(total = producer + resistant + susceptible) %>% select(-c(producer, resistant, susceptible)) %>% 
  mutate(cycle = ifelse(cycle == 299 | cycle == 499,  "second_to_last", "last")) %>%
  pivot_wider(names_from = cycle, values_from = total) %>% filter(second_to_last == last)

total_stationary_dataset <- read_csv(here::here("toxin-simulations", "simulations_liquid", "raw", "liquid_biomass_file.csv")) %>% mutate(dataset = "biomass1") %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_liquid", "raw", "liquid_biomass_file_to500.csv")) %>% mutate(dataset = "biomass2")) %>%
  dplyr::select(-c(`...1`)) %>%
  full_join(., stationary_phase_conditions, by = c("growth_rate", "production_cost", "toxin_coefficient", "resistance_cost", "metabolite_diff", "add_signal_parameter", 
                                                    "space_widths", "dataset")) %>%
  filter(!cycle %in% c(299, 499)) %>% select(-c(second_to_last, last, dataset))

write.csv(total_stationary_dataset, file = here::here("toxin-simulations", "simulations_liquid", "aggregated", "all_stationary.csv"))

# producer or susceptible - liquid
full_data <- read_csv(here::here("toxin-simulations", "simulations_liquid", "raw", "liquid_biomass_file_susceptible_increase_noresist.csv")) %>% mutate(type = "susceptible growth rate") %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_liquid", "raw", "liquid_biomass_file_susceptible_increase_noresist.csv")) %>% mutate(type = "producer growth rate"))

write.csv(full_data, file = here::here("toxin-simulations", "simulations_liquid", "aggregated", "producer_susceptible.csv"))

# all stationary simulations - spatial
stationary_phase_conditions <- read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file.csv")) %>% mutate(dataset = "biomass1") %>%
  mutate(toxin_diff = 5e-7) %>% select(-c(metabolite_diff)) %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_tostationary.csv")) %>% mutate(dataset = "biomass2") %>%
          mutate(toxin_diff = 5e-7) %>% select(-c(metabolite_diff))) %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_diffusion.csv")) %>% mutate(dataset = "biomass3")) %>%
  dplyr::select(-c(`...1`)) %>% mutate(total = producer + resistant + susceptible) %>% select(-c(producer, resistant, susceptible)) %>% 
  mutate(cycle = ifelse(cycle == 49 | cycle == 149,  "second_to_last", "last")) %>%
  pivot_wider(names_from = cycle, values_from = total) %>% filter(second_to_last == last)

total_stationary_dataset <- read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file.csv")) %>% mutate(dataset = "biomass1") %>% 
  mutate(toxin_diff = 5e-7) %>% select(-c(metabolite_diff)) %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_tostationary.csv")) %>% mutate(dataset = "biomass2") %>%
          mutate(toxin_diff = 5e-7) %>% select(-c(metabolite_diff))) %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_diffusion.csv")) %>% mutate(dataset = "biomass3")) %>%
  dplyr::select(-c(`...1`)) %>%
  inner_join(., stationary_phase_conditions, by = c("growth_rate", "production_cost", "toxin_coefficient", "resistance_cost", "metabolite_diff", "add_signal_parameter", 
                                                   "replicate", "spatial_seed", "space_widths", "dataset")) %>%
  filter(!cycle %in% c(49, 149)) %>% select(-c(second_to_last, last, dataset))

write.csv(total_stationary_dataset, file = here::here("toxin-simulations", "simulations_spatial", "aggregated", "all_stationary.csv"))

# all location files
stationary_phase_conditions <- read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file.csv")) %>% mutate(dataset = "biomass1") %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_tostationary.csv")) %>% mutate(dataset = "biomass2")) %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_diffusion.csv")) %>% mutate(dataset = "biomass3")) %>%
  dplyr::select(-c(`...1`)) %>% mutate(total = producer + resistant + susceptible) %>% select(-c(producer, resistant, susceptible)) %>% 
  mutate(cycle = ifelse(cycle == 49 | cycle == 149,  "second_to_last", "last")) %>%
  pivot_wider(names_from = cycle, values_from = total) %>% filter(second_to_last == last)

locations <- read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_founder_file.csv")) %>% mutate(dataset = "biomass1") %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_founder_file_tostationary.csv")) %>% mutate(dataset = "biomass2")) %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_founder_file_diffusion.csv")) %>% mutate(dataset = "biomass3")) %>%
  dplyr::select(-c(`...1`)) %>%
  inner_join(., stationary_phase_conditions, by = c("growth_rate", "production_cost", "toxin_coefficient", "resistance_cost", "metabolite_diff", "add_signal_parameter", 
                                                    "replicate", "spatial_seed", "space_widths", "dataset"))

write.csv(locations, file = here::here("toxin-simulations", "simulations_spatial", "aggregated", "all_stationary_start_locations.csv"))

# producer or susceptible - spatial
full_data <- read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_susceptible_increase_noresist.csv")) %>% mutate(type = "susceptible growth rate") %>%
  rbind(., read_csv(here::here("toxin-simulations", "simulations_spatial", "raw", "spatial_biomass_file_susceptible_increase_noresist.csv")) %>% mutate(type = "producer growth rate"))

write.csv(full_data, file = here::here("toxin-simulations", "simulations_spatial", "aggregated", "producer_susceptible.csv"))

