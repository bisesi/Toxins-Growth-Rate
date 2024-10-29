#load packages
library("tidyverse")
library("cowplot")
library("tidytext")
library("ranger")
library("kernlab")
library("drf")

# prepare final dataset where all conditions have reached stationary phase
spatial_dataset <- read_csv(here::here("toxin-simulations", "simulations_spatial", "aggregated", "all_stationary.csv")) %>%
  dplyr::select(-c(`...1`)) %>%
  mutate(relative_fitness = producer / (producer + resistant), 
         starting_percent = 1e-9 / (1e-9 + 1e-9 + 8e-9),
         final_percent = producer / (producer + resistant + susceptible),
         selected_for = ifelse(relative_fitness > 0.5, 1, 0),
         invaded = ifelse(final_percent > 0.1, 1, 0)) %>%
  select(-c(cycle, producer, resistant, susceptible, replicate, spatial_seed, relative_fitness, starting_percent, final_percent, metabolite_diff)) %>%
  unique()

# random forest on "selected_for"
set.seed(123)
num.trees <- 500
train_rows <- sample(1:nrow(spatial_dataset), nrow(spatial_dataset) * 0.75, replace = FALSE)
train <- spatial_dataset[train_rows,]
test <- spatial_dataset[-train_rows,]

against_cheats_forest <- ranger(factor(selected_for) ~., data = spatial_dataset %>% select(-c(invaded)), num.trees = num.trees, importance = 'impurity')
against_cheats_predictor_importance <- against_cheats_forest$variable.importance

invasion_forest <- ranger(factor(invaded) ~., data = spatial_dataset %>% select(-c(selected_for)), num.trees = num.trees, importance = 'impurity')
invasion_predictor_importance <- invasion_forest$variable.importance








predict.function <- function(model, new_observation) predict(model, new_observation, type="prob")[,2]

predictions <- data.frame()
for (i in c(1:nrow(spatial_dataset))){
  new_observation <- predict.function(model_rf, spatial_dataset[i,-7])
  br <- broken(model_rf, spatial_dataset[i,-7], data = spatial_dataset[,-7],
               predict.function = predict.function)
  new_frame <- data.frame(variable_name = br$variable_name, variable_value = br$variable_value, contribution = br$contribution) %>%
    mutate(variable_name = ifelse(variable_name == "", "final_prognosis", variable_name)) %>% filter(variable_name != "final_prognosis") %>% 
    mutate(total = sum(contribution)) %>% 
    filter(variable_name %in% c("Intercept", "production_cost", "growth_rate", "space_widths", "add_signal_parameter", "resistance_cost", "toxin_coefficient")) %>%
    mutate(observation = i)
  predictions <- rbind(predictions, new_frame)
  print(i)
}

write.csv(predictions, file = here::here("toxin-simulations", "simulations_spatial", "random_forest_predicts.csv"))

random_forest_contributions <- predictions %>%
  ggplot(aes(x = total, y = contribution, color = variable_value)) +
  geom_point() +
  facet_wrap(~variable_name, scales = "free") 

