
biomass <- read_csv(here::here("toxin-simulations", "simulations_spatial", "figure2", "starting_locations_figure2_diff.csv")) %>%
  dplyr::select(-c(`...1`))

biomass1 <- read_csv(here::here("toxin-simulations", "simulations_spatial", "figure2", "starting_locations_figure2_diff_1.csv")) %>%
  dplyr::select(-c(`...1`))

write_csv(rbind(biomass, biomass1), "starting_locations_figure2_diff.csv")


# based on predicted/observed from Rocha
rocha <- read_csv(here::here("bioinformatics-data", "rocha", "full_growth_toxin_dataset_rocha.csv")) %>%
  dplyr::select(-c(`...1`))

model <- lm(obs_rate ~ pred_rate, data = rocha %>% select(d_h, predicted_d) %>% 
              unique() %>% mutate(pred_rate = log(2) / predicted_d, obs_rate = log(2) / d_h))

new_strep <- strep %>% mutate(pred_rate = log(2) / predicted_d)

new_strep$from_model <- predict(model, newdata = new_strep %>% select(pred_rate))

binomial_models_strep <- data.frame()
no_ci <- data.frame()
strep_bcgs <- new_strep %>% dplyr::select(species_id, predicted_d) %>% unique() %>%
  inner_join(., bcg_presence_absence, by = "species_id") %>% 
  filter(present == 1) %>%
  pull(type) %>% unique()

for (i in strep_bcgs){
  print(i)
  skip_to_next <- FALSE
  data <- new_strep %>% dplyr::select(species_id, from_model) %>% unique() %>% mutate(rate = from_model) %>% 
    inner_join(., bcg_presence_absence, by = "species_id") %>% filter(type == i)
  model <- glm(cbind(present, 1 - present) ~ rate, data = data, family = "binomial")
  tryCatch({
    sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1],
                            upperci = confint(model, level = 0.95, parm = "rate")[2], lowerci = confint(model, level = 0.95, parm = "rate")[1])
    binomial_models_strep <- rbind(binomial_models_strep, sum_table)
  }, error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next){
    sum_table <- data.frame(p_value = summary(model)$coefficients[2,4], type = i, beta = summary(model)$coefficients[2,1],
                            upperci = NA, lowerci = NA)
    no_ci <- rbind(no_ci, sum_table)
  }
}

binomial_models_strep <- rbind(binomial_models_strep, no_ci)
binomial_models_strep$p_adjusted <- p.adjust(binomial_models_strep$p_value, method = "BH")

# part A
partA <- binomial_models_strep %>%
  mutate(significant = ifelse(p_adjusted <= 0.05, "sig", "ns"),
         sign = ifelse(beta >= 0, "positive\nbeta", "negative\nbeta")) %>%
  group_by(significant, sign) %>%
  summarize(n = n()) %>%
  rbind(., data.frame(significant = "sig", sign = "negative\nbeta", n = 0)) %>%
  ggplot(aes(x = sign, y = n, fill = significant)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  theme_bw(base_size = 16) +
  xlab("beta sign") +
  ylab("# of models") +
  scale_fill_manual(values = c("grey", "black")) +
  theme(legend.position = "none", axis.title.x = element_blank())

toxins <- c("RiPP-like", "NRPS", "NRPS-like", "butyrolactone", "cyanobactin",
            "amglyccycl", "blactam", "thiopeptide", "betalactone",
            "transAT-PKS-like", "furan", "LAP", "T2PKS", "T1PKS", "indole", "RRE-containing", "aminocoumarin",
            "lanthipeptide-class-ii", "lanthipeptide-class-iii", "lanthipeptide-class-iv")
# part B
partB <- binomial_models_strep %>%
  filter(p_adjusted < 0.05) %>%
  mutate(toxin = ifelse(type %in% toxins, TRUE, FALSE)) %>%
  mutate(type = case_when(type == "T1PKS" ~ "type I PKS",
                          type == "T3PKS" ~ "type III PKS",
                          type == "hserlactone" ~ "homoserine lactone",
                          type == "LAP" ~ "linear azol(in)e",
                          type == "hglE-KS" ~ "heterocyst glycolipid synthase-like PKS",
                          type == "NAPAA" ~ "non-alpha poly-amino acids",
                          type == "T2PKS" ~ "type II PKS",
                          type == "amglyccycl" ~ "aminoglycoside",
                          TRUE ~ type)) %>%
  ggplot(aes(x = fct_reorder(type, beta, .desc = TRUE), y = beta, color = toxin)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymax = upperci, ymin = lowerci)) +
  coord_flip() +
  theme_bw(base_size = 16) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ylab("beta") +
  scale_color_manual(values = c("black", "#E69F00")) +
  theme(legend.position = "none", axis.title.y = element_blank())

linear_model <- plot_grid(partA, partB)

