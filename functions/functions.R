# ATB
# Simulation analysis functions


# Function to find integer points within a Voronoi region
find_integer_points_in_neighborhood <- function(deldir_obj, point) {
  library("sp")
  tiles <- tile.list(deldir_obj)
  region <- NULL
  for (tile in tiles) {
    # Compare the tile's center point (tile$pt) with the given point
    if (round(tile$pt[1], 6) == point[1] && round(tile$pt[2], 6) == point[2]) {
      region <- tile
      break
    }
  }
  x_vertices <- region$x
  y_vertices <- region$y
  
  min_x <- floor(min(x_vertices))
  max_x <- ceiling(max(x_vertices))
  min_y <- floor(min(y_vertices))
  max_y <- ceiling(max(y_vertices))
  
  grid_points <- expand.grid(x = seq(min_x, max_x), y = seq(min_y, max_y))
  
  inside <- sp::point.in.polygon(grid_points$x, grid_points$y, x_vertices, y_vertices)
  
  return(grid_points[inside == 1, ])
}

# Function to get stationary timepoints
get_stationary_timepoint <- function(file, variable_changed){
  output <- read_csv(file) %>%
    #dplyr::select(-c(`...1`)) %>%
    mutate(biomass = dplyr::select(., matches(c("susceptible","resistant","producer"))) %>% rowSums(na.rm = TRUE)) %>%
    group_by(spatial_seed, growth_rate, {{variable_changed}}) %>%
    filter(biomass == max(biomass)) %>% slice_head() %>% ungroup() %>% 
    dplyr::select(cycle, spatial_seed, growth_rate, {{variable_changed}})
  return(output)
}

# Function to get colony biomass
get_colony_biomass <- function(file, variable_changed){
  output <- read_csv(file) %>%
    #dplyr::select(-c(`...1`)) %>% 
    pivot_longer(cols = matches(c("susceptible","resistant","producer"))) %>%
    mutate(species = case_when(grepl("susceptible", name) == TRUE ~ "susceptible",
                               grepl("resistant", name) == TRUE ~ "resistant",
                               grepl("producer", name) == TRUE ~ "producer",
                               TRUE ~ name)) %>% group_by(cycle, growth_rate, spatial_seed, species, {{variable_changed}}) %>% 
    mutate(colony_number = 1:n()) %>% ungroup()
  return(output)
}

# get all values in Voronoi neighborhoods
get_all_points_in_neighborhoods <- function(locations){
  full_data <- data.frame()
  for (i in unique(locations$spatial_seed)){
    dataset <- locations %>% filter(spatial_seed == i)
    tesselation <- deldir(dataset$x, dataset$y)
    points_inside_neighborhoods <- data.frame()
    for(j in 1:nrow(dataset)){
      x <- dataset[j,]$x
      y <- dataset[j,]$y
      in_neighborhood <- find_integer_points_in_neighborhood(tesselation, c(x,y))
      points_inside_neighborhoods <- rbind(data.frame(in_neighborhood, xi = x, yi = y, spatial_seed = i), points_inside_neighborhoods)
    }
    full_data <- rbind(full_data, points_inside_neighborhoods)
  }
  return(full_data)
}

#get all tesselations for a dataset
get_tesselations <- function(locations){
  full_data <- data.frame()
  for (i in unique(locations$spatial_seed)){
    dataset <- locations %>% filter(spatial_seed == i)
    tesselation <- deldir(dataset$x, dataset$y)
    output <- tesselation$summary %>% select(dir.area, x, y) %>% rename(xi = x, yi = y) %>% mutate(spatial_seed = i)
    full_data <- rbind(full_data, output)
  }
  return(full_data)
}

#get voronoi line segments to plot
get_voronoi_segments_for_ggplot <- function(locations){
  voronoi <- deldir(locations$x, locations$y)
  voronoi_edges <- voronoi$dirsgs
  voronoi_df <- data.frame(
    x1 = voronoi_edges$x1,
    y1 = voronoi_edges$y1,
    x2 = voronoi_edges$x2,
    y2 = voronoi_edges$y2
  )
  return(voronoi_df)
}

#get distances btween susceptibles and producers
get_susceptible_to_producer <- function(locations, metric){
  locs <- locations %>% filter(strain != "resistant") %>%
    group_by(strain) %>% ungroup() %>% mutate(label = paste(strain, colony_number, sep = "_"))
  d <- dist(as.matrix(locs)[,c("x", "y", "label")])
  d <- as.matrix(d, labels=TRUE)
  colnames(d) <- rownames(d) <- locs[['label']]
  if (metric == "mean"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("susceptible", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("producer"))) %>%
      pivot_longer(cols = matches(c("producer")), names_to = "comparison2") %>% group_by(comparison) %>% summarize(mean = mean(value)) %>% ungroup() %>% rename(value = mean)
  } else if (metric == "minimum"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("susceptible", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("producer"))) %>%
      pivot_longer(cols = matches(c("producer")), names_to = "comparison2") %>% group_by(comparison) %>% filter(value == min(value)) %>% slice_head() %>% ungroup() %>%
      dplyr::select(-c(comparison2))
  } else {print("that isn't currently supported")}
  distances <- d %>% rename(label = comparison) %>% unique()
  return(distances)
}

#get distances btween producer to susceptibles
get_producer_to_susceptibles <- function(locations, metric){
  locs <- locations %>% filter(strain != "resistant") %>%
    group_by(strain) %>% ungroup() %>% mutate(label = paste(strain, colony_number, sep = "_"))
  d <- dist(as.matrix(locs)[,c("x", "y", "label")])
  d <- as.matrix(d, labels=TRUE)
  colnames(d) <- rownames(d) <- locs[['label']]
  if (metric == "mean"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("producer", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("susceptible"))) %>%
      pivot_longer(cols = matches(c("susceptible")), names_to = "comparison2") %>% group_by(comparison) %>% summarize(mean = mean(value)) %>% ungroup() %>% rename(value = mean)
  } else if (metric == "minimum"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("producer", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("susceptible"))) %>%
      pivot_longer(cols = matches(c("susceptible")), names_to = "comparison2") %>% group_by(comparison) %>% filter(value == min(value)) %>% slice_head() %>% ungroup() %>%
      dplyr::select(-c(comparison2))
  } else {print("that isn't currently supported")}
  distances <- d %>% rename(label = comparison) %>% unique()
  return(distances)
}

#get distances btween resistants and producers
get_resistant_to_producer <- function(locations, metric){
  locs <- locations %>% filter(strain != "susceptible") %>%
    group_by(strain) %>% ungroup() %>% mutate(label = paste(strain, colony_number, sep = "_"))
  d <- dist(as.matrix(locs)[,c("x", "y", "label")])
  d <- as.matrix(d, labels=TRUE)
  colnames(d) <- rownames(d) <- locs[['label']]
  if (metric == "mean"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("resistant", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("producer"))) %>%
      pivot_longer(cols = c(producer_1:producer_5), names_to = "comparison2") %>% group_by(comparison) %>% summarize(mean = mean(value)) %>% ungroup() %>% rename(value = mean)
  } else if (metric == "minimum"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("resistant", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("producer"))) %>%
      pivot_longer(cols = c(producer_1:producer_5), names_to = "comparison2") %>% group_by(comparison) %>% filter(value == min(value)) %>% slice_head() %>% ungroup() %>%
      dplyr::select(-c(comparison2))
  } else {print("that isn't currently supported")} 
  distances <- d %>% rename(label = comparison) %>% unique()
  return(distances)
}

#get resistants to susceptibles
get_resistant_to_susceptibles <- function(locations, metric){
  locs <- locations %>% filter(strain != "producer") %>%
    group_by(strain) %>% ungroup() %>% mutate(label = paste(strain, colony_number, sep = "_"))
  d <- dist(as.matrix(locs)[,c("x", "y", "label")])
  d <- as.matrix(d, labels=TRUE)
  colnames(d) <- rownames(d) <- locs[['label']]
  if (metric == "mean"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("resistant", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("susceptible"))) %>%
      pivot_longer(cols = c(susceptible_1:susceptible_40), names_to = "comparison2") %>% group_by(comparison) %>% summarize(mean = mean(value)) %>% ungroup() %>% rename(value = mean)
  } else if (metric == "minimum"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("resistant", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("susceptible"))) %>%
      pivot_longer(cols = c(susceptible_1:susceptible_40), names_to = "comparison2") %>% group_by(comparison) %>% filter(value == min(value)) %>% slice_head() %>% ungroup() %>%
      dplyr::select(-c(comparison2))
  } else {print("that isn't currently supported")}
  distances <- d %>% rename(label = comparison) %>% unique()
  return(distances)
}

#get producer to resistants
get_producer_to_resistants <- function(locations, metric){
  locs <- locations %>% filter(strain != "susceptible") %>%
    group_by(strain) %>% ungroup() %>% mutate(label = paste(strain, colony_number, sep = "_"))
  d <- dist(as.matrix(locs)[,c("x", "y", "label")])
  d <- as.matrix(d, labels=TRUE)
  colnames(d) <- rownames(d) <- locs[['label']]
  if (metric == "mean"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("producer", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("resistant"))) %>%
      pivot_longer(cols = matches(c("resistant")), names_to = "comparison2") %>% group_by(comparison) %>% summarize(mean = mean(value)) %>% ungroup() %>% rename(value = mean)
  } else if (metric == "minimum"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("producer", comparison) == TRUE) %>% dplyr::select(c("comparison", starts_with("resistant"))) %>%
      pivot_longer(cols = matches(c("resistant")), names_to = "comparison2") %>% group_by(comparison) %>% filter(value == min(value)) %>% slice_head() %>% ungroup() %>%
      dplyr::select(-c(comparison2))
  } else {print("that isn't currently supported")}
  distances <- d %>% rename(label = comparison) %>% unique()
  return(distances)
}

#get producer to all others
get_producer_to_all <- function(locations, metric){
  locs <- locations %>% 
    group_by(strain) %>% ungroup() %>% mutate(label = paste(strain, colony_number, sep = "_"))
  d <- dist(as.matrix(locs)[,c("x", "y", "label")])
  d <- as.matrix(d, labels=TRUE)
  colnames(d) <- rownames(d) <- locs[['label']]
  if (metric == "mean"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("producer", comparison) == TRUE) %>% dplyr::select(c("comparison", matches(c("resistant", "susceptible", "producer")))) %>%
      pivot_longer(cols = matches(c("resistant", "susceptible", "producer")), names_to = "comparison2") %>% group_by(comparison) %>% summarize(mean = mean(value)) %>% ungroup() %>% rename(value = mean)
  } else if (metric == "minimum"){
    d <- data.frame(comparison = row.names(d), d) %>% filter(grepl("producer", comparison) == TRUE) %>% dplyr::select(c("comparison", matches(c("resistant", "susceptible", "producer")))) %>%
      pivot_longer(cols = matches(c("resistant", "susceptible", "producer")), names_to = "comparison2") %>% group_by(comparison) %>% filter(value == min(value)) %>% slice_head() %>% ungroup() %>%
      dplyr::select(-c(comparison2))
  } else {print("that isn't currently supported")}
  distances <- d %>% rename(label = comparison) %>% unique()
  return(distances)
}

# Function to get locations for a set of simulations
get_colony_locations <- function(file){
  output <- read_csv(file) %>%
    #dplyr::select(-c(`...1`)) %>% 
    group_by(strain, growth_rate, spatial_seed) %>%
    mutate(colony_number = 1:n()) %>% ungroup()
  return(output)
}

#get closest neighbor
number_of_closest_neighbors <- function(locations, species){
  locs <- locations %>% 
    group_by(strain) %>% ungroup() %>% mutate(label = paste(strain, colony_number, sep = "_"))
  d <- dist(as.matrix(locs)[,c("x", "y", "label")])
  d <- as.matrix(d, labels=TRUE)
  colnames(d) <- rownames(d) <- locs[['label']]
  d <- data.frame(comparison = row.names(d), d) %>% filter(grepl(species, comparison) == TRUE) %>% dplyr::select(c("comparison", matches(c("producer", "susceptible", "resistant")))) %>%
      pivot_longer(cols = matches(c("producer", "susceptible", "resistant")), names_to = "comparison2") %>% filter(value != 0) %>% 
      group_by(comparison) %>% filter(value == min(value)) %>% ungroup() %>% separate(comparison2, c("species", "number"), sep = "_") %>%
    dplyr::select(species, comparison) %>% unique() 
  return(d)
}


