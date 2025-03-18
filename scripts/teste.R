library(tidyverse)
library(effectsize)
library(here)
library(rotl)
library(ape)
library(taxize)
library(diversitree)
library(RRphylo)
library(phytools)
library(ggplot2)
library(geiger)
library(randtip)
library(stableGR)
library(OUwie)
library(corHMM)
#library(bayou)
#library(ggtree)

dados <- read.csv("data/raw/Database.csv", header = TRUE, sep = ",") 

#### ENGINEER DATA AND MODELING ####
subdata <- dados %>% 
  as.tibble() %>%
  select(paper_no, genus, species, simp_trait, T, mean) %>%
  mutate(species_complete = paste0(genus," ", species)) %>%
  group_by(paper_no, species_complete, simp_trait) %>%
  mutate(id = paste0('ID', paper_no, species_complete, simp_trait)) %>%
  select(id, paper_no, species_complete, 
         simp_trait, T, mean)  %>%
  mutate(species_complete = str_squish(species_complete))

uniques <- subdata %>%
  ungroup() %>%
  select(id) %>%
  unique() %>%
  unlist()

for(i in seq_along(uniques)){
  temporary <- subdata$id == uniques[[i]]
  subdata$id[temporary] <- i
}

subdata <- subdata %>% 
  mutate(id = as.integer(id), 
         mean = as.double(mean))

trait_frequent <- subdata %>%
  group_by(simp_trait) %>%
  summarise(len = length(simp_trait)) %>%
  arrange(desc(len)) |>
  ungroup()

more_frequent_mass <- trait_frequent[grep("ass",trait_frequent$simp_trait),]
five_more_frequent_mass <- more_frequent_mass[1, 1]
subdata$units <- dados$units
subdata_filtered_mass <- subdata[subdata$simp_trait %in% unclass(five_more_frequent_mass)$simp_trait, ]

# filter and transform data to mass #
subdata_filtered_mass$units <- gsub("\\s", "", subdata_filtered_mass$units)
subdata_filtered_mass <- subdata_filtered_mass |>
  filter(units == "g" | units == "mg")
subdata_filtered_mass[subdata_filtered_mass$units == "mg", ]$mean <- subdata_filtered_mass[subdata_filtered_mass$units == "mg", ]$mean / 1000
data <- subdata_filtered_mass

data$T <- as.numeric(gsub(",", ".", as.character(data$T)))
data$mean <- as.numeric(gsub(",", ".", as.character(data$mean)))

# Define function to calculate Hedges' g safely
hedges_g_safe <- function(mean1, mean2, sd1, sd2, n1, n2) {
  # Check for insufficient data
  if (n1 < 2 | n2 < 2 | is.na(sd1) | is.na(sd2)) return(NA)
  
  # Pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Calculate Hedges' g
  if (pooled_sd > 0) {
    g <- (mean1 - mean2) / pooled_sd
  } else {
    g <- NA
  }
  return(g)
}

# Group by id, paper_no, and species_complete, calculate min and max T and corresponding traits
hedges_g_results <- data %>%
  group_by(id, paper_no, species_complete) %>%
  summarise(
    T_min = min(T),
    T_max = max(T),
    mean_trait_min = mean[which.min(T)],
    mean_trait_max = mean[which.max(T)],
    std_trait_min = sd(mean[T == T_min]),
    std_trait_max = sd(mean[T == T_max]),
    n_trait_min = sum(T == T_min),
    n_trait_max = sum(T == T_max)
  ) %>%
  rowwise() %>%
  mutate(
    hedges_g = hedges_g_safe(mean_trait_min, mean_trait_max, std_trait_min, std_trait_max, n_trait_min, n_trait_max)
  ) %>%
  ungroup()

hedges_g_results$hedges_g <- abs(hedges_g_results$hedges_g)

species_mean_hedges_g <- hedges_g_results %>%
  group_by(species_complete) %>%
  summarise(mean_hedges_g = mean(hedges_g, na.rm = TRUE))

# View the result
print(species_mean_hedges_g)

hedgesg <- abs(species_mean_hedges_g$mean_hedges_g)

hedgesg <- setNames(hedgesg, species_mean_hedges_g$species_complete)

rates <- read.csv2("data/raw/alldat.csv", sep=',')

rates$treename <- gsub("_", " ", rates$treename)

states <- data.frame(hedgesg)
states <- rownames_to_column(states, var='species')
states$species <- tolower(states$species)

rates <- rates |>
  select(treename, bamm, clads, massRate)

names(rates)[1] <- "species"
rates$species <- tolower(rates$species)

result_final_bamm <- inner_join(states, rates, by='species')
result_final_bamm$bamm <- as.numeric(result_final_bamm$bamm)
result_final_bamm$clads <- as.numeric(result_final_bamm$clads)
result_final_bamm$massRate <- as.numeric(result_final_bamm$massRate)

plot(log(result_final_bamm$hedgesg), log(result_final_bamm$bamm))
plot(log(result_final_bamm$hedgesg), log(result_final_bamm$clads))
plot(log(result_final_bamm$hedgesg), log(result_final_bamm$massRate))

plot(result_final_bamm$hedgesg, result_final_bamm$bamm)
plot(result_final_bamm$hedgesg, result_final_bamm$clads)
plot(result_final_bamm$hedgesg, result_final_bamm$massRate)


### todos os traços ###
dados <- read.csv("data/raw/Database.csv", header = TRUE, sep = ",") 

subdata <- dados %>% 
  as.tibble() %>%
  select(paper_no, genus, species, simp_trait, T, mean) %>%
  mutate(species_complete = paste0(genus," ", species)) %>%
  group_by(paper_no, species_complete, simp_trait) %>%
  mutate(id = paste0('ID', paper_no, species_complete, simp_trait)) %>%
  select(id, paper_no, species_complete, 
         simp_trait, T, mean)  %>%
  mutate(species_complete = str_squish(species_complete))

uniques <- subdata %>%
  ungroup() %>%
  select(id) %>%
  unique() %>%
  unlist()

for(i in seq_along(uniques)){
  temporary <- subdata$id == uniques[[i]]
  subdata$id[temporary] <- i
}

subdata <- subdata %>% 
  mutate(id = as.integer(id), 
         mean = as.double(mean))

data <- subdata

data$T <- as.numeric(gsub(",", ".", as.character(data$T)))
data$mean <- as.numeric(gsub(",", ".", as.character(data$mean)))

# Define function to calculate Hedges' g safely
hedges_g_safe <- function(mean1, mean2, sd1, sd2, n1, n2) {
  # Check for insufficient data
  if (n1 < 2 | n2 < 2 | is.na(sd1) | is.na(sd2)) return(NA)
  
  # Pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Calculate Hedges' g
  if (pooled_sd > 0) {
    g <- (mean1 - mean2) / pooled_sd
  } else {
    g <- NA
  }
  return(g)
}

# Group by id, paper_no, and species_complete, calculate min and max T and corresponding traits
hedges_g_results <- data %>%
  group_by(id, paper_no, species_complete) %>%
  summarise(
    T_min = min(T),
    T_max = max(T),
    mean_trait_min = mean[which.min(T)],
    mean_trait_max = mean[which.max(T)],
    std_trait_min = sd(mean[T == T_min]),
    std_trait_max = sd(mean[T == T_max]),
    n_trait_min = sum(T == T_min),
    n_trait_max = sum(T == T_max)
  ) %>%
  rowwise() %>%
  mutate(
    hedges_g = hedges_g_safe(mean_trait_min, mean_trait_max, std_trait_min, std_trait_max, n_trait_min, n_trait_max)
  ) %>%
  ungroup()

hedges_g_results$hedges_g <- abs(hedges_g_results$hedges_g)

species_mean_hedges_g <- hedges_g_results %>%
  group_by(species_complete) %>%
  summarise(mean_hedges_g = mean(hedges_g, na.rm = TRUE))

# View the result
print(species_mean_hedges_g)

hedgesg <- abs(species_mean_hedges_g$mean_hedges_g)

hedgesg <- setNames(hedgesg, species_mean_hedges_g$species_complete)

rates <- read.csv2("data/raw/alldat.csv", sep=',')

rates$treename <- gsub("_", " ", rates$treename)

states <- data.frame(hedgesg)
states <- rownames_to_column(states, var='species')
states$species <- tolower(states$species)

rates <- rates |>
  select(treename, bamm, clads, massRate)

names(rates)[1] <- "species"
rates$species <- tolower(rates$species)

result_final_bamm <- inner_join(states, rates, by='species')
result_final_bamm$bamm <- as.numeric(result_final_bamm$bamm)
result_final_bamm$clads <- as.numeric(result_final_bamm$clads)
result_final_bamm$massRate <- as.numeric(result_final_bamm$massRate)

plot(log(result_final_bamm$hedgesg), log(result_final_bamm$bamm))
plot(log(result_final_bamm$hedgesg), log(result_final_bamm$clads))
plot(log(result_final_bamm$hedgesg), log(result_final_bamm$massRate))

plot(result_final_bamm$hedgesg, result_final_bamm$bamm)
plot(result_final_bamm$hedgesg, result_final_bamm$clads)
plot(result_final_bamm$hedgesg, result_final_bamm$massRate)
