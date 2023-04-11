#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Database Developmental plasticity thermal in reptiles 
##### Script to input data, modify and save plasticity value

### PACKAGES ####
library(tidyverse)
library(effectsize)
library(diversitree)
library(here)
library(rotl)
library(ape)
library(taxize)

### READING DATA ###

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

### To calculate Hedge's we need a group with studies, species, and traits ####

### FUNCTIONS ###
temp_ampli <- function(x, f){
   x %>%
    group_by(id, paper_no, species_complete, simp_trait) %>%
    filter(T == f(T)) %>%
    arrange(id)
}

mean_calculation <- function(x, f){
  x %>%
    select(-c(T)) %>%
    summarise(mean_variable = f(mean))# %>%
    #mutate(mean_variable = abs(mean_variable))
}

hed_g <- function(mean1, mean2, sd_p){
  round(((mean(mean1) - mean(mean2)) / sd_p), 3)
}

test_equal_temperature <- function(x, f){
  x %>%
    select(-c(mean)) %>%
    summarise(temperature = f(T))
}

#### RESULTS: TEMP, MEAN, SD POOL, AND HEDGEs G (PLASTICITY) ####

min_values <- temp_ampli(subdata, min)
max_values <- temp_ampli(subdata, max)

mean_min_values <- mean_calculation(min_values, mean)
mean_max_values <- mean_calculation(max_values, mean)

result <- mean_min_values %>%
  select(-c(mean_variable)) %>%
  mutate(sd_pool = NA)

for(i in seq_along(unique(max_values$id))){
  amostragem_1 <- min_values[min_values$id == i, ]
  amostragem_2 <- max_values[max_values$id == i, ]
  result$sd_pool[i] <- sd_pooled(amostragem_1$mean, amostragem_2$mean)
}

result$hedgesg <- mapply(hed_g, mean_min_values$mean_variable, 
                          mean_max_values$mean_variable,
                          result$sd_pool)

#### CLEAN DATA ####

min_test_values <- test_equal_temperature(min_values, mean)
max_test_values <- test_equal_temperature(max_values, mean)

is_equal <- function(x, y) x == y
result_equal <- mapply(is_equal, min_test_values$temperature, max_test_values$temperature)

result <- result[-c(which(result_equal)), ]

result <- result %>% 
  filter(!is.nan(hedgesg), !is.na(hedgesg), !is.infinite(hedgesg)) %>%
  select(-sd_pool)

#### Preparation to MUSSE ####
trait_frequent <- result %>%
  group_by(simp_trait) %>%
  summarise(len = length(simp_trait)) %>%
  arrange(desc(len))

result_all_species <- result %>% 
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))

result_mass <- result %>%
  group_by(species_complete) %>%
  filter(simp_trait == 'Mass') %>%
  summarise(hedgesg = mean(hedgesg))

#### phylogeny construction ####
species <- unique(result_all_species$species_complete)
species_mass <- unique(result_mass$species_complete)

species_info_ott <- tnrs_match_names(species)
species_tree_ott <- tol_induced_subtree(ott_ids = species_info_ott$ott_id)

#species_info_ott2 <- tnrs_match_names(species_mass)
#species_tree_ott2 <- tol_induced_subtree(ott_ids = species_info_ott2$ott_id)

species_info_ncbi <- classification(species, db='ncbi')
species_tree_ncbi <- class2tree(species_info_ncbi, check = TRUE)
species_tree_ncbi <- species_tree_ncbi$phylo

#species_info_ncbi2 <- classification(species_mass, db='ncbi')
#species_tree_ncbi2 <- class2tree(species_info_ncbi2, check = TRUE)
#species_tree_ncbi2 <- species_tree_ncbi2$phylo

plot(species_tree_ott)
plot(species_tree_ncbi)

#### MUSSE CALCULATION ###
### STATES ###
hedgesg <- abs(result_all_species$hedgesg)
hedgesg_mass <- abs(result_mass$hedgesg)

states <- function(x){
  if(x <= 0.2){
    x <- 0
  } else if (x > 0.2 & x <= 0.5){
    x <- 1
  } else if(x > 0.5 & x <= 0.8){
    x <- 2
  } else if(x > 0.8){
    x <- 3
  }
}
hedgesg <- sapply(hedgesg, states)
hedgesg_mass <- sapply(hedgesg_mass, states)

##### MUSSE ####
musse <- make.musse(species_tree_ncbi, states = hedgesg)
#resu_musse <- find.mle(musse, method="subplex")
#resu_musse$par