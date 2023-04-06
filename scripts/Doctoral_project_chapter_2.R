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
library(here)

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
         simp_trait, T, mean) 

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
    summarise(mean_variable = f(mean))
}

hed_g <- function(mean1, mean2, sd_p){
  round(((mean(mean1) - mean(mean2)) / sd_p), 3)
}

min_values <- temp_ampli(subdata, min)
max_values <- temp_ampli(subdata, max)

mean_min_values <- mean_calculation(min_values, mean)
mean_max_values <- mean_calculation(max_values, mean)

result <- mean_min_values %>%
  select(-c(mean_variable))
result$sd_pool <- NA

for(i in seq_along(unique(max_values$id))){
  amostragem_1 <- min_values[min_values$id == i, ]
  amostragem_2 <- max_values[max_values$id == i, ]
  result$sd_pool[i] <- sd_pooled(amostragem_1$mean, amostragem_2$mean)
}

result$hedgesg <- mapply(hed_g, mean_min_values$mean_variable, 
                          mean_max_values$mean_variable,
                          result$sd_pool)
result <- result %>% 
  filter(!is.nan(hedgesg), !is.na(hedgesg), !is.infinite(hedgesg)) %>%
  select(-sd_pool)
