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
subdata <- dados %>% select(paper_no, genus, species, simp_trait, T, mean)
subdata$species_complete <- paste0(subdata$genus," ", subdata$species)
subdata <- subdata %>% select(paper_no, species_complete,  simp_trait, T, mean)

### To calculate Hedge's we need a group with studies, species, and traits ####

### FUNCTIONS ###
temp_ampli <- function(x, f){
   x %>%
    as.tibble() %>%
    group_by(paper_no, species_complete, simp_trait) %>%
    filter(T == f(T))
}

length_temp <- function(x, f){
  x %>%
    summarise(length = f(T))
} 

mean_temp <- function(x, f){
  x %>%
    mutate(mean = as.double(mean)) %>%
    summarise(mean_variable = f(mean))
}

hedges <- function(mean_min, mean_max, group_1, group_2){
  abs((mean_min - mean_max) / sd_pooled(group_1, group_2))
}

### CALCULATION ###
max_values <- temp_ampli(subdata, max)
min_values <- temp_ampli(subdata, min)

length_max <- length_temp(max_values, length)
length_min <- length_temp(min_values, length)

mean_max_values <- mean_temp(max_values, mean)
mean_min_values <- mean_temp(min_values, mean)

result <- mean_max_values %>%
  select(-mean_variable)

result$hedges_g <- mapply(hedges, mean_min_values$mean_variable, mean_max_values$mean_variable, length_min$length, length_max$length)
