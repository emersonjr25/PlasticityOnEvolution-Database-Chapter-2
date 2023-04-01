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
subdata <- dados %>% select(paper_no, genus, species, population, simp_trait, T, mean)
subdata$species_complete <- paste0(subdata$genus," ", subdata$species)

### To calculate Hedge's we need a list with studies, species, traits, temperature and traits ####

max_values <- subdata %>%
  as.tibble() %>%
  group_by(paper_no, species_complete, simp_trait) %>%
  filter(T == max(T))
  
min_values <- subdata %>%
  as.tibble() %>%
  group_by(paper_no, species_complete, simp_trait) %>%
  filter(T == min(T))

length_max <- max_values %>%
  summarise(length(T))

length_min <- min_values %>%
  summarise(length(T))

mean_max_values <- max_values %>%
  select(paper_no, 
         species_complete, 
         simp_trait,
         T,
         mean) %>%
  mutate(mean = as.double(mean)) %>%
  group_by(paper_no, 
           species_complete, 
           simp_trait) %>%
  summarise(mean(mean))

mean_min_values <- min_values %>%
  select(paper_no, 
         species_complete, 
         simp_trait,
         T,
         mean) %>%
  mutate(mean = as.double(mean)) %>%
  group_by(paper_no, 
           species_complete, 
           simp_trait) %>%
  summarise(mean(mean))

result <- mean_max_values %>%
  select(-`mean(mean)`)

hedges <- function(mean_min, mean_max, group_1, group_2){
   abs((mean_min - mean_max) / sd_pooled(group_1, group_2))
}

result$hedges_g <- mapply(hedges, mean_min_values$`mean(mean)`, mean_max_values$`mean(mean)`, length_min$`length(T)`, length_max$`length(T)`)
