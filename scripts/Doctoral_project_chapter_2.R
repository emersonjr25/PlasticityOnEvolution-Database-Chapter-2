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
subdata <- subdata %>% 
  mutate(id = 1:nrow(subdata)) %>%
  select(id, paper_no, species_complete, 
         simp_trait, T, mean)

### To calculate Hedge's we need a group with studies, species, and traits ####

### FUNCTIONS ###
temp_ampli <- function(x, f){
   x %>%
    as.tibble() %>%
    mutate(mean = as.double(mean)) %>%
    group_by(paper_no, species_complete, simp_trait) %>%
    filter(T == f(T))
}

sd_temp <- function(x, f){
  x %>%
    select(-T) %>%
    summarise(sd = f(T))
} 

mean_temp <- function(x, f){
  x %>%
    mutate(mean = as.double(mean)) %>%
    summarise(mean_variable = f(mean))
}

### CALCULATION ###
min_values <- temp_ampli(subdata, min)
max_values <- temp_ampli(subdata, max)

mean_max_values <- mean_temp(max_values, mean)
mean_min_values <- mean_temp(min_values, mean)

table_result <- min_values %>%
  ungroup() %>%
  mutate(mean= as.double(mean))

max_values <- max_values %>%
  ungroup() %>%
  select(id, T, mean)

table_result <- full_join(table_result, max_values, by ='id')
table_result <- table_result %>%
  rename(temp_min = T.x, mean_min = mean.x, 
         temp_max = T.y, mean_max = mean.y)

#### TRYING CALCULATE HEDGE's G EFFECT ###
sd_max <- max_values %>%
  select(-T) %>%
  mutate(mean = as.double(mean)) %>%
  summarise(sd = sd(mean))

sd_min <- min_values %>%
  select(-T) %>%
  mutate(mean = as.double(mean)) %>%
  summarise(sd = sd(mean))



sd_pool <- function(x, y){
  sqrt((x^2 + y^2) / 2)
}

result <- mean_max_values %>%
  select(-mean_variable)

result$sd_pool <- mapply(sd_pool, sd_min$sd, sd_max$sd)
result$hedges_g <- mapply(hedgesg, mean_min_values$mean_variable, mean_max_values$mean_variable, result$sd_pool)

