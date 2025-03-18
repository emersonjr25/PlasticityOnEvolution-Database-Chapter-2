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

################################################################################
### To calculate Hedge's we need a group with studies, species, and traits ####
################################################################################

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

min_values <- temp_ampli(subdata_filtered_mass, min)
max_values <- temp_ampli(subdata_filtered_mass, max)

mean_min_values <- mean_calculation(min_values, mean)
mean_max_values <- mean_calculation(max_values, mean)

mean_min_values <- ungroup(mean_min_values)
mean_max_values <- ungroup(mean_max_values)

#write.csv2(mean_min_values, "min.csv")
#write.csv2(mean_max_values, "max.csv")
#write.csv2(subdata_filtered_mass, "mass.csv")

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

papers_to_revise <- sample(unique(result$paper_no), 8)
sort(papers_to_revise)

### carrying all data ###
load("table_and_phy_ready.RDS")

### verifying more frequent trait ###

# mass per species #
table_bms <- subdata_filtered_mass |> 
  select(-c(units)) |>
  group_by(paper_no, species_complete) |>
  summarise(mean_value_trait = mean(mean)) |>
  ungroup() |>
  group_by(species_complete) |>
  summarise(trait_value = mean(mean_value_trait))

unique(subdata_filtered_mass$species_complete)

# hedges g to bms #
result$hedgesg <- abs(result$hedgesg)
result_mass <- result %>%
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))

### putting only species mass that is present in hedgesg table #
table_bms <- table_bms[table_bms$species_complete %in% result_mass$species_complete,]
result_mass <- result_mass[result_mass$species_complete %in% table_bms$species_complete,]
table_bms$trait_value <- round(table_bms$trait_value, 3)
hedgesg_mass <- abs(result_mass$hedgesg)

### removing species with big mass (outlier) from tables ###
species_outlier <- table_bms |> 
  arrange(desc(trait_value)) |>
  select(species_complete) |>
  slice(1)
species_outlier <- unclass(species_outlier)$species_complete

hedgesg_mass <- hedgesg_mass[names(hedgesg_mass) != species_outlier]
table_bms <- table_bms[table_bms$species_complete != species_outlier, ]
result_mass <- result_mass[result_mass$species_complete != species_outlier, ]

species_mass <- unique(result_mass$species_complete)