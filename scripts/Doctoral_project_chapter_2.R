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
library(rotl)
library(ape)
library(taxize)
library(diversitree)
library(RRphylo)

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

################################################################################
##### MUSSE TO ALL TRAITS - NCBI Database ################################
################################################################################

#### Preparation to MUSSE ####
result_all_species <- result %>% 
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))

### STATES ###
hedgesg <- abs(result_all_species$hedgesg)

states <- function(x){
  if(x <= 0.2){
    x <- 1
  } else if (x > 0.2 & x <= 0.5){
    x <- 2
  } else if(x > 0.5 & x <= 0.8){
    x <- 3
  } else if(x > 0.8){
    x <- 4
  }
}

hedgesg <- sapply(hedgesg, states)
hedgesg <- setNames(hedgesg, result_all_species$species_complete)

species <- unique(result_all_species$species_complete)

#### phylogeny construction - NCBI ####
species_info_ncbi <- classification(species, db='ncbi')
species_tree_ncbi <- class2tree(species_info_ncbi, check = TRUE)
species_tree_ncbi <- species_tree_ncbi$phylo

### MUSSE CALCULATION NCBI ###
resolved_tree_ncbi <- multi2di(species_tree_ncbi)
species_wrong <- names(hedgesg)[!names(hedgesg) %in% resolved_tree_ncbi$tip.label]
species_wrong2 <- resolved_tree_ncbi$tip.label[!resolved_tree_ncbi$tip.label %in% names(hedgesg)]

### manual corrections in phylogeny ###
names(hedgesg)[grepl('ruf', names(hedgesg))] <- species_wrong2[species_wrong2 == "Lycodon rufozonatus"]
resolved_tree_ncbi$tip.label[grepl('Elaphe tae', resolved_tree_ncbi$tip.label)] <- species_wrong[species_wrong == "Elaphe taeniura"]
names(hedgesg)[grepl('maccoyi', names(hedgesg))][1] <- species_wrong2[species_wrong2 == "Anepischetosia maccoyi.x"]
names(hedgesg)[grepl('maccoyi', names(hedgesg))][2] <- species_wrong2[species_wrong2 == "Anepischetosia maccoyi.y"]
names(hedgesg)[grepl('sinensis', names(hedgesg))][1] <- species_wrong2[species_wrong2 == "Mauremys sinensis"]
names(hedgesg)[grepl('piscator', names(hedgesg))] <- species_wrong2[species_wrong2 == "Fowlea piscator"]
names(hedgesg)[grepl('lesueurii', names(hedgesg))][1] <- species_wrong2[species_wrong2 == "Amalosia lesueurii"]
names(hedgesg)[grepl('lesueurii', names(hedgesg))][2] <- species_wrong2[species_wrong2 == "Intellagama lesueurii"]

#### musse ####
resolved_tree_ncbi <- fix.poly(resolved_tree_ncbi, type='resolve')
musse <- make.musse(resolved_tree_ncbi, states = hedgesg, k = 4)
init <- starting.point.musse(resolved_tree_ncbi, k = 4)
result_musse <- find.mle(musse, x.init = init)
round(result_musse$par, 9)

# musse constrain #
musse_null <- constrain(musse, lambda2 ~ 0, lambda3 ~ 0,
                        lambda4 ~ 0, mu2 ~ 0, mu3 ~ 0,
                        mu4 ~ 0, q13 ~ 0, q31 ~ 0, 
                        q41 ~ 0, q23 ~ 0, q24 ~ 0)
result_musse_null <- find.mle(musse_null, x.init = init[argnames(musse_null)])

################################################################################
##### MUSSE TO MORE FREQUENT TRAIT - NCBI  Database #######################
################################################################################

### preparation MUSSE ###
trait_frequent <- result %>%
  group_by(simp_trait) %>%
  summarise(len = length(simp_trait)) %>%
  arrange(desc(len))

result_mass <- result %>%
  group_by(species_complete) %>%
  filter(simp_trait == 'Mass') %>%
  summarise(hedgesg = mean(hedgesg))

hedgesg_mass <- abs(result_mass$hedgesg)

states <- function(x){
  if(x <= 0.2){
    x <- 1
  } else if (x > 0.2 & x <= 0.5){
    x <- 2
  } else if(x > 0.5 & x <= 0.8){
    x <- 3
  } else if(x > 0.8){
    x <- 4
  }
}

hedgesg_mass <- sapply(hedgesg_mass, states)
hedgesg_mass <- setNames(hedgesg_mass, result_mass$species_complete)

species_mass <- unique(result_mass$species_complete)

### phylogenetic construction mass - NCBI ###
species_info_ncbi_mass <- classification(species_mass, db='ncbi')
species_info_ncbi_mass <- class2tree(species_info_ncbi_mass)
species_tree_ncbi_mass <- species_info_ncbi_mass$phylo

### MUSSE CALCULATION NCBI ###
resolved_tree_ncbi <- multi2di(species_tree_ncbi_mass)
species_wrong <- names(hedgesg_mass)[!names(hedgesg_mass) %in% resolved_tree_ncbi$tip.label]
species_wrong2 <- resolved_tree_ncbi$tip.label[!resolved_tree_ncbi$tip.label %in% names(hedgesg_mass)]

### manual corrections in phylogeny ###
names(hedgesg_mass)[grepl('piscator', names(hedgesg_mass))] <- species_wrong2[species_wrong2 == "Fowlea piscator"]
names(hedgesg_mass)[grepl('lesueurii', names(hedgesg_mass))][1] <- species_wrong2[species_wrong2 == "Amalosia lesueurii"]

#### musse ####
resolved_tree_ncbi <- fix.poly(resolved_tree_ncbi, type='resolve')
musse <- make.musse(resolved_tree_ncbi, states = hedgesg_mass, k = 4)
init <- starting.point.musse(resolved_tree_ncbi, k = 4)
result_musse <- find.mle(musse, x.init = init)
round(result_musse$par, 7)

plot(species_tree_ncbi_mass)

resolved_tree_ncbi_mass <- multi2di(species_tree_ncbi_mass)
resolved_tree_ncbi_mass$tip.label <- result_mass$species_complete
musse <- make.musse(resolved_tree_ncbi_mass, states = hedgesg_mass, k = 4)
init <- starting.point.musse(resolved_tree_ncbi_mass, k = 4)
result_musse <- find.mle(musse, x.init = init)

########################################################
#### phylogeny construction - OTT ####
#species_info_ott <- tnrs_match_names(species)
#species_tree_ott <- tol_induced_subtree(ott_ids = species_info_ott$ott_id)
# 
#plot(species_tree_ott)
# 
# #### MUSSE CALCULATION OTT - without edge.length - not ultrametric ###
#resolved_tree_ott <- multi2di(species_tree_ott)
#resolved_tree_ott$tip.label <- result_all_species$species_complete
# musse <- make.musse(resolved_tree_ott, states = hedgesg, k = 4)
# init <- starting.point.musse(resolved_tree, k = 4)
# result_musse <- find.mle(musse, x.init = init)

### phylogenetic construction mass - OTT ###
# species_info_ott_mass <- tnrs_match_names(species_mass)
# species_tree_ott_mass <- tol_induced_subtree(ott_ids = species_info_ott_mass$ott_id)
# plot(species_tree_ott_mass)
# 
# resolved_tree_ott_mass <- multi2di(species_tree_ott_mass)
# resolved_tree_ott_mass$tip.label <- result_mass$species_complete
# musse <- make.musse(resolved_tree_ott_mass, states = hedgesg_mass, k = 4)
# init <- starting.point.musse(resolved_tree_ott_mass, k = 4)
# result_musse <- find.mle(musse, x.init = init)
