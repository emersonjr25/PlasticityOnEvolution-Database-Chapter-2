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
library(phytools)
library(ggplot2)
library(geiger)
library(randtip)
library(stableGR)
library(OUwie)
library(corHMM)
library(bayou)
library(ggtree)
library(dispRity)
library(BAMMtools)
library(hisse)

###################################################################
#### FIRST WAY ####
### READING AND TREATING DATA - BAMM ###
###################################################################
load("table_and_phy_ready.RDS")

result$hedgesg <- abs(result$hedgesg)
result_all_species <- result %>%
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))
hedgesg <- abs(result_all_species$hedgesg)

reptiles_tree_time_tree <- read.newick("data/raw/species.nwk")

hedgesg <- setNames(hedgesg, result_all_species$species_complete)

### choose phylogeny - expanded or not ###
phylogeny_expanded <- "yes" # chose yes or not
seed_phy <- c(100, 101)
set.seed(seed_phy[1])
if(phylogeny_expanded == "yes"){
  reptiles_tree_time_tree$tip.label <- gsub("_", " ", reptiles_tree_time_tree$tip.label)
  info_fix_poly <- build_info(names(hedgesg), reptiles_tree_time_tree, 
                              find.ranks=TRUE, db="ncbi")
  input_fix_poly <- info2input(info_fix_poly, reptiles_tree_time_tree,
                               parallelize=F)
  tree_time_tree_ready <- rand_tip(input = input_fix_poly, 
                                   tree = reptiles_tree_time_tree,
                                   forceultrametric=TRUE,
                                   prune=TRUE)
  tree_time_tree_ready$tip.label <- gsub("_", " ", tree_time_tree_ready$tip.label)
} else if (phylogeny_expanded == "not"){
  ### manual corrections in phylogeny to use phylogeny directly ###
  lack_species <- c("Anepischetosia maccoyi", "Nannoscincus maccoyi")
  reptiles_tree_time_tree$tip.label <- gsub("_", " ", reptiles_tree_time_tree$tip.label)
  hedgesg_without_lack <- hedgesg[!names(hedgesg) %in% lack_species]
  different_species <- reptiles_tree_time_tree$tip.label[!reptiles_tree_time_tree$tip.label %in% names(hedgesg_without_lack)]
  different_species_hedgesg <- names(hedgesg_without_lack)[!names(hedgesg_without_lack) %in%  reptiles_tree_time_tree$tip.label]
  
  names(hedgesg_without_lack)[grepl('ruf', names(hedgesg_without_lack))] <- different_species[different_species == "Lycodon rufozonatus"]
  reptiles_tree_time_tree$tip.label[grepl('Elaphe tae', reptiles_tree_time_tree$tip.label)] <- different_species_hedgesg[different_species_hedgesg == "Elaphe taeniura"]
  names(hedgesg_without_lack)[grepl('sinensis', names(hedgesg_without_lack))][1] <- different_species[different_species == "Mauremys sinensis"]
  names(hedgesg_without_lack)[grepl('piscator', names(hedgesg_without_lack))] <- different_species[different_species == "Fowlea piscator"]
  names(hedgesg_without_lack)[grepl('lesueurii', names(hedgesg_without_lack))][1] <- different_species[different_species == "Amalosia lesueurii"]
  names(hedgesg_without_lack)[grepl('lesueurii', names(hedgesg_without_lack))][2] <- different_species[different_species == "Intellagama lesueurii"]
  
 # tree_time_tree_ready <- force.ultrametric(reptiles_tree_time_tree)
  hedgesg <- hedgesg_without_lack
} else {
  message("Error! Chose 'yes' or 'not' ")
}

### BAMM ###
tree_time_tree_ready <- dispRity::remove.zero.brlen(tree_time_tree_ready)
#write.tree(tree_time_tree_ready, "fil.tre")
mcmcout <- read.table("mcmc_out.txt", sep=',', fill=T, header = T)
event_data <- read.table("event_data.txt", sep=',', header = T)

plot(mcmcout$logLik ~ mcmcout$generation, pch = 19)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
conver <- sapply(postburn, effectiveSize)

tree_time_tree_ready$tip.label <- gsub(" ", "_", tree_time_tree_ready$tip.label)

results <- BAMMtools::getEventData(tree_time_tree_ready, event_data, burnin=0.25)

#plotRateThroughTime(results)
tip <- getTipRates(results)
#hist(tip$lambda.avg, xlab = "Average lambda", main = NULL)
#rates <- getCladeRates(results)
#marg_probs <- marginalShiftProbsTree(results)

### rates ###
lambda <- data.frame(tip$lambda.avg)
mu <- data.frame(tip$mu.avg)

### classifying ###
load("table_and_phy_ready.RDS")

result$hedgesg <- abs(result$hedgesg)
result_all_species <- result %>%
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))
hedgesg <- abs(result_all_species$hedgesg)

reptiles_tree_time_tree <- read.newick("data/raw/species.nwk")

### STATES ###
first_quartile <- round(summary(hedgesg)[2], 2)
second_quartile <- round(summary(hedgesg)[3], 2)

seeds <- c(2, 3, 4)
time <- 100000

states_choice <- c("one") #can be one, two, or three
if(states_choice == "one"){
  # states first way - around 0.2, 0.5, 0.8 #
  states <- function(x){
    if(x <= 0.65){
      x <- 0
    } else if (x > 0.65){
      x <- 1
    } 
  }
} else if(states_choice == "two"){
  # states second way - using 3 quartiles #
  states <- function(x){
    if(x <= second_quartile){
      x <- 0
    } else if(x > second_quartile){
      x <- 1
    }
  }
} else if(states_choice == "three"){
  # states third way - some articles #
  # can considered very low, low, and medium/high or #
  # low, medium, and high #
  states <- function(x){
    if(x <= 0.5){
      x <- 0
    }  else if(x > 0.5){
      x <- 1
    }
  } 
} else if(states_choice == "fourth"){
  # states third way - some articles #
  # can considered very low, low, and medium/high or #
  # low, medium, and high #
  states <- function(x){
    if(x <= 0.35){
      x <- 0
    }  else if(x > 0.35){
      x <- 1
    }
  }
} else {
  message("Error! Chose 'one', 'two','three' or fourth" )
}

hedgesg <- sapply(hedgesg, states)
hedgesg <- setNames(hedgesg, result_all_species$species_complete)

states <- data.frame(hedgesg)
states <- rownames_to_column(states, var='species')
states$species <- tolower(states$species)

mu <- rownames_to_column(mu, var='species')
mu$species <- gsub("_", " ", mu$species)
mu$species <- tolower(mu$species)

lambda <- rownames_to_column(lambda, var='species')
lambda$species <- gsub("_", " ", lambda$species)
lambda$species <- tolower(lambda$species)

result_final_bamm <- inner_join(states, lambda, by='species')
result_final_bamm <- inner_join(result_final_bamm, mu, by='species')
result_final_bamm$tip.diversification <- abs(result_final_bamm$tip.lambda.avg - result_final_bamm$tip.mu.avg)

plot(result_final_bamm$hedgesg, result_final_bamm$tip.lambda.avg)
plot(result_final_bamm$hedgesg, result_final_bamm$tip.mu.avg)

teste <- result_final_bamm |>
  pivot_wider(names_from=hedgesg, values_from = starts_with("tip.lambda"))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(mean(tip.lambda.avg))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(round(mean(tip.mu.avg), 9))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(mean(tip.diversification))

result_final_bamm$hedgesg1 <- result_final_bamm$hedgesg

result_final_bamm %>% 
  ggplot(aes(x = as.factor(hedgesg), y=tip.lambda.avg)) + 
  geom_boxplot()

result_final_bamm %>% 
  ggplot(aes(x = as.factor(hedgesg), y=tip.mu.avg)) + 
  geom_boxplot()

result_final_bamm %>% 
  ggplot(aes(x = as.factor(hedgesg), y=tip.diversification)) + 
  geom_boxplot()


###################################################################
#### SECOND WAY ####
### READING AND TREATING DATA - BAMM ###
###################################################################
tree_time_tree_ready <- read.tree('filogenia_repteis.tre')
tree_time_tree_ready <- read.tree('data/raw/repteis_especies_total.nwk')
species_total <- read.table('data/raw/repteis_especies_total.txt', sep=',')
tree_time_tree_ready$tip.label <- gsub("_", " ", tree_time_tree_ready$tip.label)
info_fix_poly <- build_info(unlist(species_total$V1), tree_time_tree_ready, 
                            find.ranks=TRUE, db="gbif")
#saveRDS(info_fix_poly, 'info_fix_poly.rds')
#info_fix_poly <- read_rds("info_fix_poly.rds")
input_fix_poly <- info2input(info_fix_poly, tree_time_tree_ready,
                             parallelize=F)
#saveRDS(input_fix_poly, 'input_fix_poly.rds')
#input_fix_poly <- read_rds("input_fix_poly.rds")
tree_time_tree_ready <- rand_tip(input = input_fix_poly, 
                                 tree = tree_time_tree_ready,
                                 forceultrametric=FALSE,
                                 prune=TRUE)
#saveRDS(tree_time_tree_ready, 'tree_time_tree_ready_y_prune.rds')
#tree_time_tree_ready <- read_rds("tree_time_tree_ready_y_prune.rds")
#tree_time_tree_ready <- force.ultrametric(tree_time_tree_ready, method="extend")
#tree_time_tree_ready <- dispRity::remove.zero.brlen(tree_time_tree_ready)
#is.ultrametric(tree_time_tree_ready)
#is.binary.tree(tree_time_tree_ready)
#min(tree_time_tree_ready$edge.length)
#write.tree(tree_time_tree_ready, "filogenia_repteis_prune.tre")
#write.tree(tree_time_tree_ready, "filogenia_repteis_prune_from_orig_phy.tre")

tree_time_tree_ready <- read.tree('filogenia_repteis_prune.tre')
tree_time_tree_ready2 <- read.tree('filogenia_repteis_prune_from_orig_phy.tre')

tree_time_tree_ready <- read.tree("data/raw/squamata_pyron_2013.txt") #squamatas
tree_time_tree_ready <- force.ultrametric(tree_time_tree_ready, method="extend")
is.ultrametric(tree_time_tree_ready)
is.binary.tree(tree_time_tree_ready)
min(tree_time_tree_ready$edge.length)
write.tree(tree_time_tree_ready, "filogenia_squamata.tre")

tree_time_tree_ready <- read.newick("data/raw/repteis_especies_total.nwk") #repteis
tree_time_tree_ready <- force.ultrametric(tree_time_tree_ready, method="extend")
tree_time_tree_ready <- dispRity::remove.zero.brlen(tree_time_tree_ready)
is.ultrametric(tree_time_tree_ready)
is.binary.tree(tree_time_tree_ready)
min(tree_time_tree_ready$edge.length)
write.tree(tree_time_tree_ready, "filogenia_repteis.tre")

##### reading phylogenny ######
tree_time_tree_ready <- read.tree('filogenia_repteis_prune_from_orig_phy.tre')
tree_time_tree_ready <- read.tree('filogenia_repteis_prune.tre')
mcmcout <- read.table("mcmc_out_12000.txt", sep=',', header = T)
event_data <- read.table("event_data_12000.txt", sep=',', header = T)

tree_time_tree_ready <- read.tree('filogenia_repteis.tre')
mcmcout <- read.table("mcmc_out_reptiles_2milion.txt", sep=',', header = T)
event_data <- read.table("event_data_reptiles_2milion.txt", sep=',', header = T)

tree_time_tree_ready <- read.tree('filogenia_repteis.tre')
mcmcout <- read.table("mcmc_out_reptiles.txt", sep=',', header = T)
event_data <- read.table("event_data_reptiles.txt", sep=',', header = T)

tree_time_tree_ready <- read.tree('filogenia_squamata.tre')
mcmcout <- read.table("mcmc_out_squamata.txt", sep=',', header = T)
event_data <- read.table("event_data_squamata.txt", sep=',', header = T)

plot(mcmcout$logLik ~ mcmcout$generation, pch = 19)
burnstart <- floor(0.25 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
conver <- sapply(postburn, effectiveSize)

tree_time_tree_ready$tip.label <- gsub(" ", "_", tree_time_tree_ready$tip.label)

results <- BAMMtools::getEventData(tree_time_tree_ready, 
                                   event_data, burnin=0.25,
                                   nsamples=500)

#plotRateThroughTime(results)
tip <- getTipRates(results)
#hist(tip$lambda.avg, xlab = "Average lambda", main = NULL)
#rates <- getCladeRates(results)
#marg_probs <- marginalShiftProbsTree(results)

### rates ###
lambda <- data.frame(tip$lambda.avg)
mu <- data.frame(tip$mu.avg)

### classifying ###
load("table_and_phy_ready.RDS")

result$hedgesg <- abs(result$hedgesg)
result_all_species <- result %>%
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))
hedgesg <- abs(result_all_species$hedgesg)

reptiles_tree_time_tree <- read.newick("data/raw/species.nwk")

### STATES ###
first_quartile <- round(summary(hedgesg)[2], 2)
second_quartile <- round(summary(hedgesg)[3], 2)

seeds <- c(2, 3, 4)
time <- 100000

states_choice <- c("one") #can be one, two, or three
if(states_choice == "one"){
  # states first way - around 0.2, 0.5, 0.8 #
  states <- function(x){
    if(x <= 0.65){
      x <- 0
    } else if (x > 0.65){
      x <- 1
    } 
  }
} else if(states_choice == "two"){
  # states second way - using 3 quartiles #
  states <- function(x){
    if(x <= second_quartile){
      x <- 0
    } else if(x > second_quartile){
      x <- 1
    }
  }
} else if(states_choice == "three"){
  # states third way - some articles #
  # can considered very low, low, and medium/high or #
  # low, medium, and high #
  states <- function(x){
    if(x <= 0.5){
      x <- 0
    }  else if(x > 0.5){
      x <- 1
    }
  } 
} else if(states_choice == "fourth"){
  # states third way - some articles #
  # can considered very low, low, and medium/high or #
  # low, medium, and high #
  states <- function(x){
    if(x <= 0.35){
      x <- 0
    }  else if(x > 0.35){
      x <- 1
    }
  }
} else {
  message("Error! Chose 'one', 'two','three' or fourth" )
}

hedgesg <- sapply(hedgesg, states)
hedgesg <- setNames(hedgesg, result_all_species$species_complete)

states <- data.frame(hedgesg)
states <- rownames_to_column(states, var='species')
states$species <- tolower(states$species)

mu <- rownames_to_column(mu, var='species')
mu$species <- gsub("_", " ", mu$species)
mu$species <- tolower(mu$species)

lambda <- rownames_to_column(lambda, var='species')
lambda$species <- gsub("_", " ", lambda$species)
lambda$species <- tolower(lambda$species)

result_final_bamm <- left_join(states, lambda, by='species')
result_final_bamm <- left_join(result_final_bamm, mu, by='species')
result_final_bamm$tip.diversification <- abs(result_final_bamm$tip.lambda.avg - result_final_bamm$tip.mu.avg)

result_final_bamm <- result_final_bamm[!is.na(result_final_bamm[, "tip.lambda.avg"]),]

plot(result_final_bamm$hedgesg, result_final_bamm$tip.lambda.avg)
plot(result_final_bamm$hedgesg, result_final_bamm$tip.mu.avg)

#teste <- result_final_bamm |>
#  pivot_wider(names_from=hedgesg, values_from = starts_with("tip.lambda"))

#result_final_bamm <- result_final_bamm[(result_final_bamm$tip.lambda.avg >= summary(result_final_bamm$tip.lambda.avg)[2] & 
#                    result_final_bamm$tip.lambda.avg <= summary(result_final_bamm$tip.lambda.avg)[5]), ] 

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(mean(tip.lambda.avg))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(round(mean(tip.mu.avg), 9))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(mean(tip.diversification))

result_final_bamm$hedgesg1 <- result_final_bamm$hedgesg

result_final_bamm %>% 
  ggplot(aes(x = as.factor(hedgesg), y=tip.lambda.avg)) + 
  geom_boxplot()

result_final_bamm %>% 
  ggplot(aes(x = as.factor(hedgesg), y=tip.mu.avg)) + 
  geom_boxplot()

result_final_bamm %>% 
  ggplot(aes(x = as.factor(hedgesg), y=tip.diversification)) + 
  geom_boxplot()

