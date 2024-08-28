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

### READING AND TREATING DATA - BAMM ###
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
mcmcout <- read.table("mcmc_out.txt", sep=',', header = T)
event_data <- read.table("event_data.txt", sep=',', header = T)

plot(mcmcout$logLik ~ mcmcout$generation, pch = 19)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
conver <- sapply(postburn, effectiveSize)

tree_time_tree_ready$tip.label <- gsub(" ", "_", tree_time_tree_ready$tip.label)

results <- BAMMtools::getEventData(tree_time_tree_ready, event_data, burnin=0.25, nsamples=500)

plotRateThroughTime(results)
tip <- getTipRates(results)
hist(tip$lambda.avg, xlab = "Average lambda", main = NULL)
rates <- getCladeRates(results)
marg_probs <- marginalShiftProbsTree(results)

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
  # states first way - around 0.2, 0.5, 0.8 #
  states <- function(x){
    if(x <= 0.35){
      x <- 0
    } else if (x > 0.35){
      x <- 1
    } 
  }
} else {
  message("Error! Chose 'one', 'two', or 'three")
}

# states_choice <- c("one") #can be one, two, or three
# if(states_choice == "one"){
#   # states first way - around 0.2, 0.5, 0.8 #
#   states <- function(x){
#     if(x <= 0.35){
#       x <- 1
#     } else if (x > 0.35 & x <= 0.65){
#       x <- 2
#     } else if(x > 0.65){
#       x <- 3
#     }
#   }
# } else if(states_choice == "two"){
#   # states second way - using 3 quartiles #
#   states <- function(x){
#     if(x <= first_quartile){
#       x <- 1
#     } else if (x > first_quartile & x <= second_quartile){
#       x <- 2
#     } else if(x > second_quartile){
#       x <- 3
#     }
#   }
# } else if(states_choice == "three"){
#   # states third way - some articles #
#   # can considered very low, low, and medium/high or #
#   # low, medium, and high #
#   states <- function(x){
#     if(x <= 0.2){
#       x <- 1
#     } else if (x > 0.2 & x <= 0.5){
#       x <- 2
#     } else if(x > 0.5){
#       x <- 3
#     }
#   }
# } else {
#   message("Error! Chose 'one', 'two', or 'three")
# }

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

plot(result_final_bamm$hedgesg, result_final_bamm$tip.lambda.avg)
plot(result_final_bamm$hedgesg, result_final_bamm$tip.mu.avg)

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(mean(tip.lambda.avg))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(mean(tip.mu.avg))

result_final_bamm %>% 
  ggplot(aes(tip.lambda.avg, fill = as.factor(hedgesg))) + 
  geom_density(alpha=0.7)

result_final_bamm %>% 
  ggplot(aes(tip.mu.avg, fill = as.factor(hedgesg))) + 
  geom_density(alpha=0.7)
