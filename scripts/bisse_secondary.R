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

### READING AND TREATING DATA ###
#### bisse ####
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
time <- 150000

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

hedgesg <- sapply(hedgesg, states)
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
  
  tree_time_tree_ready <- force.ultrametric(reptiles_tree_time_tree)
  hedgesg <- hedgesg_without_lack
} else {
  message("Error! Chose 'yes' or 'not' ")
}

# full #
bisse_one <- make.bisse(tree = tree_time_tree_ready, states = hedgesg)
initial <- starting.point.bisse(tree_time_tree_ready)
resu <- find.mle(bisse_one, initial)
round(resu$par, 9)

# null #
bisse_null <- constrain(bisse_one, lambda1~lambda0,mu1~mu0)
resu_null <- find.mle(bisse_null, initial[argnames(bisse_null)])
round(resu_null$par, 9)

bisseAnova<- anova(resu,
                  null=resu_null)
bisseAnova

aicw(setNames(bisseAnova$AIC,
              rownames(bisseAnova)))

save.image(paste0("output/", "phy_expanded",
                  "_", phylogeny_expanded,
                  "stat", "_",states_choice, "_", 
                  "markov", "_", seeds[1], "_", 
                  "envi.RDS"))

prior <-make.prior.exponential(1/2)

preliminar <- diversitree::mcmc(bisse_one, 
                                resu$par, 
                                nsteps=100, prior=prior,
                                w=1, print.every = 0)

w <- diff(sapply(preliminar[2:(ncol(preliminar) -1)], quantile, c(0.05, 0.95)))

bisse.mcmc <- diversitree::mcmc(bisse_one, 
                 initial[colnames(w)],
                 nsteps=time,
                 prior=prior,
                 w=w,
                 print.every=100)
write.csv2(bisse.mcmc, file=paste0("output/","phy_expanded_bisse",
                                    "_", phylogeny_expanded,
                                    "stat", "_",
                                    states_choice, "_",
                                    "markov", "_",
                                    seeds[1], "_",
                                    "mcmc.csv"),
           row.names = FALSE)

# bisse analysis #
phy_expanded_cohen_one_rep_one <- read.csv2("output/phy_expanded_bisse_yesstat_one_markov_2_mcmc.csv")
state_chosen <- phy_expanded_cohen_one_rep_one

mcmc_max <- nrow(state_chosen)
mcmc_out_burn_in <- round(nrow(state_chosen) * 0.2) + 1
mcmc_result <- state_chosen[mcmc_out_burn_in:mcmc_max, ]

n.eff(as.matrix(state_chosen[, 2:(length(mcmc_result) - 1)]))
n.eff(as.matrix(state_chosen[, 2:3]))
n.eff(as.matrix(state_chosen[, 4:5]))
n.eff(as.matrix(state_chosen[, 6:7]))

mean_posteriors <- colMeans(mcmc_result)[2:ncol(mcmc_result)]
round(mean_posteriors, 5)
