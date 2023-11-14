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

papers_to_revise <- sample(unique(result$paper_no), 8)
sort(papers_to_revise)

################################################################################
##### MUSSE TO ALL TRAITS - NCBI Database ################################
################################################################################

#### Preparation to MUSSE ####
result_all_species <- result %>% 
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))

hedgesg <- abs(result_all_species$hedgesg)
#save.image("table_and_phy_ready.RDS")

#### phylogeny time tree ####
load("table_and_phy_ready.RDS")
reptiles_tree_time_tree <- read.newick("data/raw/species.nwk")

### STATES ###
first_quartile <- round(summary(hedgesg)[2], 2)
second_quartile <- round(summary(hedgesg)[3], 2)

seeds_phylogeny_rep <- c(100, 101, 102)
seeds <- c(2, 3, 4)
time <- 100000

states_choice <- c("one") #can be one, two, or three
if(states_choice == "one"){
  # states first way - around 0.2, 0.5, 0.8 #
  states <- function(x){
    if(x <= 0.35){
      x <- 1
    } else if (x > 0.35 & x <= 0.65){
      x <- 2
    } else if(x > 0.65){
      x <- 3
    }
  }
} else if(states_choice == "two"){
  # states second way - using 3 quartiles #
  states <- function(x){
    if(x <= first_quartile){
      x <- 1
    } else if (x > first_quartile & x <= second_quartile){
      x <- 2
    } else if(x > second_quartile){
      x <- 3
    }
  }
} else if(states_choice == "three"){
  # states third way - some articles #
  # can considered very low, low, and medium/high or #
  # low, medium, and high #
  states <- function(x){
    if(x <= 0.2){
      x <- 1
    } else if (x > 0.2 & x <= 0.5){
      x <- 2
    } else if(x > 0.5){
      x <- 3
    }
  }
} else {
  message("Error! Chose 'one', 'two', or 'three")
}

hedgesg <- sapply(hedgesg, states)
hedgesg <- setNames(hedgesg, result_all_species$species_complete)

### choose phylogeny - expanded or not ###
phylogeny_expanded <- "yes" # chose yes or not
if (phylogeny_expanded == "yes"){
  reptiles_tree_time_tree$tip.label <- gsub("_", " ", reptiles_tree_time_tree$tip.label)
  info_fix_poly <- build_info(species, reptiles_tree_time_tree, 
                              find.ranks=TRUE, db="ncbi")
  input_fix_poly <- info2input(info_fix_poly, reptiles_tree_time_tree,
                               parallelize=F)
  tree_time_tree_ready <- rand_tip(input = input_fix_poly, 
                                   tree = reptiles_tree_time_tree,
                                   forceultrametric=TRUE,
                                   prune=TRUE)
} else if (phylogeny_expanded == "not"){
  ### manual corrections in phylogeny to use phylogeny directly ###
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

for(i in seq_along(seeds)){
  set.seed(seeds[i])
  ### MUSSE CALCULATION NCBI ###
  
  ### musse ###
  # 1 musse full #
  musse_full <- make.musse(tree_time_tree_ready, states = hedgesg, k = 3)
  init <- starting.point.musse(tree_time_tree_ready, k = 3)
  result_musse_full <- find.mle(musse_full, x.init = init)
  round(result_musse_full$par, 9)
  
  # 2 musse null #
  musse_null <- constrain(musse_full, lambda2 ~ lambda1, lambda3 ~ lambda1,
                          mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q21 ~ q12, 
                          q23 ~ q12, q31 ~ 0, q32 ~ q12)
  result_musse_null <- find.mle(musse_null, x.init = init[argnames(musse_null)])
  round(result_musse_null$par, 9)
  
  # 3 musse lambda #
  musse_lambda <- constrain(musse_full, mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q21 ~ q12, 
                            q23 ~ q12, q31 ~ 0, q32 ~ q12)
  result_lambda <- find.mle(musse_lambda, x.init = init[argnames(musse_lambda)])
  round(result_lambda$par, 9)
  
  # 4 musse mu #
  musse_mu <- constrain(musse_full, lambda2 ~ lambda1, lambda3 ~ lambda1, 
                        q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
  result_mu <- find.mle(musse_mu, x.init = init[argnames(musse_mu)])
  round(result_mu$par, 9)
  
  # 5 musse ordered #
  musse_ordered <- constrain(musse_full, q13 ~ 0, q31 ~ 0)
  result_musse_ordered <- find.mle(musse_ordered, x.init = init[argnames(musse_ordered)])
  round(result_musse_ordered$par, 9)
  
  # anova to see best musse model #
  anova_result <- anova(result_musse_null,
                        all.different = result_musse_full,
                        free.lambda = result_lambda,
                        free.mu = result_mu)
  
  # look results to best model #
  aicw(setNames(anova_result$AIC, row.names(anova_result)))
  round(coef(result_musse_full), 9)
  logLik(result_musse_full)
  AIC(result_musse_full)
  round(result_musse_full$par, 9)
  
  # percentage to rates #
  prop_lambda <- round(((result_musse_full$par[1:3]  / max(result_musse_full$par[1:3])) * 100), 3)
  prop_mu <- round(((result_musse_full$par[4:6]  / max(result_musse_full$par[4:6])) * 100), 3)
  prop_transi <- round(((result_musse_full$par[7:length(result_musse_full$par)]  / max(result_musse_full$par[7:length(result_musse_full$par)])) * 100), 3)
  
  save.image(paste0("output/", "phy_expanded",
                    "_", phylogeny_expanded,
                    "stat", "_",states_choice, "_", 
                    "markov", "_", seeds[i], "_", 
                    "envi.RDS"))
  
  ###### Bayesian MCMC to find posterior density #######
  prior <- make.prior.exponential(1/2)
  
  preliminar <- diversitree::mcmc(musse_full, 
                                  result_musse_full$par, 
                                  nsteps=100, prior=prior,
                                  w=1, print.every = 0)
  
  w <- diff(sapply(preliminar[2:(ncol(preliminar) -1)], quantile, c(0.05, 0.95)))
  
  mcmc_result <- diversitree::mcmc(musse_full, 
                                   init[colnames(w)], 
                                   nsteps=time,prior=prior, 
                                   w=w, print.every=100)
  
  write.csv2(mcmc_result, file=paste0("output/","phy_expanded",
                                      "_", phylogeny_expanded,
                                      "stat", "_",
                                      states_choice, "_",
                                      "markov", "_",
                                      seeds[i], "_",
                                      "mcmc.csv"),
             row.names = FALSE)
  #save.image("mcmc.rds")
}

########### BMS CALCULATION - TRAIT EVOLUTION ##########
X_to_BMS <- abs(result_all_species$hedgesg)
X_to_BMS <- setNames(X_to_BMS, result_all_species$species_complete)
Trait <- data.frame(Genus_species = names(hedgesg),
           Reg = as.numeric(hedgesg),
           X = as.numeric(X_to_BMS))
tree_to_bms <- rayDISC(resolved_tree_ncbi, 
                       Trait[,c(1,2)], model="ER",
                       node.states="marginal")

#my_map <- make.simmap(tree_to_bms$phy,hedgesg, model="ER")
#plot(my_map)
#plotRECON(tree_to_bms$phy, tree_to_bms$states)

bms <- OUwie(tree_to_bms$phy,Trait,model=c("BMS"))
OUM <- OUwie(tree_to_bms$phy,Trait,model=c("OUM"))
aicc <- c(bms$AICc, OUM$AICc)
names(aicc) <- c("BMS", "OUM")
aic.w(aicc)
