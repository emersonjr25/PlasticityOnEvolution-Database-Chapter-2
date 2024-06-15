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

### READING DATA ###

dados <- read.csv("data/raw/Database.csv", header = TRUE, sep = ",") 

################################################################################
##### MUSSE TO ALL TRAITS - NCBI Database ################################
################################################################################
#### phylogeny time tree ####
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
cat("whales rate: mean", mean(rates$lambda-rates.whales$mu),"sd", sd(rates.whales$lambda-rates.whales$mu))
cat("lamda: mean", mean(rates$lambda), "sd", sd(rates.whales$lambda))
cat("mu: mean", mean(rates$mu),"sd", sd(rates.whales$mu))
hist(rates$lambda-rates.whales$mu, main = NULL, xlab = "Diversification rate")
abline(v = mean(rates$lambda-rates.whales$mu), col = "red", lwd = 2)
marg_probs <- marginalShiftProbsTree(results)

#### bisse e hisse ####
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
    if(x <= 0.35){
      x <- 0
    } else if (x > 0.35){
      x <- 1
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
    if(x <= 0.5){
      x <- 0
    }  else if(x > 0.5){
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

bisse_one <- make.bisse(tree = tree_time_tree_ready, states = hedgesg)
initial <- starting.point.bisse(tree_time_tree_ready)
resu <- find.mle(bisse_one, initial, method = "subplex")
round(resu$par, 9)
trans.rates.hisse <- TransMatMakerHiSSE()

hisse(tree_time_tree_ready, hedgesg, hidden.states=FALSE,trans.rate=trans.rates.hisse)
hisse(tree_time_tree_ready, hedgesg)


make.hisse

library(geiger);
library(phytools);
library(phyloTop)
library(TreeSim);
library(ape);
library(diversitree)
library(apTreeshape); 
library(RPANDA)
library(BAMMtools)
library(hisse)
### muhisse ###
pars <- c(.1, .15, .2, .1, # lambda 1, 2, 3, 4
          .03, .045, .06, 0.03, # mu 1, 2, 3, 4
          .05, .05, .00, # q12, q13, q14
          .05, .00, .05, # q21, q23, q24
          .05, .00, .05, # q31, q32, q34
          .00, .05, .05)

set.seed(2)
phy <- tree.musse(pars, 30, x0=1)
states <- phy$tip.state
lik <- make.musse(phy, states, 4)
#lik <- make.musse(phy, states, 3)
diversitree.free = lik(pars)
print(diversitree.free)
states <- data.frame(phy$tip.state, phy$tip.state,
                     row.names=names(phy$tip.state))
states <- states[phy$tip.label,]
states.trans <- states
for(i in 1:Ntip(phy)){
  if(states[i,1] == 1){
    states.trans[i,1] = 0
    states.trans[i,2] = 0
  }
  if(states[i,1] == 2){
    states.trans[i,1] = 0
    states.trans[i,2] = 1
  }
  if(states[i,1] == 3){
    states.trans[i,1] = 1
    states.trans[i,2] = 0
  }
  if(states[i,1] == 4){
    states.trans[i,1] = 1
    states.trans[i,2] = 1
  }
}
pars.hisse <- c(pars[1]+pars[5],pars[2]+pars[6],pars[3]+pars[7],pars[4]+pars[8],
                pars[5]/pars[1],pars[6]/pars[2],pars[7]/pars[3],pars[8]/pars[4],
                0.05,0.05,0, 0.05,0,0.05, 0.05,0,.05, 0,0.05,.05)
model.vec = rep(0,384)
model.vec[1:20] = pars.hisse
phy$node.label = NULL
cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=TRUE,
                                         nb.tip=Ntip(phy), nb.node=Nnode(phy),
                                         bad.likelihood=exp(-300), f=c(1,1,1,1), ode.eps=0)
cache$psi <- 0
gen <- hisse:::FindGenerations(phy)
dat.tab <- hisse:::OrganizeData(states.trans, phy, f=c(1,1,1,1), hidden.states=TRUE)
hisse.constrained <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache,
                                             root.type="madfitz", condition.on.survival=TRUE,
                                             root.p=NULL)
comparison <- identical(round(hisse.constrained,4), round(diversitree.free,4))
print(comparison)
pars.hisse <- rep(c(pars[1]+pars[5],pars[2]+pars[6],pars[3]+pars[7],pars[4]+pars[8],
                    pars[5]/pars[1],pars[6]/pars[2],pars[7]/pars[3],pars[8]/pars[4],
                    0.05,0.05,0, 0.05,0,0.05, 0.05,0,.05, 0,0.05,.05, 1,rep(0,6), 1,
                    rep(0,6), 1,rep(0,6), 1,rep(0,6)),2)
model.vec = rep(0,384)
model.vec[1:96] = pars.hisse
phy$node.label = NULL
cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=TRUE,
                                         nb.tip=Ntip(phy), nb.node=Nnode(phy),
                                         bad.likelihood=exp(-300), f=c(1,1,1,1),
                                         ode.eps=0)
cache$psi <- 0
gen <- hisse:::FindGenerations(phy)
dat.tab <- hisse:::OrganizeData(states.trans, phy, f=c(1,1,1,1), hidden.states=TRUE)
hisse.constrained <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache,
                                             root.type="madfitz", condition.on.survival=TRUE,
                                             root.p=NULL)
comparison <- identical(round(hisse.constrained,4), round(diversitree.free,4))
print(comparison)
turnover <- c(1,1,1,1)
extinction.fraction <- c(1,1,1,1)
f = c(1,1,1,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
print(trans.rate)
states.trans <- cbind(phy$tip.label, states.trans)
dull.null <- MuHiSSE(phy=phy, data=states.trans, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rate)
turnover <- c(1,2,3,4)
extinction.fraction <- c(1,1,1,1)
MuSSE <- MuHiSSE(phy=phy, data=states.trans, f=f, turnover=turnover,
                 eps=extinction.fraction, hidden.states=FALSE,
                 trans.rate=trans.rate)
turnover <- c(1,2,3,4,5,6,7,8)
extinction.fraction <- rep(1, 8)
f = c(1,1,1,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1)
print(trans.rate)
MuHiSSE <- MuHiSSE(phy=phy, data=states.trans, f=f, turnover=turnover,
                   eps=extinction.fraction, hidden.states=TRUE,
                   trans.rate=trans.rate)


turnover <- c(1,1,1,1,2,2,2,2)
extinction.fraction <- rep(1, 8)
f = c(1,1,1,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1, make.null=TRUE)
print(trans.rate)
## Three hidden states
turnover <- c(rep(1,4), rep(2,4), rep(3,4))
extinction.fraction <- rep(1, 12)

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


