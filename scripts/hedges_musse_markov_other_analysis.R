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
  
  tree_time_tree_ready <- force.ultrametric(reptiles_tree_time_tree)
  hedgesg <- hedgesg_without_lack
} else {
  message("Error! Chose 'yes' or 'not' ")
}

### quasse ###
lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0,
                   single.lineage=FALSE, verbose=TRUE)
nodes <- c("nd13", "nd9", "nd5")
split.t <- Inf
pars <- c(.1, .2, 0, 2.5, .03, 0, .01)
pars4 <- unlist(rep(list(pars), 4))
sd <- 1/200
control.C.1 <- list(dt.max=1/200)
## Not run:
control.R.1 <- list(dt.max=1/200, method="fftR")
control <- list(parscale = 0.1, reltol = 0.001)
lik.C.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x, control.C.1)
init <- starting.point.quasse(phy, phy$tip.state, states.sd=NULL)
result_quasse <- find.mle(lik.C.1, x.init = init)

quasse <- make.quasse(tree_time_tree_ready, hedgesg, sd, sigmoid.x, constant.x, control.C.1)
init <- starting.point.quasse(tree_time_tree_ready, hedgesg)
result_quasse <- find.mle(quasse, x.init = init, lower = 0, control = control, verbose = 0)
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
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=2, make.null=TRUE)


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
