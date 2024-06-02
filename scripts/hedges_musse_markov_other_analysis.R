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

results <- BAMMtools::getEventData(tree_time_tree_ready, event_data, burnin=0.25, nsamples=500)

### with whales ###
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)

bamm.whales <- plot.bammdata(ed, lwd = 2, method = "phylogram", labels = TRUE, cex = 0.5);
y <- addBAMMshifts(ed, cex = 2);
addBAMMlegend(bamm.whales, direction = "vertical", location = "right", nTicks = 4, side = 2, las = 1, cex.axis = 1);
axisPhylo(side = 1, backward = T,  las = 1, cex = 1.5);
mtext("Millions of years before present", side = 1, line = 2.5, cex = 1.5)
plotRateThroughTime(ed, intervalCol = "red", avgCol = "red", ylim = c(0, 1), cex.axis = 1.5, ratetype = "speciation");
text(x = 30, y = 0.8, label = "All whales", font = 4, cex = 2.0, pos = 4)

tip.whales <- getTipRates(ed)
# Agora, explore as taxas de especiação
hist(tip.whales$lambda.avg, xlab = "Average lambda", main = NULL)
dolphins <- subtreeBAMM(ed, node = 141)
tip.dolphins <- getTipRates(dolphins)
# Explore as taxas de especiação dos golfinhos
hist(tip.dolphins$lambda.avg, xlab = "Average lambda", main = NULL)
#VocÊ pode também explorar os plots prévios de diversificação ao longo do tem
# Pegue as taxas por clado para todas as baleias
rates.whales <- getCladeRates(ed)
# Calcule as taxas para as baleias:
# Diversificação
cat("whales rate: mean", mean(rates.whales$lambda-rates.whales$mu),"sd", sd(rates.whales$lambda-rates.whales$mu))
# Especiação
cat("lamda: mean", mean(rates.whales$lambda), "sd", sd(rates.whales$lambda))
# Extinção
cat("mu: mean", mean(rates.whales$mu),"sd", sd(rates.whales$mu))
# Ou, apenas plote as taxas 
hist(rates.whales$lambda-rates.whales$mu, main = NULL, xlab = "Diversification rate")
abline(v = mean(rates.whales$lambda-rates.whales$mu), col = "red", lwd = 2)


####
data(whales, package = "geiger")
w.phy <- whales$phy
ltt.w <- ltt(w.phy,log.lineages=F)
t = max(branching.times(w.phy))
#Quantidade de linhagens de baleia no tempo inicial ("stem age"; n0):
n0 = 1
#Quantidade atual de linhagens de baleia (nt):
nt = 84
#Funcão exponencial para estimar diversificacão (l1 [lambda]; observe que nós logaritmizamos a funcão e isolamos o parâmetro lambda):
l1 = (log(nt) - log(n0))/t
l2 = seq(0.01,0.25,0.002)
#Função para gerar a verossimilhança para cada valor de lambda:
lik <- (((2.72^(l2*t))-1)^(nt-1)) / ((2.72^(l2*t))^(nt))
f.lamb1 <- function(t, y){y[1]}#especiacão
f.mu1 <- function(t, y){0}#extincão
#Determinando os valores iniciais dos parâmetros para facilitar a busca de máxima verossimilhança:
lamb_par1 <- c(0.09)#especiacão
mu_par1 <- c()#extincão
resu1 = fit_bd(phylo = w.phy, tot_time = t, f.lamb1, f.mu1, lamb_par1, mu_par1)
resu1$lamb_par#Lambda de máxima veroximilhanca
f.lamb2 <- function(t, y){y[1]}
f.mu2 <- function(t, y){y[1]}
#Determinando os valores iniciais dos parâmetros para facilitar a busca de máxima verossimilhanca:
lamb_par2 <- c(0.09)
mu_par2 <- c(0.005)
resu2 <- fit_bd(phylo = w.phy, tot_time = t, f.lamb2, f.mu2, lamb_par2, mu_par2)
resu2$lamb_par#Lambda
resu2$mu_par#Mu

cont.tr <- whales$dat[,1]
#Como algumas espécies não possuem valores de tamanho corporal, nós geramos valores de forma arbitrária. 
#Identificando quais espécies na filogenia não possuem o atributo
mis.names <- w.phy$tip.label[which(!(w.phy$tip.label%in%rownames(cont.tr)))]
#Identificando a média da distribuição de tamanho corporal
mean.val <- mean(cont.tr)
#Amostrando os valores a partir de uma distribuição normal
mis.val <- rnorm(length(which(!(w.phy$tip.label%in%rownames(cont.tr)))), mean = mean.val )
#Atribuindo os valores amostrados para as espécies
names(mis.val) <- mis.names
#Concatenando todas as espécies
cont.tr2 <- c(cont.tr,mis.val)
#Binarizando a variável de tamanho corporal para "grande" (estado 1) e "pequeno" (estado 0):
bin.tr <- ifelse(cont.tr2 < mean.val, 0, 1)
#Determinado os valores iniciais dos parâmetros para facilitar a busca de máxima verossimilhanca:
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
#Gerando a função de máxima verossimilhanca
llik.fun <- make.bisse(tree = w.phy, states = bin.tr)
#Fazendo a busca de máxima verossimilhanca:
resu3 <- find.mle(llik.fun, pars, method = "subplex")
resu3$par


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

trans.rates.hisse <- TransMatMakerHiSSE(hidden.traits=0)
hisse(tree_time_tree_ready, hedgesg, hidden.states=TRUE, turnover=c(1,2,1,2),
      eps=c(1,2,1,2), trans.rate=trans.rates.hisse)
hisse(tree_time_tree_ready, hedgesg)
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
