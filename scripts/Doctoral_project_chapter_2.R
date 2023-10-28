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

#### phylogeny construction - NCBI ####
species <- unique(result_all_species$species_complete)
species_info_ncbi <- classification(species, db='ncbi')
species_tree_ncbi <- class2tree(species_info_ncbi, check = TRUE)
tree_ncbi <- species_tree_ncbi$phylo

#save.image("table_and_phy_ready.RDS")
load("table_and_phy_ready.RDS")
### STATES ###
first_quartile <- round(summary(hedgesg)[2], 2)
second_quartile <- round(summary(hedgesg)[3], 2)

seeds <- c(2, 3, 4)
seeds_phylogeny_rep <- c(100, 101)

states_choice <- c("one") #can be one or two
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
} else {
  message("Error! Chose 'one' or 'two'")
}

hedgesg <- sapply(hedgesg, states)
hedgesg <- setNames(hedgesg, result_all_species$species_complete)

for(rep in seq_along(seeds_phylogeny_rep)){
  set.seed(seeds_phylogeny_rep[rep])
  ### manual corrections in phylogeny ###
  info_fix_poly <- build_info(species, tree_ncbi, db="ncbi")
  input_fix_poly <- info2input(info_fix_poly, tree_ncbi)
  resolved_tree_ncbi <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                                 forceultrametric=TRUE)
  resolved_tree_ncbi$tip.label <- gsub("_", " ", resolved_tree_ncbi$tip.label)
  #save.image('full_and_phy_ready.RDS')
  #load('full_and_phy_ready.RDS')
  resolved_tree_ncbi <- fix.poly(resolved_tree_ncbi, type='resolve')
  for(i in seq_along(seeds)){
    set.seed(seeds[i])
    ### MUSSE CALCULATION NCBI ###
    
    ### musse ###
    # 1 musse full #
    musse_full <- make.musse(resolved_tree_ncbi, states = hedgesg, k = 3)
    init <- starting.point.musse(resolved_tree_ncbi, k = 3)
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
    
    save.image(paste0("output/", "markov", "_", seeds[i], "_", "stat",
                      "_", states_choice, "_", "phy", "_", seeds_phylogeny_rep[rep],
                      "_", "envi.RDS"))
    
    ###### Bayesian MCMC to find posterior density #######
    prior <- make.prior.exponential(1/2)
    
    preliminar <- mcmc(musse_full, 
                       result_musse_full$par, 
                       nsteps=100, prior=prior,
                       w=1, print.every = 0)
    
    w <- diff(sapply(preliminar[2:(ncol(preliminar) -1)], quantile, c(0.05, 0.95)))
    
    mcmc_result <- mcmc(musse_full, 
                        init[colnames(w)], 
                        nsteps=100,prior=prior, 
                        w=w, print.every=10)
    
    write.csv2(mcmc_result, file=paste0("output/","markov", "_",
                                        seeds[i], "_",
                                        "stat", "_",
                                        states_choice, "_",
                                        "phy", "_", seeds_phylogeny_rep[rep],
                                        "mcmc.csv"),
               row.names = FALSE)
  }
}
read_csv2("output/markov_4_stat_one_phy_100mcmc.csv")
mcmc_max <- nrow(mcmc_result)
mcmc_out_burn_in <- round(nrow(mcmc_result) * 0.2) + 1
mcmc_result <- mcmc_result[mcmc_out_burn_in:mcmc_max, ]
#save.image("mcmc.rds")
#load("mcmc.rds")
#mcmc_result <- readRDS("rep_4_stat_one_mcmc.RDS")
####### effective size sample ########
n.eff(as.matrix(mcmc_result[, 2:(length(mcmc_result) - 1)]))
n.eff(as.matrix(mcmc_result[, 2:4]))
n.eff(as.matrix(mcmc_result[, 5:7]))
n.eff(as.matrix(mcmc_result[, 8:13]))

###### BAYES FACTOR CALCULATION ######
bf_mean <- function(x, y) x / y
bf_timestep <- function(x, y) mean(x / y)
mean_posteriors <- colMeans(mcmc_result)[2:ncol(mcmc_result)]
# lamb 1 x lamb 2 #
bf_mean(mean_posteriors[1], mean_posteriors[2])
bf_timestep(mcmc_result[, 2], mcmc_result[, 3])
# lamb 1 x lamb 3 #
bf_mean(mean_posteriors[1], mean_posteriors[3])
bf_timestep(mcmc_result[, 2], mcmc_result[, 4])
# lamb 2 x lamb 3 #
bf_mean(mean_posteriors[2], mean_posteriors[3])
bf_timestep(mcmc_result[, 3], mcmc_result[, 4])

# mu 1 x mu 2 #
bf_mean(mean_posteriors[4], mean_posteriors[5])
bf_timestep(mcmc_result[, 5], mcmc_result[, 6])
# mu 1 x mu 3 #
bf_mean(mean_posteriors[4], mean_posteriors[6])
bf_timestep(mcmc_result[, 5], mcmc_result[, 7])
# mu 2 x mu 3 #
bf_mean(mean_posteriors[5], mean_posteriors[6])
bf_timestep(mcmc_result[, 6], mcmc_result[, 7])

# transition #
# q12 x q13 #
bf_mean(mean_posteriors[7], mean_posteriors[8])
bf_timestep(mcmc_result[, 8], mcmc_result[, 9])
# q12 x q21 #
bf_mean(mean_posteriors[7], mean_posteriors[9])
bf_timestep(mcmc_result[, 8], mcmc_result[, 10])
# q12 x q23 #
bf_mean(mean_posteriors[7], mean_posteriors[10])
bf_timestep(mcmc_result[, 8], mcmc_result[, 11])
# q12 x q31 #
bf_mean(mean_posteriors[7], mean_posteriors[11])
bf_timestep(mcmc_result[, 8], mcmc_result[, 12])
# q12 x q32 #
bf_mean(mean_posteriors[7], mean_posteriors[12])
bf_timestep(mcmc_result[, 8], mcmc_result[, 13])
# q13 x q21 #
bf_mean(mean_posteriors[8], mean_posteriors[9])
bf_timestep(mcmc_result[, 9], mcmc_result[, 10])
# q13 x q23 #
bf_mean(mean_posteriors[8], mean_posteriors[10])
bf_timestep(mcmc_result[, 9], mcmc_result[, 11])
# q13 x q31 #
bf_mean(mean_posteriors[8], mean_posteriors[11])
bf_timestep(mcmc_result[, 9], mcmc_result[, 12])
# q13 x q32 #
bf_mean(mean_posteriors[8], mean_posteriors[12])
bf_timestep(mcmc_result[, 9], mcmc_result[, 13])
# q21 x q23 #
bf_mean(mean_posteriors[9], mean_posteriors[10])
bf_timestep(mcmc_result[, 10], mcmc_result[, 11])
# q21 x q31 #
bf_mean(mean_posteriors[9], mean_posteriors[11])
bf_timestep(mcmc_result[, 10], mcmc_result[, 12])
# q21 x q32 #
bf_mean(mean_posteriors[9], mean_posteriors[12])
bf_timestep(mcmc_result[, 10], mcmc_result[, 13])
# q23 x q31 #
bf_mean(mean_posteriors[10], mean_posteriors[11])
bf_timestep(mcmc_result[, 11], mcmc_result[, 12])
# q23 x q32 #
bf_mean(mean_posteriors[10], mean_posteriors[12])
bf_timestep(mcmc_result[, 11], mcmc_result[, 13])
# q31 x q32 #
bf_mean(mean_posteriors[11], mean_posteriors[12])
bf_timestep(mcmc_result[, 12], mcmc_result[, 13])

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
library(bayou)
prior <- make.prior(tree_to_bms$phy)
bayou.makeMCMC(tree=tree_a, dat=Trait_a, SE = 0, 
               model = "OU", prior=prior )

########## organizing data to graphs with pivot ########
# transitions #
transitions <- mcmc_result %>%
  pivot_longer(q12:q32)

# speciation, extinction, and diversification #
mcmc_result_pivoted <- mcmc_result %>% 
  pivot_longer(starts_with('lambda'),
               names_to="Speciation",
               values_to = 'LambdaPosterior')

temporary <- mcmc_result %>% 
  pivot_longer(starts_with('mu'),
               names_to="Extinction",
               values_to = 'ExtinctionPosterior')

mcmc_result_pivoted$Extinction <- temporary$Extinction
mcmc_result_pivoted$ExtinctionPosterior <- temporary$ExtinctionPosterior

mcmc_result_pivoted <- mcmc_result_pivoted |>
  mutate(Diversification = Speciation,
         DiversificationPosterior = LambdaPosterior - ExtinctionPosterior,
         Diversification = str_replace(Diversification, 'lambda', 'Diversif'))

#################### plots ############################
# temporal series #
save_result <- TRUE
markov <- ggplot(mcmc_result, aes(i, p)) +
  geom_line() + 
  xlab('Time') + ylab('Log(L)') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('markov_result', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(markov)
  dev.off()
}

markov_lambda <- ggplot(mcmc_result, aes(x=i)) +
  geom_line(aes(y=lambda1), color='blue') +
  geom_line(aes(y=lambda2), color='black') +
  geom_line(aes(y=lambda3), color='red') +
  xlab('Time') + ylab('Lambda') + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('markov_result_lambda', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(markov_lambda)
  dev.off()
}

# transitions #
transitions <- ggplot(transitions, aes(value, fill = name)) +
  geom_density(alpha=0.7) +
  theme_bw() +
  scale_fill_hue(name="States") +
  theme(legend.position = c(0.8, 0.75),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  xlab("Transition") + ylab('Posterior Density') 

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('transition', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(transitions)
  dev.off()
}

# speciation #
speciation <- mcmc_result_pivoted %>% 
  ggplot(aes(LambdaPosterior, fill = Speciation)) + 
  geom_density(alpha=0.7) +
  theme_bw() +
  scale_fill_hue(labels = c("Low Plasticity (1)", 
                            "Medium Plasticity (2)", 
                            "High Plasticity (3)"),
                 name="States") +
  theme(legend.position = c(0.8, 0.8),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  xlab("Speciation") + ylab('Posterior Density') 

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('speciation', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(speciation)
  dev.off()
}

# extinction #
extinction <- mcmc_result_pivoted %>% 
  ggplot(aes(ExtinctionPosterior, fill = Extinction)) + 
  geom_density(alpha=0.7) +
  theme_bw() +
  scale_fill_hue(labels = c("Low Plasticity (1)", 
                            "Medium Plasticity (2)", 
                            "High Plasticity (3)"),
                 name="States") +
  theme(legend.position = c(0.8, 0.8),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  xlab("Extinction") + ylab('Posterior Density') 

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('extinction', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(extinction)
  dev.off()
}

# diversification #
diversification <- mcmc_result_pivoted %>% 
  ggplot(aes(DiversificationPosterior, fill = Diversification)) + 
  geom_density(alpha=0.7) +
  theme_bw() +
  scale_fill_hue(labels = c("Low Plasticity (1)", 
                            "Medium Plasticity (2)", 
                            "High Plasticity (3)"),
                 name="States") +
  theme(legend.position = c(0.8, 0.8),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  xlab("Diversification") + ylab('Posterior Density') 

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('diversification', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(diversification)
  dev.off()
}

# tree per states #
color.palette = colorRampPalette(c("#ffcd74", "#ff7251", "#9b2948"))

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('tree_states', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(plotBranchbyTrait(resolved_tree_ncbi, 
                          hedgesg, mode=c("edges"),
                          palette=color.palette,
                          legend=FALSE))
  dev.off()
}

### histogram ###
table_histogram <- as.data.frame(hedgesg)
ggplot(table_histogram, aes(x = hedgesg)) +
 geom_bar() +
 scale_x_continuous(breaks = 1:16) +
 labs(x="Hedge's g effect", y= "Frequency")
