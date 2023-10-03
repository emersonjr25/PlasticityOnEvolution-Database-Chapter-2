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

# states first way #
# states <- function(x){
#   if(x <= 0.2){
#     x <- 1
#   } else if (x > 0.2 & x <= 0.5){
#     x <- 2
#   } else if(x > 0.5){
#     x <- 3
#   }
# }

# states second way #
states <- function(x){
  if(x <= 0.35){
    x <- 1
  } else if (x > 0.35 & x <= 0.65){
    x <- 2
  } else if(x > 0.65){
    x <- 3
  }
}

# states third way #
# states <- function(x){
#   if(x < 0.5){
#     x <- 1
#   } else if (x >= 0.5 & x < 0.8){
#     x <- 2
#   } else if(x >= 0.8){
#     x <- 3
#   }
# }

hedgesg <- sapply(hedgesg, states)
hedgesg <- setNames(hedgesg, result_all_species$species_complete)
species <- unique(result_all_species$species_complete)

### histogram ###
table_histogram <- as.data.frame(hedgesg)
ggplot(table_histogram, aes(x = hedgesg)) + 
  geom_bar() + 
  scale_x_continuous(breaks = 1:16) + 
  labs(x="Hedge's g effect", y= "Frequency")

#### phylogeny construction - NCBI ####
species_info_ncbi <- classification(species, db='ncbi')
species_tree_ncbi <- class2tree(species_info_ncbi, check = TRUE)
resolved_tree_ncbi <- species_tree_ncbi$phylo

### MUSSE CALCULATION NCBI ###

### manual corrections in phylogeny ###
species_wrong <- names(hedgesg)[!names(hedgesg) %in% resolved_tree_ncbi$tip.label]
species_wrong2 <- resolved_tree_ncbi$tip.label[!resolved_tree_ncbi$tip.label %in% names(hedgesg)]
names(hedgesg)[grepl('ruf', names(hedgesg))] <- species_wrong2[species_wrong2 == "Lycodon rufozonatus"]
resolved_tree_ncbi$tip.label[grepl('Elaphe tae', resolved_tree_ncbi$tip.label)] <- species_wrong[species_wrong == "Elaphe taeniura"]
names(hedgesg)[grepl('maccoyi', names(hedgesg))][1] <- species_wrong2[species_wrong2 == "Anepischetosia maccoyi.x"]
names(hedgesg)[grepl('maccoyi', names(hedgesg))][2] <- species_wrong2[species_wrong2 == "Anepischetosia maccoyi.y"]
names(hedgesg)[grepl('sinensis', names(hedgesg))][1] <- species_wrong2[species_wrong2 == "Mauremys sinensis"]
names(hedgesg)[grepl('piscator', names(hedgesg))] <- species_wrong2[species_wrong2 == "Fowlea piscator"]
names(hedgesg)[grepl('lesueurii', names(hedgesg))][1] <- species_wrong2[species_wrong2 == "Amalosia lesueurii"]
names(hedgesg)[grepl('lesueurii', names(hedgesg))][2] <- species_wrong2[species_wrong2 == "Intellagama lesueurii"]

###################### musse ##########################
resolved_tree_ncbi <- multi2di(resolved_tree_ncbi)
resolved_tree_ncbi <- fix.poly(resolved_tree_ncbi, type='resolve')

# 1 musse full #
musse_full <- make.musse(resolved_tree_ncbi, states = hedgesg, k = 3)
init <- starting.point.musse(resolved_tree_ncbi, k = 3)
result_musse_full <- find.mle(musse_full, x.init = init)
round(result_musse_full$par, 9)

# 2 musse ordered #
musse_ordered <- constrain(musse_full, q13 ~ 0, q31 ~ 0)
result_musse_ordered <- find.mle(musse_ordered, x.init = init[argnames(musse_ordered)])
round(result_musse_ordered$par, 9)

# 3 musse null #
musse_null <- constrain(musse_full, lambda2 ~ lambda1, lambda3 ~ lambda1,
                        mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q21 ~ q12, 
                        q23 ~ q12, q31 ~ 0, q32 ~ q12)
result_musse_null <- find.mle(musse_null, x.init = init[argnames(musse_null)])
round(result_musse_null$par, 9)

# 4 musse lambda #
musse_lambda <- constrain(musse_full, mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q21 ~ q12, 
                        q23 ~ q12, q31 ~ 0, q32 ~ q12)
result_lambda <- find.mle(musse_lambda, x.init = init[argnames(musse_lambda)])
round(result_lambda$par, 9)

# 5 musse mu #
musse_mu <- constrain(musse_full, lambda2 ~ lambda1, lambda3 ~ lambda1, 
          q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
result_mu <- find.mle(musse_mu, x.init = init[argnames(musse_mu)])
round(result_mu$par, 9)

# anova to see best musse model #
anova_result <- anova(result_musse_null,
                      all.different = result_musse_full,
                      free.ordered = result_musse_ordered,
                      free.lambda = result_lambda,
                      free.mu = result_mu)

# look results to best model #
aicw(setNames(anova_result$AIC, row.names(anova_result)))
round(coef(result_musse_ordered), 9)
logLik(result_musse_ordered)
AIC(result_musse_ordered)

#### Bayesian MCMC to find posterior density ####
prior <- make.prior.exponential(1/2)

preliminar <- mcmc(musse_ordered, 
                   result_musse_ordered$par, 
                   nsteps=100, prior=prior,
                   w=1, print.every = 0)

w <- diff(sapply(preliminar[2:(ncol(preliminar) -1)], quantile, c(0.05, 0.95)))

mcmc_result <- mcmc(musse_ordered, 
                    init[colnames(w)], 
                    nsteps=3000,prior=prior, 
                    w=w, print.every=10)
#saveRDS(mcmc_result, file="mcmc.rds")

# remove burn-in #
#mcmc_result <- readRDS("mcmc.rds")
mcmc_max <- nrow(mcmc_result)
mcmc_out_burn_in <- round(nrow(mcmc_result) * 0.2) + 1
mcmc_result <- mcmc_result[mcmc_out_burn_in:mcmc_max, ]

# mean result per variable #
colMeans(mcmc_result)[2:ncol(mcmc_result)]

### organizing data to graphs with pivot ###

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

### plots ###

# temporal series #
save_result <- FALSE
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
ggplot(transitions, aes(value, fill = name)) +
  geom_density(alpha=0.7) +
  theme_bw() +
  scale_fill_hue(name="States") +
  theme(legend.position = c(0.8, 0.75),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  xlab("Transition") + ylab('Posterior Density') 

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

################################################################################
##### MUSSE TO MORE FREQUENT TRAIT #######################
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
  } else if(x > 0.5){
    x <- 3
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
resolved_tree_ncbi_mass <- multi2di(species_tree_ncbi_mass)
species_wrong <- names(hedgesg_mass)[!names(hedgesg_mass) %in% resolved_tree_ncbi_mass$tip.label]
species_wrong2 <- resolved_tree_ncbi_mass$tip.label[!resolved_tree_ncbi_mass$tip.label %in% names(hedgesg_mass)]

### manual corrections in phylogeny ###
names(hedgesg_mass)[grepl('piscator', names(hedgesg_mass))] <- species_wrong2[species_wrong2 == "Fowlea piscator"]
names(hedgesg_mass)[grepl('lesueurii', names(hedgesg_mass))][1] <- species_wrong2[species_wrong2 == "Amalosia lesueurii"]

#### musse ####
resolved_tree_ncbi_mass <- fix.poly(resolved_tree_ncbi_mass, type='resolve')

# 1 musse full #
musse_full <- make.musse(resolved_tree_ncbi_mass, states = hedgesg_mass, k = 3)
init <- starting.point.musse(resolved_tree_ncbi_mass, k = 3)
result_musse <- find.mle(musse_full, x.init = init)
round(result_musse$par, 7)

# 2 musse ordered #
musse_ordered <- constrain(musse_full, q13 ~ 0, q31 ~ 0)
result_musse_ordered <- find.mle(musse_ordered, x.init = init[argnames(musse_ordered)])
round(result_musse_ordered$par, 9)

# 3 musse null #
musse_null <- constrain(musse_full, lambda2 ~ lambda1, lambda3 ~ lambda1,
                        mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q21 ~ q12, 
                        q23 ~ q12, q31 ~ 0, q32 ~ q12)
result_musse_null <- find.mle(musse_null, x.init = init[argnames(musse_null)])
round(result_musse_null$par, 9)

# 4 musse lambda #
musse_lambda <- constrain(musse_full, mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q21 ~ q12, 
                          q23 ~ q12, q31 ~ 0, q32 ~ q12)
result_lambda <- find.mle(musse_lambda, x.init = init[argnames(musse_lambda)])
round(result_lambda$par, 9)

# 5 musse mu #
musse_mu <- constrain(musse_full, lambda2 ~ lambda1, lambda3 ~ lambda1, 
                      q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
result_mu <- find.mle(musse_mu, x.init = init[argnames(musse_mu)])
round(result_mu$par, 9)

# 6 musse lamba-mu #
musse_lambda_mu <- constrain(musse_full, q13 ~ 0, q21 ~ q12, q23 ~ q12, 
                             q31 ~ 0, q32 ~ q12)
result_lambda_mu <- find.mle(musse_lambda_mu, init[argnames(musse_lambda_mu)])
round(result_lambda_mu$par, 9)

# 7 musse transitions #
musse_transitions <- constrain(musse_full, lambda2 ~ lambda1, lambda3 ~ lambda1, 
                               mu2 ~ mu1, mu3 ~ mu1)
result_transitions <- find.mle(musse_transitions, init[argnames(musse_transitions)])
round(result_transitions$par, 9)

# anova to see best musse model #
anova_result <- anova(result_musse_null,
                      all.different = result_musse_full,
                      free.ordered = result_musse_ordered,
                      free.lambda = result_lambda,
                      free.mu = result_mu,
                      free.lambda.mu = result_lambda_mu,
                      free.q = result_transitions)

# look results to best model #
aicw(setNames(anova_result$AIC, row.names(anova_result)))
round(coef(result_musse_ordered), 9)
logLik(result_musse_ordered)
AIC(result_musse_ordered)
plotTree(resolved_tree_ncbi)

#### Bayesian MCMC to best model ####
prior <- make.prior.exponential(1/2)

preliminar <- mcmc(musse_ordered, 
                   result_musse_ordered$par, 
                   nsteps=100, prior=prior,
                   w=1, print.every = 0)

w <- diff(sapply(preliminar[2:(ncol(preliminar) -1)], quantile, c(0.05, 0.95)))

mcmc_result <- mcmc(musse_ordered, 
                    init[colnames(w)], 
                    nsteps=500,prior=prior, 
                    w=w, print.every=10)

plot(mcmc_result$i, mcmc_result$p, type='l', 
     xlab='generation', ylab='log(L)')

mcmc_max <- nrow(mcmc_result)
mcmc_out_burn_in <- nrow(mcmc_result) * 0.2
mcmc_result <- mcmc_result[mcmc_out_burn_in:mcmc_max, ]
colMeans(mcmc_result)[2:ncol(mcmc_result)]
colors <- setNames(c('yellow', 'green', 'red'), 1:3)
profiles.plot(mcmc_result[, grep('lambda', colnames(mcmc_result))],
              col.line=colors, las=1, legend.pos = 'topright')
profiles.plot(mcmc_result[, grep('mu', colnames(mcmc_result))],
              col.line=colors, las=1, legend.pos = 'topright')
net_div <- mcmc_result[, grep('lambda', colnames(mcmc_result))] -
  mcmc_result[, grep('mu', colnames(mcmc_result))]
colnames(net_div) <- paste('lambda-mu(', 1:3, ")", sep="")
profiles.plot(net_div,
              xlab="Net diversification rate",
              ylab="Probability density",
              legend.pos='topleft', col.line=setNames(colors, colnames(net_div)),
              lty=1)

##################################################################
resolved_tree_ncbi <- species_tree_ncbi
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

info_fix_poly <- build_info(species, species_tree_ncbi, db="ncbi")
input_fix_poly <- info2input(info_fix_poly, species_tree_ncbi)
tree_solved <- rand_tip(input = input_fix_poly, tree = species_tree_ncbi,
                        forceultrametric=TRUE)
test <- fix.poly(tree_solved, type='resolve')
test2 <- fix.poly(species_tree_ncbi, type='resolve')
musse <- make.musse(resolved_tree_ncbi, states = hedgesg, k = 3)
init <- starting.point.musse(resolved_tree_ncbi, k = 3)
result_musse <- find.mle(musse, x.init = init)
round(result_musse$par, 9)