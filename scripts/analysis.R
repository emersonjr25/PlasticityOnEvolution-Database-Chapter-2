#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Database Developmental plasticity thermal in reptiles 
##### Script to analysis data, modify and save plasticity value
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

#######################################################
########### phylogeny with 3 cohen ways ###############
#######################################################
### state one - 3 rep mcmc ###
phy_expanded_cohen_one_rep_one <- read.csv2("output/markov_files/cohen 1/rep1/phy_expanded_yesstat_one_markov_2_mcmc.csv")
#phy_expanded_cohen_one_rep_two <- read.csv2("output/phy_expanded_yesstat_one_markov_3_mcmc.csv")
#phy_expanded_cohen_one_rep_three <- read.csv2("output/phy_expanded_yesstat_one_markov_4_mcmc.csv")

### state three - 3 rep mcmc ###
#phy_expanded_cohen_three_rep_one <- read.csv2("output/phy_expanded_yesstat_three_markov_2_mcmc.csv")
#phy_expanded_cohen_three_rep_two <- read.csv2("output/phy_expanded_yesstat_three_markov_3_mcmc.csv")
#phy_expanded_cohen_three_rep_three <- read.csv2("output/phy_expanded_yesstat_three_markov_4_mcmc.csv")

state_chosen <- phy_expanded_cohen_one_rep_one

mcmc_max <- nrow(state_chosen)
mcmc_out_burn_in <- round(nrow(state_chosen) * 0.2) + 1
mcmc_result <- state_chosen[mcmc_out_burn_in:mcmc_max, ]

n.eff(as.matrix(state_chosen[, 2:(length(mcmc_result) - 1)]))
n.eff(as.matrix(state_chosen[, 2:4]))
n.eff(as.matrix(state_chosen[, 5:7]))
n.eff(as.matrix(state_chosen[, 8:13]))

###### BAYES FACTOR CALCULATION ######
bf_mean <- function(x, y) x / y
bf_timestep <- function(x, y) mean(x / y)
mean_posteriors <- colMeans(mcmc_result)[2:ncol(mcmc_result)]
round(mean_posteriors, 5)
# lamb 1 x lamb 2 #
bf_mean(mean_posteriors[1], mean_posteriors[2])
# lamb 1 x lamb 3 #
bf_mean(mean_posteriors[1], mean_posteriors[3])
# lamb 2 x lamb 3 #
bf_mean(mean_posteriors[2], mean_posteriors[3])

# mu 1 x mu 2 #
bf_mean(mean_posteriors[4], mean_posteriors[5])
# mu 1 x mu 3 #
bf_mean(mean_posteriors[4], mean_posteriors[6])
# mu 2 x mu 3 #
bf_mean(mean_posteriors[5], mean_posteriors[6])

# transition #
# q12 x q13 #
bf_mean(mean_posteriors[7], mean_posteriors[8])
# q12 x q21 #
bf_mean(mean_posteriors[7], mean_posteriors[9])
# q12 x q23 #
bf_mean(mean_posteriors[7], mean_posteriors[10])
# q12 x q31 #
bf_mean(mean_posteriors[7], mean_posteriors[11])
# q12 x q32 #
bf_mean(mean_posteriors[7], mean_posteriors[12])
# q13 x q21 #
bf_mean(mean_posteriors[8], mean_posteriors[9])
# q13 x q23 #
bf_mean(mean_posteriors[8], mean_posteriors[10])
# q13 x q31 #
bf_mean(mean_posteriors[8], mean_posteriors[11])
# q13 x q32 #
bf_mean(mean_posteriors[8], mean_posteriors[12])
# q21 x q23 #
bf_mean(mean_posteriors[9], mean_posteriors[10])
# q21 x q31 #
bf_mean(mean_posteriors[9], mean_posteriors[11])
# q21 x q32 #
bf_mean(mean_posteriors[9], mean_posteriors[12])
# q23 x q31 #
bf_mean(mean_posteriors[10], mean_posteriors[11])
# q23 x q32 #
bf_mean(mean_posteriors[10], mean_posteriors[12])
# q31 x q32 #
bf_mean(mean_posteriors[11], mean_posteriors[12])

mean_posteriors[7:12]

# speciation, extinction, and diversification #
# transitions #
transitions <- mcmc_result %>%
  pivot_longer(q12:q32)

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

### bf diversification ###
bf_diversif <- mcmc_result_pivoted |>
  select(Diversification, DiversificationPosterior) |>
  group_by(Diversification) |>
  mutate(row = row_number()) |>
  pivot_wider(names_from=Diversification,
              values_from=DiversificationPosterior) |>
  select(-row)
diversification_means <- colMeans(bf_diversif)
# diversification #
# diversi 1 x diversi 2 #
bf_mean(diversification_means[1], diversification_means[2])
# diversi 1 x diversi 3 #
bf_mean(diversification_means[1], diversification_means[3])
# diversi 2 x diversi 3 #
bf_mean(diversification_means[2], diversification_means[3])

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
  scale_x_log10() +
  theme_bw() +
  scale_fill_hue(name="States") +
  theme(legend.position = c(0.1, 0.7),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) +
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
  scale_x_log10() +
  theme_bw() +
  scale_fill_hue(labels = c("Low Plasticity (1)", 
                            "Medium Plasticity (2)", 
                            "High Plasticity (3)"),
                 name="States") +
  theme(legend.position = c(0.2, 0.75),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) +
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
  scale_x_log10() +
  theme_bw() +
  scale_fill_hue(labels = c("Low Plasticity (1)", 
                            "Medium Plasticity (2)", 
                            "High Plasticity (3)"),
                 name="States") +
  theme(legend.position = c(0.2, 0.75),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) +
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
  scale_x_log10() +
  theme_bw() +
  scale_fill_hue(labels = c("Low Plasticity (1)", 
                            "Medium Plasticity (2)", 
                            "High Plasticity (3)"),
                 name="States") +
  theme(legend.position = c(0.2, 0.75),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) +
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

### carrying envi to musse results ###
load("output/markov_files/cohen 1/rep3/phy_expanded_yesstat_one_markov_4_envi.RDS")
aicw(setNames(anova_result$AIC, row.names(anova_result)))
round(coef(result_musse_full), 9)

#### tree general ####
tiff("output/tree_original.jpg",
     width = 800, height = 600, res=100)
ggtree(tree_time_tree_ready, aes(hedgesg)) + 
  geom_tiplab(align=TRUE, 
              linesize=.1, size = 2.5) + 
  geom_tree() + theme_tree() + 
  hexpand(.5)
dev.off()

# tree per states #
color.palette = colorRampPalette(c("#ffcd74", "#ff7251", "#9b2948"))

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('tree_states_original', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(plotBranchbyTrait(tree_time_tree_ready, 
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
