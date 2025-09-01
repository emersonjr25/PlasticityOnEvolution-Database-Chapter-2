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
#library(bayou)
#library(ggtree)
library(dispRity)
library(BAMMtools)
library(hisse)

### saving phylogeny ####
#tree_time_tree_ready <- read.newick("data/raw/repteis_especies_total.nwk") #repteis
#tree_time_tree_ready <- force.ultrametric(tree_time_tree_ready, method="extend")
#tree_time_tree_ready <- dispRity::remove.zero.brlen(tree_time_tree_ready)
#is.ultrametric(tree_time_tree_ready)
#is.binary.tree(tree_time_tree_ready)
#min(tree_time_tree_ready$edge.length)
#write.tree(tree_time_tree_ready, "filogenia_repteis.tre")

##### reading phylogeny ######
tree_time_tree_ready <- read.tree('filogenia_repteis.tre')
mcmcout <- read.table("mcmc_out_reptiles_2milion.txt", sep=',', header = T)
event_data <- read.table("event_data_reptiles_2milion.txt", sep=',', header = T)

plot(mcmcout$logLik ~ mcmcout$generation, pch = 19)
burnstart <- floor(0.8 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, pch = 19)
#conver <- sapply(postburn, effectsize::effectiveSize)

tree_time_tree_ready$tip.label <- gsub(" ", "_", tree_time_tree_ready$tip.label)

results <- BAMMtools::getEventData(tree_time_tree_ready, 
                                   event_data, burnin=0.8,
                                   nsamples=1000)

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

#teste <- result_final_bamm |>
#  pivot_wider(names_from=hedgesg, values_from = starts_with("tip.lambda"))

result_final_bamm <- result_final_bamm[(result_final_bamm$tip.lambda.avg >= summary(result_final_bamm$tip.lambda.avg)[2] & 
                                         result_final_bamm$tip.lambda.avg <= summary(result_final_bamm$tip.lambda.avg)[5]), ] 

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(round(mean(tip.mu.avg), 9))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(round(mean(tip.lambda.avg), 9))

result_final_bamm |>
  group_by(hedgesg) |>
  summarise(mean(tip.diversification))

p1 <- result_final_bamm %>% 
  ggplot(aes(x = hedgesg, y = tip.lambda.avg)) + 
  geom_jitter(alpha = 0.7, size = 2) + 
  theme_bw() +
  xlab("Hedges' g") +
  ylab("Speciation rate") +
  theme(
    plot.title = element_text(size = 14, face = 2, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "bottom"
  )

p2 <- result_final_bamm %>% 
  ggplot(aes(x = hedgesg, y = tip.mu.avg)) + 
  geom_jitter(alpha = 0.7, size = 2) + 
  theme_bw() +
  xlab("Hedges' g") +
  ylab("Extinction rate") +
  theme(
    plot.title = element_text(size = 14, face = 2, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "bottom"
  )

p3 <- result_final_bamm %>% 
  ggplot(aes(x = hedgesg, y = tip.diversification)) + 
  geom_jitter(alpha = 0.7, size = 2) + 
  theme_bw() +
  xlab("Hedges' g") +
  ylab("Diversification rate") +
  theme(
    plot.title = element_text(size = 14, face = 2, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "bottom"
  )

# Save as TIFF (high resolution, e.g., 300 dpi)
ggsave("speciation_rate.tiff", plot = p1, width = 6, height = 4, dpi = 300)
ggsave("extinction_rate.tiff", plot = p2, width = 6, height = 4, dpi = 300)
ggsave("diversification_rate.tiff", plot = p3, width = 6, height = 4, dpi = 300)
write.csv(result_final_bamm, 
          file = "result_final_bamm.csv", 
          row.names = FALSE)
