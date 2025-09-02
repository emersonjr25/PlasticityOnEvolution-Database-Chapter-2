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
result_final_bamm$tip.diversification <- result_final_bamm$tip.lambda.avg - result_final_bamm$tip.mu.avg

result_final_bamm <- result_final_bamm[!is.na(result_final_bamm[, "tip.lambda.avg"]),]

# ---- calcular limites ----
q1_lambda <- quantile(result_final_bamm$tip.lambda.avg, 0.25, na.rm = TRUE)
q3_lambda <- quantile(result_final_bamm$tip.lambda.avg, 0.75, na.rm = TRUE)
iqr_lambda <- q3_lambda - q1_lambda
lower_lambda <- q1_lambda - 1.5 * iqr_lambda
upper_lambda <- q3_lambda + 1.5 * iqr_lambda

q1_mu <- quantile(result_final_bamm$tip.mu.avg, 0.25, na.rm = TRUE)
q3_mu <- quantile(result_final_bamm$tip.mu.avg, 0.75, na.rm = TRUE)
iqr_mu <- q3_mu - q1_mu
lower_mu <- q1_mu - 1.5 * iqr_mu
upper_mu <- q3_mu + 1.5 * iqr_mu

# ---- filtrar todos os outliers ----
result_final_bamm <- result_final_bamm %>%
  filter(
    tip.lambda.avg >= lower_lambda & tip.lambda.avg <= upper_lambda,
    tip.mu.avg >= lower_mu & tip.mu.avg <= upper_mu
  )

subdata <- dados %>% 
  as.tibble() %>%
  select(paper_no, genus, order, family, species, simp_trait, T, mean) %>%
  mutate(species_complete = paste0(genus," ", species)) %>%
  group_by(paper_no, species_complete, simp_trait) %>%
  mutate(id = paste0('ID', paper_no, species_complete, simp_trait)) %>%
  select(id, paper_no, species_complete, order, family,
         simp_trait, T, mean)  %>%
  mutate(species_complete = str_squish(species_complete))

subdata$species <- tolower(subdata$species_complete)
subdata <- subdata |>
  ungroup() |>
  select(species, order) |> 
  distinct()

result_final_bamm <- left_join(result_final_bamm, subdata, by='species')

result_final_bamm <- result_final_bamm %>%
  filter(order != "Crocodilia")

result_final_bamm |>
  group_by(hedgesg, order) |>
  summarise(round(mean(tip.mu.avg), 9))

result_final_bamm |>
  group_by(hedgesg, order) |>
  summarise(round(mean(tip.lambda.avg), 9))

result_final_bamm |>
  group_by(hedgesg, order) |>
  summarise(mean(tip.diversification))

p1 <- result_final_bamm %>% 
  ggplot(aes(x = hedgesg, y = tip.lambda.avg, color = order)) +
  geom_jitter(alpha = 0.7, size = 2) + 
  geom_smooth(aes(group = order), method = "lm", se = FALSE) +  # regression line per order
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
  ggplot(aes(x = hedgesg, y = tip.mu.avg, color = order)) +
  geom_jitter(alpha = 0.7, size = 2) + 
  geom_smooth(aes(group = order), method = "lm", se = FALSE) +  # regression line per order
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
  ggplot(aes(x = hedgesg, y = tip.diversification, color = order)) +
  geom_jitter(alpha = 0.7, size = 2) + 
  geom_smooth(aes(group = order), method = "lm", se = FALSE) +  # regression line per order
  theme_bw() +
  xlab("Hedges' g") +
  ylab("Diversification rate") +
  theme(
    plot.title = element_text(size = 14, face = 2, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "bottom"
  )

ggsave("speciation_rate_order.tiff", plot = p1, width = 6, height = 4, dpi = 100)
ggsave("extinction_rate_order.tiff", plot = p2, width = 6, height = 4, dpi = 100)
ggsave("diversification_rate_order.tiff", plot = p3, width = 6, height = 4, dpi = 100)

write.csv(result_final_bamm, 
          file = "result_final_bamm.csv", 
          row.names = FALSE)
