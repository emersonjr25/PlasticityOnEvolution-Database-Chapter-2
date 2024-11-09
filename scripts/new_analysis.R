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
library(readxl)

###################################################################
#### FIRST WAY ####
### READING AND TREATING DATA - BAMM ###
###################################################################
load("table_and_phy_ready.RDS")

result$hedgesg <- abs(result$hedgesg)
result_all_species <- result %>%
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))
hedgesg <- abs(result_all_species$hedgesg)

hedgesg <- setNames(hedgesg, result_all_species$species_complete)

rates <- read.csv2("data/raw/alldat.csv", sep=',')

rates$treename <- gsub("_", " ", rates$treename)

states <- data.frame(hedgesg)
states <- rownames_to_column(states, var='species')
states$species <- tolower(states$species)

rates <- rates |>
  select(treename, bamm, clads, massRate)

names(rates)[1] <- "species"
rates$species <- tolower(rates$species)

result_final_bamm <- inner_join(states, rates, by='species')
result_final_bamm$bamm <- as.numeric(result_final_bamm$bamm)
result_final_bamm$clads <- as.numeric(result_final_bamm$clads)
result_final_bamm$massRate <- as.numeric(result_final_bamm$massRate)

plot(log(result_final_bamm$hedgesg), log(result_final_bamm$bamm))
plot(log(result_final_bamm$hedgesg), log(result_final_bamm$clads))
plot(log(result_final_bamm$hedgesg), log(result_final_bamm$massRate))

plot(result_final_bamm$hedgesg, result_final_bamm$bamm)
plot(result_final_bamm$hedgesg, result_final_bamm$clads)
plot(result_final_bamm$hedgesg, result_final_bamm$massRate)


rates <- read_excel("data/raw/REP_BAMM.xlsx", sheet = 2)
rates <- read_excel("data/raw/REP_DivRate_ES.xlsx", sheet = 2)
rates <- read_excel("data/raw/REP_DivRate_FP.xlsx", sheet = 2)
rates <- read_excel("data/raw/REP_BRANCH LENGTHS.xlsx", sheet = 2)

rates$Species<- gsub("_", " ", rates$Species)

states <- data.frame(hedgesg)
states <- rownames_to_column(states, var='species')
states$species <- tolower(states$species)

rates <- rates |>
  select(Species, Mean)

names(rates)[1] <- "species"
rates$species <- tolower(rates$species)

result_final_bamm <- inner_join(states, rates, by='species')

plot(result_final_bamm$hedgesg, result_final_bamm$Mean)
plot(log(result_final_bamm$hedgesg), log(result_final_bamm$Mean))

