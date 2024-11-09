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

### carrying all data ###
load("table_and_phy_ready.RDS")

### verifying more frequent trait ###
trait_frequent <- result %>%
  group_by(simp_trait) %>%
  summarise(len = length(simp_trait)) %>%
  arrange(desc(len))

# selecting data to more frequent trait - mass #
more_frequent_mass <- trait_frequent[grep("ass",trait_frequent$simp_trait),]
five_more_frequent_mass <- more_frequent_mass[1, 1]
subdata$units <- dados$units
subdata_filtered_mass <- subdata[subdata$simp_trait %in% unclass(five_more_frequent_mass)$simp_trait, ]

# filter and transform data to mass #
subdata_filtered_mass$units <- gsub("\\s", "", subdata_filtered_mass$units)
subdata_filtered_mass <- subdata_filtered_mass |>
  filter(units == "g" | units == "mg")
subdata_filtered_mass[subdata_filtered_mass$units == "mg", ]$mean <- subdata_filtered_mass[subdata_filtered_mass$units == "mg", ]$mean / 1000

# mass per species #
table_bms <- subdata_filtered_mass |> 
  select(-c(units)) |>
  group_by(paper_no, species_complete) |>
  summarise(mean_value_trait = mean(mean)) |>
  ungroup() |>
  group_by(species_complete) |>
  summarise(trait_value = mean(mean_value_trait))
  
unique(subdata_filtered_mass$species_complete)

# hedges g to bms #
result$hedgesg <- abs(result$hedgesg)
result_mass <- result %>%
  group_by(species_complete) %>%
  summarise(hedgesg = mean(hedgesg))

### putting only species mass that is present in hedgesg table #
table_bms <- table_bms[table_bms$species_complete %in% result_mass$species_complete,]
result_mass <- result_mass[result_mass$species_complete %in% table_bms$species_complete,]
table_bms$trait_value <- round(table_bms$trait_value, 3)
hedgesg_mass <- abs(result_mass$hedgesg)

### STATES ###
first_quartile <- round(summary(hedgesg)[2], 2)
second_quartile <- round(summary(hedgesg)[3], 2)

seeds <- c(2, 3, 4)

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

hedgesg_mass <- sapply(hedgesg_mass, states)
hedgesg_mass <- setNames(hedgesg_mass, result_mass$species_complete)

### removing species with big mass (outlier) from tables ###
species_outlier <- table_bms |> 
  arrange(desc(trait_value)) |>
  select(species_complete) |>
  slice(1)
species_outlier <- unclass(species_outlier)$species_complete

hedgesg_mass <- hedgesg_mass[names(hedgesg_mass) != species_outlier]
table_bms <- table_bms[table_bms$species_complete != species_outlier, ]
result_mass <- result_mass[result_mass$species_complete != species_outlier, ]

species_mass <- unique(result_mass$species_complete)

### reading and correcting phylogeny ###
#write(species_mass, "species_bms_general.txt")
reptiles_tree_time_tree <- read.newick("data/raw/species_bms_general.nwk")

seed_phy <- c(100, 101)
set.seed(seed_phy[1])
info_fix_poly <- build_info(names(hedgesg_mass), reptiles_tree_time_tree, 
                            find.ranks=TRUE, db="ncbi")
input_fix_poly <- info2input(info_fix_poly, reptiles_tree_time_tree,
                             parallelize=F)
tree_time_tree_ready <- rand_tip(input = input_fix_poly, 
                                 tree = reptiles_tree_time_tree,
                                 forceultrametric=TRUE,
                                 prune=TRUE)
tree_time_tree_ready$tip.label <- gsub("_", " ", tree_time_tree_ready$tip.label)

### seeing distribuition of the data ###
boxplot(table_bms$trait_value)
boxplot(result_mass$hedgesg)

### Preparing table data to BMS ###
set.seed(seeds[3])
hedgesg_mass <- as.factor(hedgesg_mass)
X_to_BMS <- abs(table_bms$trait_value)
X_to_BMS <- setNames(X_to_BMS, table_bms$species_complete)
Trait <- data.frame(Genus_species = names(hedgesg_mass),
                    Reg = as.numeric(hedgesg_mass),
                    X = as.numeric(X_to_BMS))

tree_to_ouwie <- make.simmap(tree_time_tree_ready, hedgesg_mass, model="ER")
plot(tree_to_ouwie)

# BMS #
bms <- OUwie(tree_to_ouwie,Trait,model="BMS", simmap.tree=TRUE, root.station=FALSE)
bms
