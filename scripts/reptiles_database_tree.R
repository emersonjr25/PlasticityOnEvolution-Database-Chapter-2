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

# reading #
reptiles_names <- readxl::read_excel("data/raw/Usados/reptile_checklist_2023_07.xlsx")
species_names <- reptiles_names$Species

### species focal and complete phylogeny ###
tree_one <- read.tree("data/raw/Usados/phyliptree.phy")
tree_two <- readRDS('resolved_tree_ncbi.RDS')
tree_three <- readRDS('testINFO.RDS')
tree_four <- readRDS('tree_general.RDS')

### species total and incomplete phylogeny ###
info_fix_poly <- readRDS('info_fix_poly.RDS')
input_fix_poly <- readRDS('input_fix_poly.RDS')

info_fix_poly <- build_info(unique(species_names), 
                            tree_ncbi,
                            find.ranks=TRUE, db="ncbi")
input_fix_poly <- info2input(info_fix_poly, tree_ncbi)

resolved_tree_ncbi <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                               forceultrametric=TRUE,
                               prune=TRUE)


info_fix_poly <- build_info(species, 
                            tree_one,
                            find.ranks=TRUE, db="ncbi")
saveRDS(info_fix_poly, "testINFO.RDS")
input_fix_poly <- info2input(info_fix_poly, tree_one)
resolved_tree_ncbi <- rand_tip(input = input_fix_poly, tree = tree_one,
                               forceultrametric=TRUE)

saveRDS(resolved_tree_ncbi, "tree_general.RDS")

info_fix_poly <- build_info(species, 
                            resolved_tree_ncbi,
                            find.ranks=TRUE, db="ncbi",
                            mode="list")
input_fix_poly <- info2input(info_fix_poly, resolved_tree_ncbi)

tree_final <- rand_tip(input = input_fix_poly, tree = resolved_tree_ncbi,
                               forceultrametric=TRUE,
                               prune=TRUE)
tree <- force.ultrametric(tree_one)
