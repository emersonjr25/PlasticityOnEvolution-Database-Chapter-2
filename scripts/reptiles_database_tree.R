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
reptiles_names <- readxl::read_excel("data/raw/reptile_checklist_2023_07.xlsx")
species_names <- reptiles_names$Species

info_fix_poly <- readRDS('info_fix_poly.RDS')

info_fix_poly <- build_info(unique(species_names), 
                            tree_ncbi,
                            find.ranks=TRUE, db="ncbi", 
                            mode="backbone")
input_fix_poly <- info2input(info_fix_poly, tree_ncbi)

resolved_tree_ncbi <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                               forceultrametric=FALSE)
