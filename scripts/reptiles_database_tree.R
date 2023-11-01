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

# reading table #
reptiles_names <- readxl::read_excel("data/raw/Usados/reptile_checklist_2023_07.xlsx")
species_names <- reptiles_names$Species
load("full_and_phy_ready.RDS")

### species total and incomplete phylogeny - backbone + without prune ###
# verifying input #
info_fix_poly <- readRDS('data/raw/input_rds_randtip/info_fix_poly.RDS')
input_fix_poly <- readRDS('data/raw/input_rds_randtip/input_fix_poly.RDS')

info_fix_poly <- build_info(unique(species_names), 
                            tree_ncbi,
                            find.ranks=TRUE, db="ncbi",
                            mode="backbone")
input_fix_poly <- info2input(info_fix_poly, tree_ncbi)

resolved_tree_ncb1 <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                               forceultrametric=TRUE,
                               prune=FALSE)
resolved_tree_ncb2 <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                               forceultrametric=FALSE,
                               prune=FALSE)
resolved_tree_ncbi3 <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                                forceultrametric=FALSE,
                                prune=TRUE)
resolved_tree_ncbi4 <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                                forceultrametric=TRUE,
                                prune=TRUE)

### species focal and complete phylogeny - list + prune ###
tree <- readRDS('data/raw/input_rds_randtip/resolved_tree_ncbi4.RDS')

info_fix_poly <- build_info(species, 
                            tree,
                            find.ranks=TRUE, db="ncbi",
                            mode="list")
input_fix_poly <- info2input(info_fix_poly, tree)
tree_final <- rand_tip(input = input_fix_poly, tree = tree,
                       forceultrametric=TRUE,
                       prune=TRUE)

tree_quali <- readRDS("output/tree/tree_quali.RDS")
tree_prune <- readRDS("output/tree/tree_prune.RDS")

tiff("output/tree_prune.jpg",
     width = 800, height = 600, res=100)
ggtree(tree_final) + 
  geom_tiplab(align=TRUE, 
              linesize=.1, size = 2.5) + 
  geom_tree() + theme_tree() + 
  hexpand(.5)
dev.off()
