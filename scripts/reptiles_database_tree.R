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

# reading #
reptiles_names <- readxl::read_excel("data/raw/Usados/reptile_checklist_2023_07.xlsx")
species_names <- reptiles_names$Species
load("full_and_phy_ready.RDS")

tree_quali <- readRDS("output/tree/tree_quali.RDS")
tree_prune <- readRDS("output/tree/tree_prune.RDS")

tiff("tree_prune.jpg",
     width = 800, height = 600, res=100)
ggtree(tree_prune) + 
  geom_tiplab(align=TRUE, 
              linesize=.1, size = 2.5) + 
  geom_tree() + theme_tree() + 
  hexpand(.5)
dev.off()

### species focal and complete phylogeny - list + prune ###
tree <- readRDS('resolved_tree_ncbi.RDS')
info_fix_poly <- build_info(species, 
                            tree,
                            find.ranks=TRUE, db="ncbi",
                            mode="list")
input_fix_poly <- info2input(info_fix_poly, tree)
tree_final <- rand_tip(input = input_fix_poly, tree = tree,
                       forceultrametric=TRUE,
                       prune=TRUE)

### species total and incomplete phylogeny - backbone + without prune ###
info_fix_poly <- readRDS('info_fix_poly.RDS')
input_fix_poly <- readRDS('input_fix_poly.RDS')

resolved_tree_ncb2 <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                               forceultrametric=FALSE,
                               prune=FALSE)
resolved_tree_ncbi3 <- rand_tip(input = input_fix_poly, tree = tree_ncbi,
                                forceultrametric=FALSE,
                                prune=TRUE)
saveRDS(resolved_tree_ncbi1, file="tree_without_prune.RDS")
