#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Database Developmental plasticity thermal in reptiles 
##### Script to input data, modify and save plasticity value

### PACKAGES ####
library(dplyr)
#library(stringr)
#library(effectsize)
#library(purrr)
library(here)

### READING DATA ###

dados <- read.csv("data/raw/Database.csv", header = TRUE, sep = ",")
subdata <- dados %>% select(paper_no, genus, species, population, simp_trait, T, mean)
subdata$species_complete <- paste0(subdata$genus," ", subdata$species)

### To calculate Hedge's we need a list with studies, species, traits, temperature and traits ####

index_per_study <- unique(subdata$paper_no)

list_information_studies <- vector("list", length(index_per_study))
species_per_study <- vector("list", length(index_per_study))
name_traits_per_study <- vector("list", length(index_per_study))
pos_species_in_study <-  vector("list", length(index_per_study))

for(i in seq_along(index_per_study)){
  list_information_studies[[i]] <- subdata[index_per_study[i] == subdata$paper_no, ]
}

for(i in seq_along(list_information_studies)){
  for(j in seq_along(list_information_studies[[i]])){
    species_per_study[[i]] <- unique(list_information_studies[[i]]$species_complete)
    name_traits_per_study[[i]] <- unique(list_information_studies[[i]]$simp_trait)
  }
}

for(i in 1:length(name_traits_per_study)){
  pos_species_in_study[[i]] <- vector("list", length(species_per_study[[i]]))
}

for(i in seq_along(pos_species_in_study)){
  for(j in seq_along(pos_species_in_study[[i]])){
    pos_species_in_study[[i]][[j]] <- lapply(name_traits_per_study[[i]], function(x) x)
    names(pos_species_in_study[[i]][[j]]) <- name_traits_per_study[[i]]
  }
}

##### LIST READY TO SAVE VALUES ####

only_trait <- vector("list", length(index_per_study))

for(i in 1:length(name_traits_per_study)){
  only_trait[[i]] <- vector("list", length(name_traits_per_study[[i]]))
}


only_trait <- lapply(list_information_studies, FUN = function(dados){
  logic <- names(dados) == "simp_trait"
 return(dados[logic])
})

pos_trait <- list()

for(i in 1:length(name_traits_per_study)){
  for(j in seq_along(name_traits_per_study[[i]])){
    pos_trait[[i]][[j]] <- which(name_traits_per_study[[i]][[j]] == only_trait[[i]])
  }
}

which(name_traits_per_study[[1]][[1]] == only_trait[[1]])



for(i in seq_along(pos_species_in_study)){
  for(j in seq_along(pos_species_in_study[[i]])){
    pos_each_paper <- which(index_per_study[i] == subdata$paper_no)
    pos_species_in_study[[i]][[j]] <- which(subdata$species_complete[pos_each_paper] == species_per_study[[i]][[j]])
  } 
}

pos_trait_per_specie <- list()

for(i in 1:length(species_per_study)){
  pos_trait_per_specie[[i]] <- vector("list", length(species_per_study[[i]]))
}
names(pos_trait_per_specie) <- unique(subdata$paper_no)

for(i in seq_along(pos_trait_per_specie)){
  names(pos_trait_per_specie[[i]]) <- species_per_study[[i]]
}

for(i in seq_along(pos_trait_per_specie)){
  for(j in seq_along(pos_trait_per_specie[[i]])){
    pos_trait_per_specie[[i]][[j]] <- lapply(name_traits_per_study[[i]], names)
  }
}

for(i in seq_along(pos_trait_per_specie)){
  for(j in seq_along(pos_trait_per_specie[[i]])){
    names(pos_trait_per_specie[[i]][[j]]) <- name_traits_per_study[[i]]
  }
}

#species_in_time <- 0

# opa <- list()
# for(i in seq_along(pos_trait)){
#   for(j in seq_along(pos_trait[[i]])){
#     for(k in seq_along(pos_trait_per_specie[[i]][[j]])){
#         opa[[k]] <- pos_trait[[i]][[j]][pos_species_in_study[[i]][[j]]]
#     }
#   }
# }
# 
# 
# pos_trait[[6]][[1]][pos_species_in_study[[6]][[1]]]


### FILTERING TO FIND HEDGE's per species in studies ###

### MINIMUM ###

minimo <- list()

for(i in 1:length(name_traits)){
  minimo[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    minimo[[i]][[j]] <- min(list_information_studies[[i]]$T[pos_trait[[i]][[j]]])
  }
}

### MAXIMUM ###

maximo <- list()

for(i in 1:length(name_traits)){
  maximo[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    maximo[[i]][[j]] <- max(list_information_studies[[i]]$T[pos_trait[[i]][[j]]])
  }
}

### MIN POSITION ###

posmin <- list()

for(i in 1:length(name_traits)){
  posmin[[i]] <- vector("list", length(name_traits[[i]]))
}

separada <- 0

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    separada <- list_information_studies[[i]]$T
    separada[pos_trait[[i]][[j]]] <- 1
    separada[separada != 1] <- 0
    separada[separada == 1] <- list_information_studies[[i]]$T[pos_trait[[i]][[j]]]
    posmin[[i]][[j]] <- which(separada == minimo[[i]][[j]])
  }
}

### MAX POSITION ###

posmax <- list()

for(i in 1:length(name_traits)){
  posmax[[i]] <- vector("list", length(name_traits[[i]]))
}

separada <- 0

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    separada <- list_information_studies[[i]]$T
    separada[pos_trait[[i]][[j]]] <- 1
    separada[separada != 1] <- 0
    separada[separada == 1] <- list_information_studies[[i]]$T[pos_trait[[i]][[j]]]
    posmax[[i]][[j]] <- which(separada == maximo[[i]][[j]])
  }
}

### MIN MEAN ###

mean_min <- list()

for(i in 1:length(name_traits)){
  mean_min[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    mean_min[[i]][[j]] <- mean(as.numeric(list_information_studies[[i]]$mean[posmin[[i]][[j]]]))
  }
}

### MAX MEAN ###

mean_max <- list()

for(i in 1:length(name_traits)){
  mean_max[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    mean_max[[i]][[j]] <- mean(as.numeric(list_information_studies[[i]]$mean[posmax[[i]][[j]]]))
  }
}

### GROUP 1 ###

group_1 <- list()

for(i in 1:length(name_traits)){
  group_1[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    group_1[[i]][[j]] <- as.numeric(list_information_studies[[i]]$mean[posmin[[i]][[j]]])
  }
}

### GROUP 2 ###

group_2 <- list()

for(i in 1:length(name_traits)){
  group_2[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    group_2[[i]][[j]] <- as.numeric(list_information_studies[[i]]$mean[posmax[[i]][[j]]])
  }
}

### HEDGES's g EFFECT###

hedge_effect <- list()

for(i in 1:length(name_traits)){
  hedge_effect[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    hedge_effect[[i]][[j]] <- abs((mean_min[[i]][[j]] - mean_max[[i]][[j]]) / sd_pooled(group_1[[i]][[j]], group_2[[i]][[j]]))
  }
}

hedge_effect_mean <- list()

for(i in 1:length(name_traits)){
  hedge_effect_mean[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(hedge_effect)){
    hedge_effect_mean[[i]] <- lapply(hedge_effect[[i]], is.na)
}

hedge_effect_out_na <- list()

for(i in 1:length(name_traits)){
  hedge_effect_out_na[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(hedge_effect)){
  pos <- which(hedge_effect_mean[[i]] == FALSE)
  hedge_effect_out_na[[i]] <- as.numeric(hedge_effect[[i]][pos])
  pos <- 0
}

hedge_effect_mean_mean <- list()

for(i in seq_along(hedge_effect)){
  hedge_effect_mean_mean <- lapply(hedge_effect_out_na, mean)
}

pos <- 0
posicao_final <- list()

for(i in 1:length(index_per_study)){
  pos <- which(index_per_study[i] == subdata$paper_no)
  posicao_final[[i]] <- pos[1]
}

posicao_final <- unlist(posicao_final)

names_species <- subdata$species_complete[posicao_final]

names(hedge_effect_mean_mean) <- names_species

#### ORGANIZING TO USE HEDGE's G ###

data_hedges <- do.call(cbind.data.frame, hedge_effect_mean_mean)

L <- split.default(data_hedges, f = gsub("(.*)\\..*", "\\1", names(data_hedges)))

L2 <- list()

L2 <- lapply(L, as.numeric)

L3 <- list()

L3 <- lapply(L2, is.infinite)

L4 <- list()

for(i in seq_along(L2)){
  pos <- which(L3[[i]] == FALSE)
  L4[[i]] <- L2[[i]][pos]
}

L5 <- list()

L5 <- lapply(L2, FUN = function(dados){
    return(is.na(dados))
})

L6 <- list()

for(i in seq_along(L3)){
  pos <- which(L5[[i]] == FALSE)
  L6[[i]] <- L4[[i]][pos]
}

L7 <- list()
L7 <- lapply(L6, mean)

L8 <- list()
L8 <- lapply(L7, is.na)

final_list_hedges <- list()

for(i in seq_along(L3)){
  pos <- which(L8[[i]] == FALSE)
  final_list_hedges[[i]] <- L7[[i]][pos]
}

names(final_list_hedges) <- names(L2)

final_list_hedges2 <- final_list_hedges[lapply(final_list_hedges, length) > 0]

final_list_hedges3 <- final_list_hedges2[final_list_hedges2 != 0]

hedges_g_to_use <- do.call(cbind.data.frame, final_list_hedges3)

path <- here("output")

write.csv2(hedges_g_to_use, file.path(path, "hedges_g.csv"), row.names = FALSE)
