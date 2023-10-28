reptiles_names <- readxl::read_excel("data/raw/reptile_checklist_2023_07.xlsx")

species_names <- reptiles_names$Species
species_info_ott <- tnrs_match_names(species_names)
species_info_ott_without_na <- species_info_ott[!is.na(species_info_ott$ott_id), ]
lista <- list(ott132336 = "pruned_ott_id", ott168143 = "pruned_ott_id", ott2897028 = "pruned_ott_id", ott306228 = "pruned_ott_id", ott3129916 = "pruned_ott_id", ott3247652 = "pruned_ott_id", ott3611538 = "pruned_ott_id", ott4120447 = "pruned_ott_id", ott4120808 = "pruned_ott_id", ott4120816 = "pruned_ott_id", ott4121615 = "pruned_ott_id", ott4122099 = "pruned_ott_id", ott4122151 = "pruned_ott_id", ott4122152 = "pruned_ott_id", ott4122218 = "pruned_ott_id", ott4122219 = "pruned_ott_id", ott4123558 = "pruned_ott_id", 
              ott4123560 = "pruned_ott_id", ott4123561 = "pruned_ott_id", ott4124107 = "pruned_ott_id", ott4124142 = "pruned_ott_id", ott4124143 = "pruned_ott_id", ott4124144 = "pruned_ott_id", ott4124145 = "pruned_ott_id", ott4124146 = "pruned_ott_id", ott4124147 = "pruned_ott_id", ott4124148 = "pruned_ott_id", ott4124149 = "pruned_ott_id", ott4124415 = "pruned_ott_id", ott4124422 = "pruned_ott_id", ott4125009 = "pruned_ott_id", ott4146684 = "pruned_ott_id", ott4434900 = "pruned_ott_id", ott4527868 = "pruned_ott_id", 
              ott457023 = "pruned_ott_id", ott4642289 = "pruned_ott_id", ott4945797 = "pruned_ott_id", ott4945915 = "pruned_ott_id", ott4945916 = "pruned_ott_id", ott4945918 = "pruned_ott_id", ott4945919 = "pruned_ott_id", ott4945920 = "pruned_ott_id", ott4945921 = "pruned_ott_id", ott4945946 = "pruned_ott_id", ott4954046 = "pruned_ott_id", ott5029544 = "pruned_ott_id", ott5034751 = "pruned_ott_id", ott5044891 = "pruned_ott_id", ott5077641 = "pruned_ott_id", ott5078556 = "pruned_ott_id", ott5081239 = "pruned_ott_id", 
              ott5084956 = "pruned_ott_id", ott6156862 = "pruned_ott_id", ott6157467 = "pruned_ott_id", ott6157468 = "pruned_ott_id", ott6157469 = "pruned_ott_id", ott6157880 = "pruned_ott_id", ott6338650 = "pruned_ott_id", ott7069206 = "pruned_ott_id", ott733420 = "pruned_ott_id", ott740044 = "pruned_ott_id", ott7665092 = "pruned_ott_id", ott7665381 = "pruned_ott_id", ott7765249 = "pruned_ott_id", ott7981613 = "pruned_ott_id")
ids <- vector(mode="numeric", length(lista))
ids <- as.numeric(ids)
for(i in seq_along(lista)){
  ids[i] <- gsub("[a-z]", "", names(lista[i]))
}
pos <- which(species_info_ott_without_na$ott_id %in% ids)
species_info_ott_without_na <- species_info_ott_without_na[-pos, ]
species_tree_ott <- tol_induced_subtree(ott_ids = species_info_ott_without_na$ott_id)

species_info_ncbi <- classification(species_names, db='ncbi')
species_tree_ncbi <- class2tree(species_info_ncbi, check = TRUE)
tree_ncbi <- species_tree_ncbi$phylo
