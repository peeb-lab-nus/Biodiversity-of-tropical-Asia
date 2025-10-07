######################################
###            METADATA            ###
######################################
# Version: 4.0
#
# Author: Alexander Skeels
#
# Date: 4.07.2024
#
# Projects: BTAS
#
# Contact: alexander.skeels@gmail.com
#
######################################

#libraries
library(fishtree)
library(rfishbase)

# directory
rabosky_dir <- "doi_10_5061_dryad_fc71cp4__v20190529/Archive"

#trees
trees <- fishtree_complete_phylogeny()

#spatial <- read.csv(file.path(rabosky_dir,"dataFiles", "species_by_cell_PAmat_merged.csv"))
meow <- read.csv(file.path(rabosky_dir,"dataFiles", "species_endemicity_MEOW.csv"))
tropasaian_species <- meow$X[which(meow$Central_Indo.Pacific > 1)]
tropasian_genera <- unique(sapply(tropasaian_species, FUN=function(x)strsplit(x, "_")[[1]][1]))

# set up data storage
asia_data <- data.frame(genus=tropasian_genera, 
                        total_sp=NA, tropasia_tot=NA, 
                        crown_age_median =NA, crown_age_max =NA, crown_age_min =NA,
                        stem_age_median =NA, stem_age_max =NA, stem_age_min =NA)



all_tips <- meow$X[grepl(paste(tropasian_genera, collapse="|"), meow$X)]
trees <- lapply(trees, FUN=function(x)keep.tip(x, x$tip.label[which(x$tip.label %in% all_tips)]))

# loop over each genus to get crown and stem ages
for(i in 1:length(tropasian_genera)){
  
  genus_i <- tropasian_genera[i]
  print(paste(i, genus_i ))
  
  tropasian_species_i <- tropasaian_species[grepl( genus_i, tropasaian_species)]
  allspecies_i <- meow$X[grepl( genus_i, meow$X)]
  allspecies_i <- allspecies_i[which(allspecies_i %in% trees[[1]]$tip.label)]
  if(!any(trees[[1]]$tip.label %in%allspecies_i)){next}
  
  prop_tropasian <-  length(tropasian_species_i) /  length(allspecies_i)
  crown_ages_i <- vector("numeric", 100)
  stem_ages_i <- vector("numeric", 100)
  
  if(length(allspecies_i) <2){
    trees_i <- lapply(trees, FUN=function(x) {
      return(keep.tip(x, allspecies_i))}
    )
  } else{
    trees_i <- lapply(trees, FUN=function(x) {
      crown_node <- getMRCA(x,allspecies_i[which(allspecies_i %in% x$tip.label)]);
      sister_node <-  Siblings(x, crown_node);
      stem_node <- getMRCA(x, c(crown_node, sister_node));
      return(extract.clade(x,  stem_node))}
    )
    
  }
  
  for(j in 1:100){
    tree_j <- trees_i[[j]]
    tree_j <- drop.tip(tree_j, tree_j$tip.label[which(!tree_j$tip.label %in% allspecies_i)])
    
    crown_ages_i[j] <- max(nodeHeights(tree_j))
    if(length(tree_j$tip.label)==1){stem_ages_i[j] <- crown_ages_i[j]; crown_ages_i[j] <- NA; next}
    
    crown_node <- getMRCA(trees_i[[j]],allspecies_i[which(allspecies_i %in% trees_i[[j]]$tip.label)])
    sister_node <-  Siblings(trees_i[[j]], crown_node)
    stem_node <- getMRCA(trees_i[[j]], c(crown_node, sister_node))

    total_tree_age <- max(c(nodeHeights(trees_i[[j]])))
    stem_ages_i[j] <-total_tree_age- nodeheight(trees_i[[j]], stem_node)
  }
  asia_data$total_sp[i] <-  length(allspecies_i)
  asia_data$tropasia_tot[i] <-  length(tropasian_species_i)
  asia_data$crown_age_median[i]   <- median(crown_ages_i, na.rm=T)
  asia_data$crown_age_max[i]    <- max(crown_ages_i, na.rm=T)
  asia_data$crown_age_min[i]      <- min(crown_ages_i, na.rm=T)
  asia_data$stem_age_median[i] <- median(stem_ages_i, na.rm=T)
  asia_data$stem_age_max[i] <- max(stem_ages_i, na.rm=T)
  asia_data$stem_age_min[i] <- min(stem_ages_i, na.rm=T)
  
}

write.csv(asia_data, paste0("fish_crown_ages_complete.csv"))

# now get habitat association data for each clade from fishbase
fishb <- read.csv(file.path(rabosky_dir,"dataFiles", "taxonConversion.csv"))
species_list <- fishb$fishbase[which(fishb$treeTaxon %in% all_tips)]
all_dat <- species_list(Genus=tropasian_genera)
all_ecology <- ecology(species_list = all_dat )

reef_associated <- all_ecology$Species[which(all_ecology$CoralReefs %in% -1)]
pelagic <- all_ecology$Species[which(  all_ecology$Pelagic %in% -1 | 
                                       all_ecology$Mesopelagic == -1 | 
                                       all_ecology$Bathypelagic == -1 |
                                       all_ecology$Hadopelagic == -1 |
                                       all_ecology$Epipelagic == -1|
                                       all_ecology$Abyssopelagic == -1)]

asia_data$n_reef_associated=NA
asia_data$n_pelagic=NA

for(i in 1:length(tropasian_genera)){
  
  genus_i <- tropasian_genera[i]
  print(paste(i, genus_i ))
  
  tropasian_species_i <- tropasaian_species[grepl( genus_i, tropasaian_species)]
  allspecies_i <- meow$X[grepl( genus_i, meow$X)]
  allspecies_i <- allspecies_i[which(allspecies_i %in% trees[[1]]$tip.label)]
  asia_data$n_reef_associated[i] <-  length(which(allspecies_i %in% gsub(" ", "_", reef_associated)))
  asia_data$n_pelagic[i] <-  length(which(allspecies_i %in% gsub(" ", "_", pelagic)))
  
}

write.csv(asia_data, "fish_crown_ages_complete_ecology.csv", row.names = F)
