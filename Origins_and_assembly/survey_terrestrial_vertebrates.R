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

#functions
endemic <- function(data, region, regions, endemic_prop){
  
  other_regions <- regions[which(!regions %in% region)]
  
  ifelse(data[,region] >= (endemic_prop *data$total_sp) &
           rowSums(data[, other_regions], na.rm=T) == 0, TRUE, FALSE)
  
}

#libraries
library(sp)
library(phytools)
library(raster)
library(ape)
library(terra)
library(ggplot2)
library(phangorn)

# get a detailed map of the region using Botanical Countries
regions <- vect("../wgsrpd-master/wgsrpd-master/level2/level2.shp") 
regions <- regions[which(regions$LEVEL2_NAM == "Malesia" | regions$LEVEL2_NAM == "Indo-China")]

# load in the key tectonic plates and give a CRS (Bird, P. (2003), An updated digital model of plate boundaries, Geochem. Geophys. Geosyst., 4, 1027, doi:10.1029/2001GC000252, 3. )
setwd(file.path("plate_boundaries", "tectonicplates-master"))

# read in shape files of all plates
plates <- vect( "PB2002_plates.shp")

# create new vector of Major plate names
India       <- plates[which(plates$PlateName == "India")]
TropicalAsia <- regions

# set CRS
crs(India ) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
crs(TropicalAsia) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

# simplify poolygons (speeds up calculations)
India <- buffer(simplifyGeom(India, tol=0.2), width=100000)
TropicalAsia <- buffer(simplifyGeom(TropicalAsia, tol=0.2), width=100000)

# set CRS
crs(India ) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
crs(TropicalAsia) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

# load koppen bands raster (Beck, H., Zimmermann, N., McVicar, T. et al. Present and future Köppen-Geiger climate classification maps at 1-km resolution. Sci Data 5, 180214 (2018). https://doi.org/10.1038/sdata.2018.214)
koppen <- rast("Beck_KG_v1_present_0p5.tif")

# load elevation raster (STRM; https://www.earthdata.nasa.gov)
elevation <- rast("elevation.tif")

# get regional masks for rasters
india_ras <- mask(koppen, India )
values(india_ras)[which(!is.na(values(india_ras)))] <- 1
values(india_ras)[which(!is.na(values(india_ras)))] <- 1

tropical_asia_ras <- mask(koppen, TropicalAsia )
values(tropical_asia_ras)[which(!is.na(values(tropical_asia_ras)))] <- 1

# set up data storage
asia_data <- data.frame(class=NA, order=NA, family=NA, genus=NA, 
                        total_sp=NA,
                        tropical_seasonal_tot=NA,tropical_rainforest_tot=NA,
                        india_tot =NA, tropasia_tot=NA, 
                        montane_tot=NA,
                        crown_age_median =NA, crown_age_max =NA, crown_age_min =NA,
                        stem_age_median =NA, stem_age_max =NA, stem_age_min =NA)


asia_data_family <- data.frame(family=NA,
                        total_sp_family=NA, 
                        india_tot_family =NA,  tropasia_tot_family=NA, 
                        tropical_seasonal_tot_family=NA,tropical_rainforest_tot_family=NA,
                        montane_tot_family=NA,
                        crown_age_family_median =NA, crown_age_family_max =NA, crown_age_family_min =NA,
                        stem_age_family_median =NA, stem_age_family_max =NA, stem_age_family_min =NA)



# loop over classes
indx <- 1
indx_family <- 1

for(class_i in 1:5){
  
  # load data
  setwd(data_dir)
  if(class_i ==1) {emp <- readRDS("cleaned_mammal_data.rds"); families <- unique(emp$taxonomy$fam); class <- "mammals"}
  if(class_i ==2) {emp <- readRDS("cleaned_bird_data.rds"); families <- unique(emp$taxonomy$BLFamilyLatin); class <- "birds"}
  if(class_i ==3) {emp <- readRDS("cleaned_amphibian_data.rds"); families <- unique(emp$taxonomy$Family); class <- "amphibians"}
  if(class_i ==4) {emp <- readRDS("cleaned_snake_data.rds");  families <- unique(emp$traits$family); class <- "reptiles"}
  if(class_i ==5) {emp <- readRDS("cleaned_lizard_data.rds"); families <- unique(emp$traits$Family); class <- "reptiles"}


  # loop through families1:length(families)
  for(family_j in 1:length(families)){
    
    family <- families[family_j]
  
     print(paste(family_j, family, sep="   :   "))
     
     # slightly different operation depednding on data source
     if(class_i ==1){
       species <-  unique(emp$taxonomy$Species_Name[which(emp$taxonomy$fam %in% family)])
       order <- unique(emp$taxonomy$ord[which(emp$taxonomy$fam %in% family)])
     }
     if(class_i ==2){
       species <-  unique(emp$taxonomy$Scientific[which(emp$taxonomy$BLFamilyLatin %in% family)])
       order <- unique(emp$taxonomy$IOCOrder[which(emp$taxonomy$BLFamilyLatin %in% family)])
     }
     if(class_i ==3){
       species <-  unique(emp$taxonomy$Scientific.Name[which(emp$taxonomy$Family %in% family)])
       species <- gsub(" ", "_", species)
       order <- unique(emp$taxonomy$Taxon[which(emp$taxonomy$Family %in% family)])
       
     }
     if(class_i ==4){
       species <-  unique(emp$traits$species[which(emp$traits$family %in% family)])
       order <- "Squamata"
     }
     if(class_i ==5){
       species <-  unique(emp$traits$species[which(emp$traits$Family %in% family)])
       order <- "Squamata"
     }

     asia_data_family[indx_family, ]$family <- family
     asia_data_family[indx_family, ]$total_sp_family <- length(species)
     
     # subset relavant phylogeny
     phy <- drop.tip(emp$phy[[1]], emp$phy[[1]]$tip.label[which(!emp$phy[[1]]$tip.label %in% species)])
     # if its empty skip
     if(is.null(phy)){indx_family <- indx_family + 1; indx <- indx+1; next}
     
     # get crown and stem ages from 100 different phylogenies from posterior
     crown_age <- c()
     stem_age <- c()
     for(phy_l in 1:100){
       phy <- drop.tip(emp$phy[[phy_l]], emp$phy[[phy_l]]$tip.label[which(!emp$phy[[phy_l]]$tip.label %in% species)])
       crown_age[phy_l] <- max(c(nodeHeights(phy)))
       
       if(length(phy$tip.label)==1){stem_age[phy_l] <- crown_age[phy_l]; next}
       crown_node <- getMRCA(emp$phy[[phy_l]], species[which(species %in% emp$phy[[phy_l]]$tip.label)])
       sister_node <-  Siblings(emp$phy[[phy_l]], crown_node)
       stem_node <- getMRCA(emp$phy[[phy_l]], c(crown_node, sister_node))
       
       total_tree_age <- max(c(nodeHeights(emp$phy[[phy_l]])))
       stem_age[phy_l] <-total_tree_age- nodeheight(emp$phy[[phy_l]], stem_node)
       
     }
     
     # save crown and stem ages for families
     asia_data_family[indx_family, ]$crown_age_family_median <- median(crown_age, na.rm=T)
     asia_data_family[indx_family, ]$crown_age_family_max <- max(crown_age, na.rm=T)
     asia_data_family[indx_family, ]$crown_age_family_min <- min(crown_age, na.rm=T)
     asia_data_family[indx_family, ]$stem_age_family_median <- median(stem_age, na.rm=T)
     asia_data_family[indx_family, ]$stem_age_family_max <- max(stem_age, na.rm=T)
     asia_data_family[indx_family, ]$stem_age_family_min <- min(stem_age, na.rm=T)
     
     # get spatial data for family
     sp <- emp$sp[, c(1:2, which(colnames(emp$sp) %in% phy$tip.label))]
     
     # subset
     sp <- sp[which(rowSums(sp[, 3:ncol(sp), drop=F]) > 0),, drop=F ]
     coords_sp <- sp[, 1:2, drop=F]
     sp <- sp[, 3:ncol(sp), drop=F]
     if(class_i %in% c(4,5)){
       coords_sp <- SpatialPoints(coords_sp)
       crs(coords_sp) <-      '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs'
       coords_sp<- spTransform(coords_sp, "+proj=longlat +datum=WGS84 +no_defs +type=crs")@coords
       
     }
     
     #coolies
     sp_df <- SpatialPointsDataFrame(as.data.frame(coords_sp), as.data.frame(sp))
     
     crs(sp_df) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
     sp_df_v <- vect(sp_df)
     ind_sp  <- data.frame(intersect(India,sp_df_v))
     tropasia_sp  <- data.frame(intersect(TropicalAsia,sp_df_v))
     tropasia_sp  <- tropasia_sp[,-which(colnames(tropasia_sp) %in% c("LEVEL2_COD",
                                                                      "LEVEL1_COD",
                                                                      "LEVEL1_NAM",
                                                                      "LEVEL2_NAM")), drop=FALSE]
     
     # get species names found in each region
     ind_species_names  <- which(colSums(ind_sp[,4:ncol(ind_sp), drop=F]) > 0)
     tropasia_species_names  <- which(colSums(tropasia_sp) > 0)
     
     asia_data_family[indx_family, ]$india_tot_family    <- length(ind_species_names)
     asia_data_family[indx_family, ]$tropasia_tot_family       <- length(tropasia_species_names)
     
     species_koppen <- extract(koppen, sp_df_v)
     values(elevation)[which(is.na(values(elevation)))] <- 0
     species_elevation <- extract(elevation, sp_df_v)
     
     species_koppen_list <- sapply(sp_df@data,  FUN=function(x){table(species_koppen[which(x==1), "Beck_KG_V1_present_0p5"])})
     species_elevation_list <- sapply(sp_df@data,  FUN=function(x){species_elevation[which(x==1), "srtm_1km"]})
     
     #if 10% of a species range is in an environment, lets class it as present in that environment
     tropical_rainforest <- species_koppen_list[sapply(species_koppen_list, FUN=function(x){if(any(names(x) %in% c( "1"))){
       sum(x[which(names(x) %in% c( "1"))]) >= (sum(x)*0.1)
     } else{ FALSE }
     })]
     
     tropical_seasonal <- species_koppen_list[sapply(species_koppen_list, FUN=function(x){if(any(names(x) %in% c("2", "3"))){
       sum(x[which(names(x) %in% c("2","3"))]) >= (sum(x)*0.1)
     } else{ FALSE }
     })]
     
     montane <- species_elevation_list[sapply(species_elevation_list, FUN=function(x){x[which(is.na(x))] <- 0; if(any(x >1000)){
       length(x[which(x > 1000)]) >= (length(x)*0.1)
     } else{ FALSE }
     })]
     
     tropical_rainforest <- tropical_rainforest[names(tropical_rainforest) %in% names(tropasia_species_names)]
     tropical_seasonal <- tropical_seasonal[names(tropical_seasonal) %in% names(tropasia_species_names)]
     montane <- montane[names(montane) %in% names(tropasia_species_names)]

     asia_data_family[indx_family, ]$tropical_rainforest_tot_family <- length(tropical_rainforest)
     asia_data_family[indx_family, ]$tropical_seasonal_tot_family   <- length(tropical_seasonal)
     asia_data_family[indx_family, ]$montane_tot_family   <- length(montane)
     
     genera <- unique(sapply(species, FUN=function(x)strsplit(x, "_")[[1]][1]))
     indx_family <- indx_family +1
     
     # now for each family, loop over each genus and do the same as before
     for(genus_k in 1:length(genera)){
       genus <- genera[genus_k]
       species_in_genus <- species[grepl(genus, species)]
       
       asia_data[indx, ]$class <- class
       asia_data[indx, ]$order <- order
       asia_data[indx, ]$family <- family
       asia_data[indx, ]$genus <- genus
       
       # get total number of species
       asia_data[indx, ]$total_sp <- length(species_in_genus)
       
       # get phylo
       phy <- drop.tip(emp$phy[[1]], emp$phy[[1]]$tip.label[which(!emp$phy[[1]]$tip.label %in% species_in_genus)])
       
       # if its empty skip
       if(is.null(phy)){indx <- indx + 1; next}
       
       crown_age <- c()
       stem_age <- c()
       for(phy_l in 1:100){
         phy <- drop.tip(emp$phy[[phy_l]], emp$phy[[phy_l]]$tip.label[which(!emp$phy[[phy_l]]$tip.label %in% species_in_genus)])
         
         crown_age[phy_l] <- max(c(nodeHeights(phy)))
         if(length(phy$tip.label)==1){stem_age[phy_l] <- crown_age[phy_l]; next}
         
         crown_node <- getMRCA(emp$phy[[phy_l]], species_in_genus[which(species_in_genus %in% emp$phy[[phy_l]]$tip.label)])
         sister_node <-  Siblings(emp$phy[[phy_l]], crown_node)
         stem_node <- getMRCA(emp$phy[[phy_l]], c(crown_node, sister_node))
         
         total_tree_age <- max(c(nodeHeights(emp$phy[[phy_l]])))
         stem_age[phy_l] <-total_tree_age- nodeheight(emp$phy[[phy_l]], stem_node)

       }
  
       asia_data[indx, ]$crown_age_median <- median(crown_age)
       asia_data[indx, ]$crown_age_max <- max(crown_age)
       asia_data[indx, ]$crown_age_min <- min(crown_age)
       asia_data[indx, ]$stem_age_median <- median(stem_age)
       asia_data[indx, ]$stem_age_max <- max(stem_age)
       asia_data[indx, ]$stem_age_min <- min(stem_age)
       
       # get spatial data
       sp <- emp$sp[, c(1:2, which(colnames(emp$sp) %in% phy$tip.label))]
       
       # subset
       sp <- sp[which(rowSums(sp[, 3:ncol(sp), drop=F]) > 0),, drop=F ]
       
       # fix up
       coords_sp <- sp[, 1:2, drop=F]
       sp <- sp[, 3:ncol(sp), drop=F]
       if(class_i %in% c(4,5)){
         coords_sp <- SpatialPoints(coords_sp)
         crs(coords_sp) <-      '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs'
         coords_sp<- spTransform(coords_sp, "+proj=longlat +datum=WGS84 +no_defs +type=crs")@coords
         
       }
       #coolies
       sp_df <- SpatialPointsDataFrame(as.data.frame(coords_sp), as.data.frame(sp))
       
       crs(sp_df) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
       sp_df_v <- vect(sp_df)
       ind_sp  <- data.frame(intersect(India,sp_df_v))[,,drop=F]
       tropasia_sp  <- data.frame(intersect(TropicalAsia,sp_df_v))
       tropasia_sp  <- tropasia_sp[,-which(colnames(tropasia_sp) %in% c("LEVEL2_COD",
                                                                        "LEVEL1_COD",
                                                                        "LEVEL1_NAM",
                                                                        "LEVEL2_NAM")), drop=FALSE]
       


       
       # get species names found in each region
       ind_species_names  <- which(colSums(ind_sp[,4:ncol(ind_sp), drop=F]) > 0)
       tropasia_species_names  <- which(colSums(tropasia_sp) > 0)
       
       
       asia_data[indx, ]$india_tot    <- length(ind_species_names)
       asia_data[indx, ]$tropasia_tot <- length( tropasia_species_names)
       
       asia_data[indx, ]$tropical_rainforest_tot <- length(tropical_rainforest[grepl(genus, names(tropical_rainforest))])
       asia_data[indx, ]$tropical_seasonal_tot   <- length(tropical_seasonal[grepl(genus, names(tropical_seasonal))])
       asia_data[indx, ]$montane_tot   <- length(montane[grepl(genus, names(montane))])
       
       
       indx <- indx + 1
       
     }
  }
  saveRDS(asia_data, file.path(save_dir,"asia_vertebrate_survey_realms_habitat.rds"))
  saveRDS(asia_data_family, file.path(save_dir, "asia_vertebrate_family_survey_realms_habitat.rds"))
}

# merge genus and family level data
vertebrate_data <- merge(asia_data, asia_data_family, by="family", all=T)

saveRDS(vertebrate_data, file.path(save_dir, "asia_vertebrate_family_and_genus_survey_realms_habitat.rds"))

# data and parameters
endemic_prop<- 0.7

# get which are endemic
vertebrate_data$india_endemic <- ifelse(vertebrate_data$india_tot >= (endemic_prop *vertebrate_data$total_sp), TRUE, FALSE)
vertebrate_data$tropasia_endemic <- ifelse(vertebrate_data$tropasia_tot >= (endemic_prop *vertebrate_data$total_sp), TRUE, FALSE)
vertebrate_data$troprain_endemic <- ifelse(vertebrate_data$tropical_rainforest_tot >= (endemic_prop *vertebrate_data$tropasia_tot), TRUE, FALSE)
vertebrate_data$tropseas_endemic <- ifelse(vertebrate_data$tropical_seasonal_tot >= (endemic_prop *vertebrate_data$tropasia_tot), TRUE, FALSE)
vertebrate_data$montane_endemic <- ifelse(vertebrate_data$montane_tot >= (endemic_prop *vertebrate_data$tropasia_tot), TRUE, FALSE)

# subset only asian clades
asia_data <- vertebrate_data[which(vertebrate_data$tropasia_endemic==1),]

# get main region
asia_data$region <- "TropicalAsia"
asia_data$region[which(asia_data$india_endemic)] <- "India"

write.csv(asia_data, "tropical_asian_vertebrate_crownstem_ages_habitat.csv")