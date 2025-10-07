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

# load plant coordinates
emp <- read.csv(file.path("CoordinateCleaner_cleanCoordinates.csv"))
emp <- emp[-which(emp$species == ""),]
families <- unique(emp$family)
class <- "plants"

# loop through families
for(family_j in 1:length(families)){
    
    family <- families[family_j]
    
    print(paste(family_j, family, sep="   :   "))

    species <-  unique(emp$species[which(emp$family %in% family)])
    order <- "Plantae"

    asia_data_family[indx_family, ]$family <- family
    asia_data_family[indx_family, ]$total_sp_family <- length(species)
    
    # get spatial data
    sp <- emp[which(emp$species %in% species), ]

    species_ranges <- c(unlist(sapply(unique(sp$species), FUN=function(x){vect(c(c(vect(sp[which(sp$species %in% x), ], geom= c("decimalLongitude", "decimalLatitude")))))})))
    
    india_species <- sapply(species_ranges, FUN=function(x){any(extract(india_ras,x )$Beck_KG_V1_present_0p5 %in% 1)})
    tropical_asian_species <- sapply(species_ranges, FUN=function(x){any(extract(tropical_asia_ras,x )$Beck_KG_V1_present_0p5 %in% 1)})
    
    ind_species_names  <- names(india_species)[which(india_species)]
    tropasia_species_names  <- names(tropical_asian_species)[which(tropical_asian_species)]
    
    asia_data_family[indx_family, ]$india_tot_family     <- length(ind_species_names)
    asia_data_family[indx_family, ]$tropasia_tot_family       <- length(tropasia_species_names)

    tropical_rainforest <- sapply(species_ranges, FUN=function(x){length(which(extract(koppen,x )$Beck_KG_V1_present_0p5 %in% 1)) >= 0.1*length(extract(koppen,x )$Beck_KG_V1_present_0p5)})
    tropical_seasonal <- sapply(species_ranges, FUN=function(x){length(which(extract(koppen,x )$Beck_KG_V1_present_0p5 %in% c(2,3))) >= 0.1*length(extract(koppen,x )$Beck_KG_V1_present_0p5)})
    montane <- sapply(species_ranges, FUN=function(x){length(which(extract(elevation,x )$srtm_1km >=1000)) >= 0.1*length(extract(elevation,x )$srtm_1km)})
    
    tropical_rainforest <- tropical_rainforest[names(tropical_rainforest) %in% tropasia_species_names]
    tropical_seasonal <- tropical_seasonal[names(tropical_seasonal) %in% tropasia_species_names]
    montane <- montane[names(montane) %in% tropasia_species_names]
    
    
    asia_data_family[indx_family, ]$tropical_rainforest_tot_family <- length(which(tropical_rainforest))
    asia_data_family[indx_family, ]$tropical_seasonal_tot_family   <- length(which(tropical_seasonal))
    asia_data_family[indx_family, ]$montane_tot_family   <- length(which(montane))
    
    genera <- na.omit(unique(sapply(species, FUN=function(x)strsplit(x, " ")[[1]][1])))
    indx_family <- indx_family +1
    
    for(genus_k in 1:length(genera)){
      genus <- genera[genus_k]
      species_in_genus <- species[grepl(genus, species)]
      
      asia_data[indx, ]$class <- class
      asia_data[indx, ]$order <- order
      asia_data[indx, ]$family <- family
      asia_data[indx, ]$genus <- genus
      
      # get total number of species
      asia_data[indx, ]$total_sp <- length(species_in_genus)
      
      # get spatial data
      sp <- emp[which(emp$species %in% species_in_genus), ]
      
      species_ranges <- c(unlist(sapply(unique(sp$species), FUN=function(x){vect(c(c(vect(sp[which(sp$species %in% x), ], geom= c("decimalLongitude", "decimalLatitude")))))})))
      
      india_species <- sapply(species_ranges, FUN=function(x){any(extract(india_ras,x )$Beck_KG_V1_present_0p5 %in% 1)})
      tropical_asian_species <- sapply(species_ranges, FUN=function(x){any(extract(tropical_asia_ras,x )$Beck_KG_V1_present_0p5 %in% 1)})
      
      ind_species_names  <- names(india_species)[which(india_species)]
      tropasia_species_names  <- names(tropical_asian_species)[which(tropical_asian_species)]
      
      
      asia_data[indx, ]$india_tot    <- length(ind_species_names)
      asia_data[indx, ]$tropasia_tot       <- length(tropasia_species_names)
      
      tropical_rainforest <- sapply(species_ranges, FUN=function(x){length(which(extract(koppen,x )$Beck_KG_V1_present_0p5 %in% 1)) >= 0.1*length(extract(koppen,x )$Beck_KG_V1_present_0p5)})
      tropical_seasonal <- sapply(species_ranges, FUN=function(x){length(which(extract(koppen,x )$Beck_KG_V1_present_0p5 %in% c(2,3))) >= 0.1*length(extract(koppen,x )$Beck_KG_V1_present_0p5)})
      montane <- sapply(species_ranges, FUN=function(x){length(which(extract(elevation,x )$srtm_1km >=1000)) >= 0.1*length(extract(elevation,x )$srtm_1km)})
      
      tropical_rainforest <- tropical_rainforest[names(tropical_rainforest) %in% tropasia_species_names]
      tropical_seasonal <- tropical_seasonal[names(tropical_seasonal) %in% tropasia_species_names]
      montane <- montane[names(montane) %in% tropasia_species_names]
      
      
      asia_data[indx, ]$tropical_rainforest_tot <- length(which(tropical_rainforest))
      asia_data[indx, ]$tropical_seasonal_tot   <- length(which(tropical_seasonal))
      asia_data[indx, ]$montane_tot   <- length(which(montane))
      indx <- indx + 1
      
    }
  }
saveRDS(asia_data, file.path(save_dir,"asia_plant_survey_realms_habitat.rds"))
saveRDS(asia_data_family, file.path(save_dir, "asia_plant_family_survey_realms_habitat.rds"))

plant_data <- merge(asia_data, asia_data_family, by="family", all=T)
saveRDS(plant_data, file.path(save_dir, "asia_plant_family_and_genus_survey_habitat.rds"))

