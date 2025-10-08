####################################################################################################
### Run intersections between mammal rangemaps (MDD v1.2) and bioregions
### Charlie Marsh
### charliem2003@github
### 03/2024
###
### Regions map used: Bioregions_rangemaps
### Taxonomy used:    MDD v1.2 (https://zenodo.org/records/4139818)
### Rangemaps used:   MDD v1.2 (exclusively marine species removed)
###
### CITATION: Mammal Diversity Database. (2020). Mammal Diversity Database (Version 1.2) [Data set].
###  Zenodo. http://doi.org/10.5281/zenodo.4139818.
###  Map of Life. (2021). Mammal range maps harmonised to the Mammals Diversity Database [Data set].
###  Map of Life. https://doi.org/10.48600/MOL-48VZ-P413
###
### saves csv with the range area (km2) occurring within each region/grid cell along with total
### range size (range_area)
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(sf)
library(dplyr)
library(tidyr)
library(units)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity_and_endemism"                        # project dir
rangeDir  <- file.path("MDD", "rangemaps", "directory")      # dir that contains the rangemaps data
taxonDir  <- file.path("MDD", "assessments", "directory")    # dir that contains the taxonomy data (for description year)
gadmDir   <- file.path("GADM", "directory")                  # dir that contains GADM in equal-area projection (called 'GADM_410_land_Equal_Area.gpkg')

### You shouldn't need to adjust these folders
regionDir <- file.path(projDir, "Data")                      # dir that contains the subregions data
funDir    <- file.path(projDir, "Analysis_functions")        # dir that contains the function scripts
resDir    <- file.path(projDir, "Intersections")             # dir to save results to
spDir     <- file.path(resDir,  "Intersections", "Mammals")  # dir to save species intersections to

### Create folders for storing species intersections
if(!dir.exists(resDir)) { dir.create(resDir, recursive = TRUE) }
if(!dir.exists(spDir))  { dir.create(spDir,  recursive = TRUE) }
if(!dir.exists(file.path(spDir, "Errors")))         { dir.create(file.path(spDir, "Errors")) }
if(!dir.exists(file.path(spDir, "Non-native")))     { dir.create(file.path(spDir, "Non-native")) }
if(!dir.exists(file.path(spDir, "Outside_extent"))) { dir.create(file.path(spDir, "Outside_extent")) }
if(!dir.exists(file.path(spDir, "Running")))        { dir.create(file.path(spDir, "Running")) }

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Intersections function
source(file.path(funDir, "Intersections", "intersections_rangemaps_parallel.R"))
source(file.path(funDir, "Intersections", "fast_intersections_helper_functions.R"))
source(file.path(funDir, "Intersections", "convertToPA.R"))

### Function for cleaning dates and extracting region info from inside brackets
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))

### regions data to calculate polygons for
regions <- st_read(file.path(regionDir, "Bioregions_rangemaps.gpkg"))

### global GADM for masking out dodgy coastlines in IUCN range maps on global range maps
gadm <- st_read(file.path(gadmDir, "GADM_410_land_Equal_Area.gpkg"))

### range maps
done <- c(list.files(spDir, include.dirs = FALSE, recursive = FALSE),
          list.files(file.path(spDir, "Outside_extent"), include.dirs = FALSE, recursive = FALSE),
          list.files(file.path(spDir, "Non-native"), include.dirs = FALSE, recursive = FALSE))
done <- gsub("_", " ", done)

# MDD maps are arranged by order
# - loop through orders (loading in range maps for all species is v. large)
orders <- list.files(file.path(rangeDir, "orders"))
for(order in orders) {
  rangemaps <- st_read(file.path(rangeDir, "orders", order))

  # remove strictly marine species
  marine <- read.csv(file.path(rangeDir, "marineSpList_mdd.csv")) %>%
    filter(other_freshwater == 0 & terrestrial == 0) %>%
    select(sciname)
  marine <- gsub("_", " ", marine$sciname)
  rangemaps <- filter(rangemaps, !sciname %in% marine)

  ### total species list extracted from range maps - 6238 species for full range maps
  spList <- unique(rangemaps$sciname)
  spList <- spList[!spList %in% gsub(".csv", "", done)]
  spList <- sample(spList, length(spList), replace = FALSE)
  print(length(spList))
  rangemaps <- filter(rangemaps, sciname %in% spList)

  # reproject to equal-area projection
  rangemaps <- rangemaps %>%
    st_transform(st_crs(regions))

  #==================================================================================================#
  #--------------------------------------- Run intersections ----------------------------------------#
  #==================================================================================================#
  
  ### NOTE - parallelisation only works on linux or mac. For windows set threads = 1
  intersections_rangemaps_parallel(rangemaps        = rangemaps,    # species rangemaps
                                   spNameCol        = "sciname",    # column name where species names are kept
                                   regions          = regions,      # shapefile with regions to summarise by
                                   regionNameCol    = "Bioregion",  # column name where region names are kept
                                   gadm             = gadm,         # global gadm for initial masking out where coastlines don't overlap
                                   crop_to_region   = TRUE,         # exclude species outside region
                                   areaThresh       = 1000000,     # area threshold in km2 for when to split the shapefile in to multiple parts
                                   spList           = spList,       # species list to subset rangemaps with
                                   spDir            = spDir,        # directory to save individual species intersections
                                   overwriteSpFiles = FALSE,        # overwrite existing species intersections, otherwise skips
                                   threads          = 5)            # no cores for parallelisation
}

#==================================================================================================#
#---------------------------------- Merge together species files ----------------------------------#
#==================================================================================================#

### read in all trop Asia files
spList <- list.files(spDir, recursive = FALSE, full.names = TRUE, include.dirs = FALSE, pattern = ".csv")
intersections <- data.frame()
for(i in 1:length(spList)) {
  sp <- read.csv(spList[i])
  intersections <- rbind(intersections, sp)
}

### pivot results data frame so each regional polygon is a column
intersections <- intersections %>%
  pivot_wider(names_from  = Region,
              values_from = Intersect_area,
              values_fn   = sum,
              values_fill = 0) %>%
  mutate(Species = gsub(" ", "_", Species))

### clean taxonomy data frame to get taxonomic info and year described
taxonomy <- read.csv(file.path(taxonDir, "MDD_v1.2_6485species.csv"))  %>%
  select(sciName, genus, family, authoritySpeciesYear) %>%
  rename(Species = sciName, Genus = genus, Family = family, Year = authoritySpeciesYear) %>%
  select(Species, Genus, Family, Year)

### their is a name discrepancy with the range maps in one species
taxonomy$Species[taxonomy$Species == "Paradoxurus_philippensis"] <- "Paradoxurus_philippinensis"

### merge intersections with taxonomy data frame
intersections <- intersections %>%
  left_join(taxonomy, by = "Species") %>%
  mutate(Genus = sapply(Species, function(x) strsplit(x, "_")[[1]][1])) %>%
  relocate(Species, Genus, Family, Year)

#==================================================================================================#
#---------------------------------- Convert to presence-absence -----------------------------------#
#==================================================================================================#

### 0.0% - 0%
pa0_0 <- convertToPA(rangeIntersects = intersections,    # data frame of occupied area for species range maps x bioregions
                     threshAsia      = NULL,             # prop of species rangemap that needs to occur within trop. asia to be considered
                     threshBioregion = 0,                # prop of bioregion area that needs to be occupied to be considered present in bioregion
                     threshPropRange = 0,                # if < threshBioregion, prop of range to occur within bioregion to still be considered present
                     threshAsiaEnd   = 0.999,             # prop. total species range that occurs in tropical asia to be considered regional endemic (Australia is removed)
                     rangeAreaCol    = "Range_area",     # col name that contains species total range areas
                     regionNameCol   = "Bioregion",   # column name where unique region names are kept
                     regions         = regions,          # shapefile with bioregions
                     verbose         = TRUE)
write.csv(pa0_0, file.path(resDir, "Intersections", "Intersections_bioregions_mammals_0_0.csv"),
          quote = FALSE, row.names = FALSE)

### 2.5% - 15%
pa2.5_15 <- convertToPA(rangeIntersects = intersections, # data frame of occupied area for species range maps x bioregions
                        threshAsia      = NULL,          # prop of species rangemap that needs to occur within trop. asia to be considered
                        threshBioregion = 0.025,         # prop of bioregion area that needs to be occupied to be considered present in bioregion
                        threshPropRange = 0.15,          # if < threshBioregion, prop of range to occur within bioregion to still be considered present
                        threshAsiaEnd   = 0.90,          # prop. total species range that occurs in tropical asia to be considered regional endemic (Australia is removed)
                        rangeAreaCol    = "Range_area",  # col name that contains species total range areas
                        regionNameCol   = "Bioregion",   # column name where unique region names are kept
                        regions         = regions,       # shapefile with bioregions
                        verbose         = TRUE)
write.csv(pa2.5_15, file.path(resDir, "Intersections", "Intersections_bioregions_mammals_2.5_15.csv"),
          quote = FALSE, row.names = FALSE)

### 5% - 25%
pa5.0_25 <- convertToPA(rangeIntersects = intersections, # data frame of occupied area for species range maps x bioregions
                        threshAsia      = NULL,          # prop of species rangemap that needs to occur within trop. asia to be considered
                        threshBioregion = 0.050,         # prop of bioregion area that needs to be occupied to be considered present in bioregion
                        threshPropRange = 0.25,          # if < threshBioregion, prop of range to occur within bioregion to still be considered present
                        threshAsiaEnd   = 0.95,          # prop. total species range that occurs in tropical asia to be considered regional endemic (Australia is removed)
                        rangeAreaCol    = "Range_area",  # col name that contains species total range areas
                        regionNameCol   = "Bioregion",   # column name where unique region names are kept
                        regions         = regions,       # shapefile with bioregions
                        verbose         = TRUE)
write.csv(pa5.0_25, file.path(resDir, "Intersections", "Intersections_bioregions_mammals_5.0_25.csv"),
          quote = FALSE, row.names = FALSE)

### 7.5% - 35%
pa7.5_35 <- convertToPA(rangeIntersects = intersections, # data frame of occupied area for species range maps x bioregions
                        threshAsia      = NULL,          # prop of species rangemap that needs to occur within trop. asia to be considered
                        threshBioregion = 0.075,         # prop of bioregion area that needs to be occupied to be considered present in bioregion
                        threshPropRange = 0.35,          # if < threshBioregion, prop of range to occur within bioregion to still be considered present
                        threshAsiaEnd   = 0.99,          # prop. total species range that occurs in tropical asia to be considered regional endemic (Australia is removed)
                        rangeAreaCol    = "Range_area",  # col name that contains species total range areas
                        regions         = regions,       # shapefile with bioregions
                        verbose         = TRUE)
write.csv(pa7.5_35, file.path(resDir, "Intersections", "Intersections_bioregions_mammals_7.5_35.csv"),
          quote = FALSE, row.names = FALSE)

### Overall species list
pa0_0$Species[!pa0_0$Species %in% pa2.5_15$Species]
pa2.5_15$Species[!pa2.5_15$Species %in% pa5.0_25$Species]
pa5.0_25$Species[!pa5.0_25$Species %in% pa7.5_35$Species]

allSp <- select(pa0_0, Species, Genus, Family, Year, Range_area, AsiaEndemic) %>%
  mutate(thresh0_0 = 1) %>%
  rename(AsiaEndemic_0.999 = AsiaEndemic) %>%
  left_join(select(pa2.5_15, Species, AsiaEndemic) %>% mutate(thresh2.5_15 = 1), by = "Species") %>%
  rename(AsiaEndemic_0.90 = AsiaEndemic) %>%
  left_join(select(pa5.0_25, Species, AsiaEndemic) %>% mutate(thresh5.0_25 = 1), by = "Species") %>%
  rename(AsiaEndemic_0.95 = AsiaEndemic) %>%
  left_join(select(pa7.5_35, Species, AsiaEndemic) %>% mutate(thresh7.5_35 = 1), by = "Species") %>%
  rename(AsiaEndemic_0.99 = AsiaEndemic) %>%
  mutate(across(all_of(c("thresh0_0", "thresh2.5_15", "thresh5.0_25", "thresh7.5_35")),
                ~ case_when(. == 1   ~ 1,
                            is.na(.) ~ 0))) %>%
  relocate(c(thresh0_0, thresh2.5_15, thresh5.0_25, thresh7.5_35), .before = AsiaEndemic_0.999)
write.csv(allSp, file.path(resDir, "Intersections", "Intersections_bioregions_mammals_spList.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#--------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

### removing rangemaps doesn't clear up the mem used.
### We need to restart the session (at least in RStudio and Linux)
rm(rangemaps)
.rs.restartR()
