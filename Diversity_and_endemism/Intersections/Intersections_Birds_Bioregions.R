####################################################################################################
### Run intersections between bird rangemaps (BirdLife International) and bioregions
### Charlie Marsh
### charliem2003@github
### 04/2024
###
### Regions map used: Bioregions_rangemaps
### Taxonomy used:    HBW_BirdLife List of Birds v.81 - downloaded 06/12/2023
### Rangemaps used:   BirdLife international 2023.1 - downloaded 06/12/2023
###
### CITATION: BirdLife International and Handbook of the Birds of the World (2023) Bird species
###   distribution maps of the world. Version 2023.1. Available at http://datazone.birdlife.org/species/requestdis.
### HBW and BirdLife International (2023) Handbook of the Birds of the World and BirdLife International
###   digital checklist of the birds of the world. Version 7. Available at:
###   http://datazone.birdlife.org/userfiles/file/Species/Taxonomy/HBW-BirdLife_Checklist_v8_Dec23.zip
###
### saves csv with the range area (km2) occurring within each region/grid cell along with total
### range size (range_area)
####################################################################################################

### NOTE: BirdLife Interational has maps missing for the following sensitive species. We have
###       included them in this analysis, but you must request them separately
# Cacatua sulphurea – Yellow-crested Cockatoo (native distribution only)
# Chloropsis cochinchinensis – Javan Leafbird
# Cissa thalassina – Javan Green-magpie
# Cyornis banyumas – Javan Blue-flycatcher
# Garrulax bicolor – Sumatran Laughingthrush
# Gracula robusta - Nias Hill Myna
# Leucopsar rothschildi– Bali Myna
# Padda oryzivora – Java Sparrow (native distribution only)
# Pterorhinus courtoisi– Blue-crowned Laughingthrush
# Pycnonotus bimaculatus – Orange-spotted Bulbul
# Pycnonotus zeylanicus – Straw-headed Bulbul
# Rhinoplax vigil – Helmeted Hornbill
# Tanygnathus everetti - Blue-backed Parrot
# Thalasseus bernsteini – Chinese Crested Tern (non-breeding distribution only)
# Zosterops flavus – Javan White-eye

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
projDir   <- "Diversity_and_endemism"                         # project dir
rangeDir  <- file.path("BirdLife", "rangemaps", "directory")  # dir that contains the rangemaps data
taxonDir  <- file.path("BirdLife", "taxonomy", "directory")   # dir that contains the taxonomy data ('HBW_BirdLife List of Birds v.81.csv')
gadmDir   <- file.path("GADM", "directory")                   # dir that contains GADM in equal-area projection (called 'GADM_410_land_Equal_Area.gpkg')

### You shouldn't need to adjust these folders
regionDir <- file.path(projDir, "Data")                       # dir that contains the subregions data
funDir    <- file.path(projDir, "Analysis_functions")         # dir that contains the function scripts
resDir    <- file.path(projDir, "Intersections")              # dir to save results to
spDir     <- file.path(resDir,  "Intersections", "Birds")     # dir to save species intersections to

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
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))

### regions data to calculate polygons for
regions <- st_read(file.path(regionDir, "Bioregions_rangemaps.gpkg"))

### global GADM for masking out dodgy coastlines in IUCN range maps on global range maps
gadm <- st_read(file.path(gadmDir, "GADM_410_land_Equal_Area.gpkg"))

### range maps - 17,522 polygons, 11,184 species
rangemaps <- st_read(file.path(rangeDir, "BOTW.gdb"))

### add in sensitive species and append to range maps - 17,549 polygons, 11,196 species
sensitiveSp <- st_read(file.path(rangeDir, "Sensitive_species", "SppDataRequest.shp"))
names(sensitiveSp) <- tolower(names(sensitiveSp))
sensitiveSp <- rename(sensitiveSp, Shape = "geometry")
rangemaps <- bind_rows(rangemaps, sensitiveSp)

### filter out non-relevant range polygons - 13,640 polygons, 11,009 species
rangemaps <- rangemaps %>%
  filter(!presence %in% c(4, 5, 6)) %>% # 1 = Extant; 2 = Prob. extant; 3 = Poss. extant; 4 = Poss. extinct; 5 = Extinct; 6 = Presence uncertain
  filter(!origin   %in% c(3, 4, 6)) %>% # 1 = Native; 2 = Reintroduced; 3 = Introduced; 4 = Vagrant; 5 = Origin uncertain; 6 = Assisted colonisation
  filter(!seasonal %in% c(3, 4))        # 1 = Resident; 2 = Breeding season; 3 = Non-breeding season; 4 = Passage; 5 = Seasonal occurrence uncertain

### add in taxonomy details to get family info
taxonomy <- read.csv(file.path(taxonDir, "HBW_BirdLife List of Birds v.81.csv"))
rangemaps <- rangemaps %>%
  left_join(select(taxonomy, Scientific.name, Family.name), join_by("sci_name" == "Scientific.name"))

### filter out pelagics - 13,226 polygons, 10,798 species
rangemaps <- rangemaps %>%
  filter(!Family.name %in% c("Phaethontidae", "Spheniscidae", "Oceanitidae", "Hydrobatidae", "Diomedeidae",
                             "Procellariidae", "Fregatidae", "Sulidae", "Stercorariidae", "Alcidae"))

### total species list extracted from range maps - 10,798 species
spList <- unique(rangemaps$sci_name)

done <- c(list.files(spDir),
          list.files(file.path(spDir, "Outside_extent")),
          list.files(file.path(spDir, "Non-native")))
done <- gsub("_", " ", done)
spList <- spList[!spList %in% gsub(".csv", "", done)]
length(spList)

### reproject to equal-area projection
rangemaps <- rangemaps %>%
  st_transform(st_crs(regions))

#==================================================================================================#
#--------------------------------------- Run intersections ----------------------------------------#
#==================================================================================================#

### NOTE - parallelisation only works on linux or mac. For windows set threads = 1
intersections_rangemaps_parallel(rangemaps        = rangemaps,    # species rangemaps
                                 spNameCol        = "sci_name",   # column name where species names are kept
                                 regions          = regions,      # shapefile with regions to summarise by
                                 regionNameCol    = "Bioregion",  # column name where region names are kept
                                 gadm             = gadm,         # global gadm for initial masking out where coastlines don't overlap
                                 crop_to_region   = TRUE,         # exclude species outside region
                                 areaThresh       = 10000000,     # area threshold in km2 for when to split the shapefile in to multiple parts
                                 spList           = spList,       # species list to subset rangemaps with
                                 spDir            = spDir,        # directory to save individual species intersections
                                 overwriteSpFiles = FALSE,        # overwrite existing species intersections, otherwise skips
                                 threads          = 3)            # no cores for parallelisation

### One species, Rhyticeros narcondami (Narcondam Hornbill), the range map is shifted off the island
### so that no part of the range overlaps land. Input it manually
sp <- "Rhyticeros narcondami"
fileName <- file.path(spDir, "Non-native", paste0(gsub(" ", "_", sp, ".csv")))
if(file.exists(fileName)) { file.remove(fileName) }
spMap <- filter(rangemaps, sci_name == sp)
overlap <- data.frame(Species        = gsub(" ", "_", sp),
                      Region         = "Andamans",
                      Range_area     = round(drop_units(sum(set_units(st_area(spMap), km^2))), 4),
                      Intersect_area = round(drop_units(sum(set_units(st_area(spMap), km^2))), 4))
write.csv(overlap, fileName, quote = FALSE, row.names = FALSE)

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
  group_by(Species) %>%
  pivot_wider(names_from  = Region,
              values_from = Intersect_area,
              values_fn   = sum,
              values_fill = 0) %>%
  mutate(Species = gsub(" ", "_", Species))

### clean taxonomy data frame to get higher taxonomic info
taxonomy <- read.csv(file.path(taxonDir, "HBW_BirdLife List of Birds v.81.csv"))  %>%
  select(Scientific.name, Family.name) %>%
  rename(Species = Scientific.name, Family = Family.name) %>%
  mutate(Species = gsub(" ", "_", Species))

### extract year described from shapefiles
year <- st_read(file.path(rangeDir, "BOTW.gdb"), layer = "Checklist_v8_txt") %>%
  mutate(Year = clean_publication_dates(Authority)) %>%
  mutate(Species = gsub(" ", "_", ScientificName)) %>%
  select(Species, Year)

### merge intersections with taxonomy data frame
intersections <- intersections %>%
  left_join(taxonomy, by = "Species") %>%
  left_join(year, by = "Species") %>%
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
write.csv(pa0_0, file.path(resDir, "Intersections", "Intersections_bioregions_birds_0_0.csv"),
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
write.csv(pa2.5_15, file.path(resDir, "Intersections", "Intersections_bioregions_birds_2.5_15.csv"),
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
write.csv(pa5.0_25, file.path(resDir, "Intersections", "Intersections_bioregions_birds_5.0_25.csv"),
          quote = FALSE, row.names = FALSE)

### 7.5% - 35%
pa7.5_35 <- convertToPA(rangeIntersects = intersections, # data frame of occupied area for species range maps x bioregions
                        threshAsia      = NULL,          # prop of species rangemap that needs to occur within trop. asia to be considered
                        threshBioregion = 0.075,         # prop of bioregion area that needs to be occupied to be considered present in bioregion
                        threshPropRange = 0.35,          # if < threshBioregion, prop of range to occur within bioregion to still be considered present
                        threshAsiaEnd   = 0.99,          # prop. total species range that occurs in tropical asia to be considered regional endemic (Australia is removed)
                        rangeAreaCol    = "Range_area",  # col name that contains species total range areas
                        regionNameCol   = "Bioregion",   # column name where unique region names are kept
                        regions         = regions,       # shapefile with bioregions
                        verbose         = TRUE)
write.csv(pa7.5_35, file.path(resDir, "Intersections", "Intersections_bioregions_birds_7.5_35.csv"),
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
write.csv(allSp, file.path(resDir, "Intersections", "Intersections_bioregions_birds_spList.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

### removing rangemaps doesn't clear up the mem used.
### We need to restart the session (at least in RStudio and Linux)
rm(rangemaps)
.rs.restartR()
