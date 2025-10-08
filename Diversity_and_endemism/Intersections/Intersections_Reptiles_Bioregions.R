####################################################################################################
### Run intersections between reptile rangemaps (GARD v1.7) and bioregions
### Charlie Marsh
### charliem2003@github
### 04/2024
###
### Regions map used: Bioregions_rangemaps
### Taxonomy used:    GARD v1.7 (follows www.reptile-database.org/)
###                   Uses reptile_checklist_2021_05, reptile_checklist_2015_03 and reptile_synonyms_2023_04
### Rangemaps used:   GARD v1.7
###
### CITATION:
### Roll, Uri; Meiri, Shai (2022). GARD 1.7 - updated global distributions for all
###  terrestrial reptiles [Dataset]. Dryad. https://doi.org/10.5061/dryad.9cnp5hqmb
### Roll U., et al. 2017. The global distribution of tetrapods reveals a need for targeted reptile
###  conservation. Nature Ecol. Evol. 1:1677–1682. (GARD 1.7). https://www.nature.com/articles/s41559-017-0332-2
### Caetano G.H.O., et al. 2022. Automated assessment reveals that the extinction risk of reptiles
###  is widely underestimated across space and phylogeny. PLoS Biology
### Uetz, P., Freed, P. & Hošek, J. (eds.) (2015) The Reptile Database, http://www.reptile-database.org, accessed April 2015
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
library(readr)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity_and_endemism"                         # project dir
rangeDir  <- file.path("IUCN", "rangemaps", "directory")      # dir that contains the rangemaps data
taxonDir  <- file.path("IUCN", "assessments", "directory")    # dir that contains the taxonomy data (for description year)
gadmDir   <- file.path("GADM", "directory")                   # dir that contains GADM in equal-area projection (called 'GADM_410_land_Equal_Area.gpkg')

### You shouldn't need to adjust these folders
regionDir <- file.path(projDir, "Data")                       # dir that contains the subregions data
funDir    <- file.path(projDir, "Analysis_functions")         # dir that contains the function scripts
resDir    <- file.path(projDir, "Intersections")              # dir to save results to
spDir     <- file.path(resDir,  "Intersections", "Reptiles")  # dir to save species intersections to

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
rangemaps <- st_read(file.path(rangeDir, "Gard_1_7_ranges.shp"))

# reproject to equal-area projection
rangemaps <- rangemaps %>%
  st_transform(st_crs(regions))

### total species list extracted from range maps - 10,914 species
spList <- unique(rangemaps$binomial)

#==================================================================================================#
#--------------------------------------- Run intersections ----------------------------------------#
#==================================================================================================#

### NOTE - parallelisation only works on linux or mac. For windows set threads = 1
intersections_rangemaps_parallel(rangemaps        = rangemaps,    # species range maps
                                 spNameCol        = "binomial",   # column name where species names are kept
                                 regions          = regions,      # shapefile with regions to summarise by
                                 regionNameCol    = "Bioregion",  # column name where region names are kept
                                 gadm             = gadm,         # global gadm for initial masking out where coastlines don't overlap
                                 crop_to_region   = TRUE,         # exclude species outside region
                                 areaThresh       = 10000000,     # area threshold in km2 for when to split the shapefile in to multiple parts
                                 spList           = spList,       # species list to subset rangemaps with
                                 spDir            = spDir,        # directory to save individual species intersections
                                 overwriteSpFiles = FALSE,        # overwrite existing species intersections, otherwise skips
                                 threads          = 7)            # no cores for parallelisation

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
spList <- gsub("_", " ", intersections$Species)

### clean taxonomy data frame to get taxonomic info and year described.
### Range maps have some species missing from the reptile database so we will glean information from
### the database at the 1st GARD version (2015) and 2nd GARD version (2021) and the latest version (2022)
tax2021 <- read.csv(file.path(taxonDir, "reptile_checklist_2021_05.csv")) %>%
  select(Species, Author)
tax2015 <- read.csv(file.path(taxonDir, "reptile_checklist_2015_03.csv")) %>%
  select(Species, Author) %>%
  filter(!Species %in% tax2021$Species)
taxonomy <- bind_rows(tax2015, tax2021) %>%
  mutate(Year = clean_publication_dates(Author))

### 22 missing species - take from the reptile database synonym list
synonyms <- read.csv(file.path(taxonDir, "reptile_synonyms_2023_04.csv"))
missing <- data.frame(Species = spList[!spList %in% taxonomy$Species]) %>%
  left_join(synonyms, join_by(Species == synonym)) %>%
  left_join(taxonomy, join_by("current.name" == "Species"))
taxonomy <- bind_rows(taxonomy, missing) %>%
  select(Species, Year) %>%
  mutate(Species = gsub(" ", "_", Species))

### merge intersections with taxonomy data frame
intersections <- intersections %>%
  left_join(taxonomy, by = "Species")

### add in genus and family (from rangemaps)
rangemaps <- st_read(file.path(rangeDir, "Gard_1_7_ranges.shp")) %>%
  st_drop_geometry() %>%
  select(binomial, family) %>%
  rename(Species = binomial, Family = family) %>%
  mutate(Species = gsub(" ", "_", Species))

intersections <- intersections %>%
  left_join(rangemaps, by = "Species") %>%
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
write.csv(pa0_0, file.path(resDir, "Intersections", "Intersections_bioregions_reptiles_0_0.csv"),
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
write.csv(pa2.5_15, file.path(resDir, "Intersections", "Intersections_bioregions_reptiles_2.5_15.csv"),
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
write.csv(pa5.0_25, file.path(resDir, "Intersections", "Intersections_bioregions_reptiles_5.0_25.csv"),
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
write.csv(pa7.5_35, file.path(resDir, "Intersections", "Intersections_bioregions_reptiles_7.5_35.csv"),
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
write.csv(allSp, file.path(resDir, "Intersections", "Intersections_bioregions_reptiles_spList.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#--------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

### removing rangemaps doesn't clear up the mem used.
### We need to restart the session (at least in RStudio and Linux)
rm(rangemaps)
.rs.restartR()
 