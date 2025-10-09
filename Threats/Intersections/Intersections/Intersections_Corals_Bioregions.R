####################################################################################################
### Run intersections between hard corals (IUCN) and marine ecoregions
### Charlie Marsh
### charliem2003@github
### 05/2025
###
### Regions map used: MEOW_BTAS
### Taxonomy used:    IUCN 2025-1 - downloaded 09/05/2025
### Rangemaps used:   IUCN 2025-1 - downloaded 09/05/2025
###
### CITATION: IUCN 2025. The IUCN Red List of Threatened Species. 2025-1.
###   https://www.iucnredlist.org. Downloaded on 09/05/2025.
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
library(parallel)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Threats"                                          # project dir
rangeDir  <- file.path("IUCN", "rangemaps", "directory")        # dir that contains the rangemaps data
taxonDir  <- file.path("IUCN", "assessments", "directory")      # dir that contains the iucn assessments
gadmDir   <- file.path("GADM", "directory")                     # dir that contains GADM in equal-area projection (called 'GADM_410_land_Equal_Area.gpkg')

### You shouldn't need to adjust these folders
regionDir <- file.path(projDir, "Data")                         # dir that contains the subregions data
funDir    <- file.path(projDir, "Analysis_functions")           # dir that contains the function scripts
resDir    <- file.path(projDir, "Intersections")                # dir to save results to
spDir     <- file.path(resDir, "Intersections", "Corals")       # dir to save species intersections to

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
regions <- st_read(file.path(regionDir, "MEOW_BTAS.gpkg"))

### global GADM for masking out dodgy coastlines in IUCN range maps on global range maps
gadm <- st_read(file.path(gadmDir, "GADM_410_land_Equal_Area.gpkg"))

### range maps - 892 polygons, 892 species - 1 species with no range map (Pseudocyathoceras avis)
rangemaps <- rbind(st_read(file.path(rangeDir, "REEF_FORMING_CORALS_PART1.shp")),
                   st_read(file.path(rangeDir, "REEF_FORMING_CORALS_PART2.shp")),
                   st_read(file.path(rangeDir, "REEF_FORMING_CORALS_PART3.shp")))

## filter out non-relevant range polygons - 903 polygons, 903 species
## Not relevant for corals - all rangemps are 1, 1, 1
# presence: 1 = Extant; 2 = Prob. extant; 3 = Poss. extant; 4 = Poss. extinct; 5 = Extinct; 6 = Presence uncertain
# origin: 1 = Native; 2 = Reintroduced; 3 = Introduced; 4 = Vagrant; 5 = Origin uncertain; 6 = Assisted colonisation
# seasonal: 1 = Resident; 2 = Breeding season; 3 = Non-breeding season; 4 = Passage; 5 = Seasonal occurrence uncertain

### total species list extracted from range maps - 892 species
spList <- unique(rangemaps$sci_name)

done <- c(list.files(spDir),
          list.files(file.path(spDir, "Outside_extent")),
          list.files(file.path(spDir, "Non-native")))
done <- gsub("_", " ", done)
spList <- spList[!spList %in% gsub(".csv", "", done)]
length(spList)

rangemaps <- filter(rangemaps, sci_name %in% spList)

### reproject to equal-area projection
rangemaps <- rangemaps %>%
  st_transform(st_crs(regions))

#==================================================================================================#
#--------------------------------------- Run intersections ----------------------------------------#
#==================================================================================================#

intersections_rangemaps_parallel(rangemaps        = rangemaps,  # species range maps
                                 spNameCol        = "sci_name", # column name where species names are kept
                                 regions          = regions,    # shapefile with regions to summarise by
                                 regionNameCol    = "PROVINCE", # column name where region names are kept
                                 gadm             = gadm,       # global gadm for initial masking out where coastlines don't overlap
                                 marine           = TRUE,       # if TRUE masks out land areas. If FALSE masks out marine areas
                                 crop_to_region   = TRUE,       # exclude species outside region
                                 areaThresh       = 10000000,   # area threshold in km2 for when to split the shapefile in to multiple parts
                                 spList           = spList,     # species list to subset rangemaps with
                                 spDir            = spDir,      # directory to save individual species intersections
                                 overwriteSpFiles = FALSE,      # overwrite existing species intersections, otherwise skips
                                 threads          = 20)          # no cores for parallelisation

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

### Get taxonomic and year info
taxonomy <- read.csv(file.path(taxonDir, "taxonomy.csv")) %>%
  rename(Species = scientificName, Genus = genusName, Family = familyName) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  mutate(Year = clean_publication_dates(authority, max = TRUE)) %>%
  select(Species, Genus, Family, Year)

### merge intersections with taxonomy data frame
intersections <- intersections %>%
  left_join(taxonomy, by = "Species") %>%
  relocate(Species, Genus, Family, Year)

#==================================================================================================#
#---------------------------------- Convert to presence-absence -----------------------------------#
#==================================================================================================#

### 5% - 25%
pa5.0_25 <- convertToPA(rangeIntersects = intersections, # data frame of occupied area for species range maps x bioregions
                        threshAsia      = NULL,          # prop of species rangemap that needs to occur within trop. asia to be considered
                        threshBioregion = 0.050,         # prop of bioregion area that needs to be occupied to be considered present in bioregion
                        threshPropRange = 0.25,          # if < threshBioregion, prop of range to occur within bioregion to still be considered present
                        threshAsiaEnd   = 0.95,          # prop. total species range that occurs in tropical asia to be considered regional endemic (Australia is removed)
                        rangeAreaCol    = "Range_area",  # col name that contains species total range areas
                        regionNameCol   = "PROVINCE",    # column name where unique region names are kept
                        regions         = regions,       # shapefile with bioregions
                        verbose         = TRUE)
write.csv(pa5.0_25, file.path(resDir, "Intersections_bioregions_corals_5.0_25.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

### removing rangemaps doesn't clear up the mem used.
### We need to restart the session (at least in RStudio and Linux)
rm(rangemaps)
.rs.restartR()
