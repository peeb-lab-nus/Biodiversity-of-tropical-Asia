####################################################################################################
### Run intersections between marine bony fish (class: Actinopterygii) and marine ecoregions
### Note: this isn't comprehensive but also the definition has changed from the 2024-2 versions to
### now be "Species from the classes Actinopterygii, Chondrichthyes, Myxini, Petromyzonti and 
### Sarcopterygii which are found in marine habitats". Chondrichthyes will be removed as they
### are covered in sharks and rays.
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
spDir     <- file.path(resDir, "Intersections", "Bony_fish")    # dir to save species intersections to

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

### range maps - 15,383 polygons, 14,281 species names
rangemaps <- rbind(st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART1.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART2.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART3.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART4.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART5.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART6.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART7.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART8.shp")),
                   st_read(file.path(rangeDir, "polygons", "MARINEFISH_PART9.shp")))

### Read in point data. The IUCN is messed up and has text strings not in quotes that contain commas
### This makes it very tricky to read in. We'll have to go through line by line and select the first
### six columns and the final ten valid columns
pts <- readr::read_csv(file.path(rangeDir, "points", "MARINEFISH_points.csv"),
                       col_types = c("d", "d", "c", "d", "d", "d", "c", "d", rep("c", 18), "d", "d"))
ptsClean <- pts[1, c(1:6, 19:28)]
ptsClean$dec_lat  <- as.numeric(ptsClean$dec_lat)
ptsClean$dec_long <- as.numeric(ptsClean$dec_long)
pb <- txtProgressBar(2, nrow(pts), style = 3)
for(i in 2:nrow(pts)) {
  setTxtProgressBar(pb, i)
  ptsRow <- pts[i, ]
  validCells <- which(!is.na(ptsRow))
  # ptsRow <- ptsRow[c(1:6, validCells[(length(validCells) - 9):length(validCells)])]
  ptsRow <- ptsRow[c(1:6, (max(validCells) - 9):max(validCells))]
  names(ptsRow) <- names(ptsClean)
  ptsRow$dec_lat  <- as.numeric(ptsRow$dec_lat)
  ptsRow$dec_long <- as.numeric(ptsRow$dec_long)
  ptsClean <- rbind(ptsClean, ptsRow)
}
filter(ptsClean, is.na(dec_long) | is.na(dec_lat))
ptsClean <- filter(ptsClean, !is.na(dec_long) & !is.na(dec_lat))

table(ptsClean$presence, useNA = "always")
table(ptsClean$origin,   useNA = "always")
table(ptsClean$seasonal, useNA = "always")
table(ptsClean$legend,   useNA = "always")
table(ptsClean$class,    useNA = "always")
table(ptsClean$order_,   useNA = "always")
table(ptsClean$category, useNA = "always")
filter(ptsClean, is.na(legend))

### Species with point data - buffer by 25km (using equal-area projection) - 89,755 polygons, 158 species
ptsClean <- ptsClean %>%
  st_as_sf(coords = c("dec_long", "dec_lat")) %>%
  st_set_crs(st_crs(rangemaps)) %>%
  st_transform(st_crs(regions)) %>%
  st_buffer(dist = 25000) %>%
  st_transform(st_crs(rangemaps))

### Merge all maps together - 105,144 polygons, 14,283 species
rangemaps <- bind_rows(rangemaps, ptsClean)

### Separate out the bony fish (Actinopterygii) - 103,886 polygons, 13,103 species
rangemaps <- filter(rangemaps, class != "CHONDRICHTHYES")

### filter out non-relevant range polygons - 83,911 polygons, 13,084 species
rangemaps <- rangemaps %>%
  filter(!presence %in% c(4, 5, 6)) %>% # 1 = Extant; 2 = Prob. extant; 3 = Poss. extant; 4 = Poss. extinct; 5 = Extinct; 6 = Presence uncertain
  filter(!origin   %in% c(3, 4, 6)) %>% # 1 = Native; 2 = Reintroduced; 3 = Introduced; 4 = Vagrant; 5 = Origin uncertain; 6 = Assisted colonisation
  filter(!seasonal %in% c(3, 4))        # 1 = Resident; 2 = Breeding season; 3 = Non-breeding season; 4 = Passage; 5 = Seasonal occurrence uncertain

### total species list extracted from range maps - 13,084 species
spList <- unique(rangemaps$sci_name)

done <- c(list.files(spDir),
          list.files(file.path(spDir, "Outside_extent")),
          list.files(file.path(spDir, "Non-native")))
done <- gsub("_", " ", done)
spList <- spList[!spList %in% gsub(".csv", "", done)]
spList <- sample(spList, length(spList), replace = FALSE)
length(spList)

rangemaps <- filter(rangemaps, sci_name %in% spList)
gc()

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
                                 areaThresh       = 100000000,   # area threshold in km2 for when to split the shapefile in to multiple parts
                                 spList           = spList,     # species list to subset rangemaps with
                                 spDir            = spDir,      # directory to save individual species intersections
                                 overwriteSpFiles = FALSE,      # overwrite existing species intersections, otherwise skips
                                 threads          = 16)          # no cores for parallelisation

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
write.csv(pa5.0_25, file.path(resDir, "Intersections_bioregions_bony_fish_5.0_25.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

### removing rangemaps doesn't clear up the mem used.
### We need to restart the session (at least in RStudio and Linux)
rm(rangemaps)
.rs.restartR()
