####################################################################################################
### Run intersections between freshwater fish (IUCN) and bioregions
### Charlie Marsh
### charliem2003@github
### 08/2024
###
### Regions map used: Bioregions_rangemaps
### Taxonomy used:    IUCN 2024.1 - downloaded 23/08/2024
### Rangemaps used:   IUCN 2024.1 - downloaded 01/08/2024
###
### CITATION: IUCN 2024. The IUCN Red List of Threatened Species. 2024.1.
###      https://www.iucnredlist.org. Downloaded on 01/08/2024.
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
projDir   <- "Diversity_and_endemism"                               # project dir
rangeDir  <- file.path("IUCN", "rangemaps", "directory")            # dir that contains the rangemaps data
hydroDir  <- file.path("IUCN", "HydroBASINS", "directory")          # dir containing hydroBASINs data
taxonDir  <- file.path("IUCN", "assessments", "directory")          # dir that contains the taxonomy data (for description year)
gadmDir   <- file.path("GADM", "directory")                         # dir that contains GADM in equal-area projection (called 'GADM_410_land_Equal_Area.gpkg')

### You shouldn't need to adjust these folders
regionDir <- file.path(projDir, "Data")                             # dir that contains the subregions data
funDir    <- file.path(projDir, "Analysis_functions")               # dir that contains the function scripts
resDir    <- file.path(projDir, "Intersections")                    # dir to save results to
spDir     <- file.path(resDir,  "Intersections", "Freshwater_fish") # dir to save species intersections to

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

### hydroBASINs data (IUCN versions)
hydroBasin <- bind_rows(st_read(file.path(hydroDir, "HydroBASINS_level08", "HydroBASINS_level08_w_attributes_2023_08.shp")),
                        st_read(file.path(hydroDir, "HydroBASINS_level10", "HydroBASINS_level10_w_attributes_2023_08.shp")),
                        st_read(file.path(hydroDir, "HydroBASINS_level12", "HydroBASINS_level12_w_attributes_2023_08.shp"))) %>%
  select(hybas_id)

### range maps - 21,234 polygons, 14,569 species names
rangemaps <- rbind(st_read(file.path(rangeDir, "polygons", "FW_FISH_PART1.shp")),
                   st_read(file.path(rangeDir, "polygons", "FW_FISH_PART2.shp")))

### For some species different subpopulations have their own species name e.g. "Salmo salar Nunavik subpopulation"
### Rename to just the species name - 21,234 polygons, 14,484 species
subpops <- unique(rangemaps$sci_name[grepl("subpopulation", rangemaps$sci_name)])
for(i in subpops) {
  rangemaps$sci_name[rangemaps$sci_name == i] <- paste(strsplit(i, " ")[[1]][1:2], collapse = " ")
}

### Species with hydrobasin level data - 9,461,951 polygons, 14,011 species
basins <- rbind(read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part1.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part2.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part3.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part4.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part5.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part6.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part7.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part8.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part9.csv")),
                read.csv(file.path(rangeDir, "hydrobasins", "fish_hybas_table_part10.csv")))

### Again, rename subpopulations to just the species name - 21,234 polygons, 13,926 species
subpops <- unique(basins$sci_name[grepl("subpopulation", basins$sci_name)])
for(i in subpops) {
  basins$sci_name[basins$sci_name == i] <- paste(strsplit(i, " ")[[1]][1:2], collapse = " ")
}

### 36 species don't seem to intersect with hybroBASINs but they are all featured in the range maps data
missing <- unique(basins$sci_name[basins$hybas_id %in% basins$hybas_id[!basins$hybas_id %in% hydroBasin$hybas_id]])
missing[!missing %in% rangemaps$sci_name]

### Merge with hydroBASINS shapefile and remove species already in range maps - 48,777 polygons, 132 species
basins <- left_join(hydroBasin, basins, by = "hybas_id", relationship = "many-to-many") %>%
  filter(!is.na(sci_name)) %>%
  filter(sci_name %in% basins$sci_name[!basins$sci_name %in% rangemaps$sci_name]) %>%
  select(-hybas_id, -objectid)

### Species with point data - buffer by 25km (using equal-area projection) - 255,128 polygons, 3,235 species
pts <- read.csv(file.path(rangeDir, "points", "FW_FISH_points.csv")) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(st_crs(rangemaps)) %>%
  st_transform(st_crs(regions)) %>%
  st_buffer(dist = 25000) %>%
  st_transform(st_crs(rangemaps))

### Merge all maps together - 325,139 polygons, 14,621 species
rangemaps <- bind_rows(rangemaps, basins, pts)

# ### Remove marine species - 323,198 polygons, 13,721 species
rangemaps <- filter(rangemaps, marine == "false" | is.na(marine))

### filter out non-relevant range polygons - 317,703 polygons, 13,519 species
rangemaps <- rangemaps %>%
  filter(!presence %in% c(4, 5, 6)) %>% # 1 = Extant; 2 = Prob. extant; 3 = Poss. extant; 4 = Poss. extinct; 5 = Extinct; 6 = Presence uncertain
  filter(!origin   %in% c(3, 4, 6)) %>% # 1 = Native; 2 = Reintroduced; 3 = Introduced; 4 = Vagrant; 5 = Origin uncertain; 6 = Assisted colonisation
  filter(!seasonal %in% c(3, 4))        # 1 = Resident; 2 = Breeding season; 3 = Non-breeding season; 4 = Passage; 5 = Seasonal occurrence uncertain

### total species list extracted from range maps - 13,519 species
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

### NOTE - parallelisation only works on linux or mac. For windows set threads = 1
intersections_rangemaps_parallel(rangemaps        = rangemaps,    # species rangemaps
                                 spNameCol        = "sci_name",   # column name where species names are kept
                                 regions          = regions,      # shapefile with regions to summarise by
                                 regionNameCol    = "Bioregion",  # column name where region names are kept
                                 gadm             = gadm,         # global gadm for initial masking out where coastlines don't overlap
                                 marine           = FALSE,        # if TRUE masks out land areas. If FALSE masks out marine areas
                                 crop_to_region   = TRUE,         # exclude species outside region
                                 areaThresh       = 10000000,     # area threshold in km2 for when to split the shapefile in to multiple parts
                                 spList           = spList,       # species list to subset rangemaps with
                                 spDir            = spDir,        # directory to save individual species intersections
                                 overwriteSpFiles = FALSE,        # overwrite existing species intersections, otherwise skips
                                 threads          = 4)            # no cores for parallelisation

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

### There are five species with no info in the IUCN taxonomy file. Add in manually
intersections$Genus[intersections$Species  == "Arothron_manilensis"] <- "Arothron"         # https://www.iucnredlist.org/species/193607/2246630
intersections$Family[intersections$Species == "Arothron_manilensis"] <- "Tetraodontidae"
intersections$Year[intersections$Species   == "Arothron_manilensis"] <- 1822
intersections$Genus[intersections$Species  == "Callogobius_hasseltii"] <- "Callogobius"    # https://www.iucnredlist.org/species/193034/2186850
intersections$Family[intersections$Species == "Callogobius_hasseltii"] <- "Gobiidae"
intersections$Year[intersections$Species   == "Callogobius_hasseltii"] <- 1851
intersections$Genus[intersections$Species  == "Mugilogobius_abei"] <- "Mugilogobius"       # https://www.iucnredlist.org/species/181269/1715700
intersections$Family[intersections$Species == "Mugilogobius_abei"] <- "Gobiidae"
intersections$Year[intersections$Species   == "Mugilogobius_abei"] <- 1901
intersections$Genus[intersections$Species  == "Pastinachus_sephen"] <- "Pastinachus"       # https://www.iucnredlist.org/species/70682503/109922153
intersections$Family[intersections$Species == "Pastinachus_sephen"] <- "Dasyatidae"
intersections$Year[intersections$Species   == "Pastinachus_sephen"] <- 1775
intersections$Genus[intersections$Species  == "Urogymnus_acanthobothrium"] <- "Urogymnus"  # https://www.iucnredlist.org/species/200322937/200322952 
intersections$Family[intersections$Species == "Urogymnus_acanthobothrium"] <- "Dasyatidae"
intersections$Year[intersections$Species   == "Urogymnus_acanthobothrium"] <- 2016

### One species is currently being described -> put as 2024
intersections$Year[intersections$Species == "Melanotaenia_sp._nov._'Malanda'"] <- 2024

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
write.csv(pa0_0, file.path(resDir, "Intersections", "Intersections_bioregions_freshwater_fish_0_0.csv"),
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
write.csv(pa2.5_15, file.path(resDir, "Intersections", "Intersections_bioregions_freshwater_fish_2.5_15.csv"),
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
write.csv(pa5.0_25, file.path(resDir, "Intersections", "Intersections_bioregions_freshwater_fish_5.0_25.csv"),
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
write.csv(pa7.5_35, file.path(resDir, "Intersections", "Intersections_bioregions_freshwater_fish_7.5_35.csv"),
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
write.csv(allSp, file.path(resDir, "Intersections", "Intersections_bioregions_freshwater_fish_spList.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

### removing rangemaps doesn't clear up the mem used.
### We need to restart the session (at least in RStudio and Linux)
rm(rangemaps)
.rs.restartR()
