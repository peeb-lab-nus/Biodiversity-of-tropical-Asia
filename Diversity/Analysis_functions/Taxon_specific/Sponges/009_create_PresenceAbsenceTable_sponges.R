#######################################

# Create presence/absence table for sponges
# Buntaro Kusumoto

#######################################
# Packages-----------------------------
library(readr) # to handle big table
library(tidyverse) # to code efficiently
library(data.table)
library(sf)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir <- "Diversity"                                                       # project dir

### You shouldn't need to adjust these folders
resDir  <- file.path(projDir, "Intersections", "Intersections", "Sponges") # dir that contains obis data for sponges
dataDir <- file.path(projDir, "Data")                                      # dir that contains MEOW ecoregions

# Biogeographic unit (MEOW)----------------
# https://doi.org/10.1641/B570707
# MEOW at provincial level overlapping BTAS zone
meow	<- st_read(file.path(dataDir, "MEOW_BTAS.gpkg"))

# Species occurrence points (OBIS) --------------------
# Only target sponges (Porifera)
file_name	<- file.path(resDir, "obis_worms_sponges.csv")
pt	<- read_csv(file_name)

# Check the contents--------------------
pt %>% colnames()

# Number of species
sp_all <- pt %>% pull(species) %>% unique() 
sp_all %>% length()
### 8741 in original

# Lon/lat coordinates
pt %>% pull(decimalLongitude) %>% summary()
###    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
###-180.000  -61.639   -4.897   -1.829    6.475  180.000 

pt %>% pull(decimalLatitude) %>% summary()
###   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### -78.87   10.65   49.97   27.56   52.83   88.78

# coordinate uncertainty
## coordinateUncertaintyInMeters: The horizontal distance (in meters) from the given dwc:decimalLatitude and dwc:decimalLongitude describing the smallest circle containing the whole of the dcterms:Location. Leave the value empty if the uncertainty is unknown, cannot be estimated, or is not applicable (because there are no coordinates). Zero is not a valid value for this term.
## coordinatePrecision: A decimal representation of the precision of the coordinates given in the dwc:decimalLatitude and dwc:decimalLongitude.
pt %>% pull(coordinateUncertaintyInMeters) %>% unique() %>% hist()
pt %>% pull(coordinateUncertaintyInMeters) %>% unique() %>% quantile(0.95, na.rm=T)
pt %>% pull(coordinatePrecision) %>% unique() 
pt %>% pull(coordinatePrecision) %>% hist() 

# Filter the points based on coordinate uncertainty (horiz dist) 
# Removing points whose uncertainty is lager than 100 km
pt2 <- pt %>% filter( coordinateUncertaintyInMeters < 100000 | is.na(coordinateUncertaintyInMeters) )

# Habitat-------------------------------
pt %>% pull(marine) %>% unique()
pt %>% pull(brackish) %>% unique()
pt %>% pull(freshwater) %>% unique()
pt %>% pull(terrestrial) %>% unique()
# Keep only points with marine==TRUE (marine species)
pt3	<- pt2 %>% filter(marine==TRUE)

# Occurrence state-----------------------
# not useful for filtering
pt3 %>% pull(occurrenceStatus) %>% unique()

#absence---------------------------------
# no absence data
pt3 %>% pull(absence) %>% unique()

# Number of occurrence points -----------
## original: 320,991
## filtered: 319,878 (-1113)

# Number of species ---------------------
sp_marine	<- pt3 %>% pull(species) %>% unique() 
sp_marine %>% length()
## original: 8741
## filtered: 8731 (-10)
## removed species: 
sp_all[!sp_all %in% sp_marine ]
### "Protosuberites collaris"
### "Amorphinopsis mollis"         
### "Racodiscula spiroglyphi"
### "Metschnikowia tuberculata"    
### "Ircinia solida"
### "Dactylia elegans"             
### "Spongilla alba"
### "Acheliderma fistulatum"       
### "Protosuberites lacustris"
### "Protosuberites aquaedulcioris"

# Intersect occurrence points with MEOW----
# Convert pt3 to spatial data
pt_spat	<- st_as_sf(pt3, coords = c("decimalLongitude", "decimalLatitude"))
st_crs(pt_spat)	<- st_crs(4326)
pt_spat <- st_transform(pt_spat, st_crs(meow))

# Intersect ------------------------
meow_selected	<- meow %>% select("PROV_CODE", "PROVINCE")

pt_meow	<- st_join(pt_spat,
                   meow_selected,
                   join = st_intersects)
pt_meow_tib	<- pt_meow %>% as_tibble()

# Select necessary columns and remove duplicates
# Number of species: 132,270 records at province level
# col_select	<- c(colnames(pt3)[2:10], "PROV_CODE")
pt_selected	<- pt_meow_tib %>% distinct()

# Transform long-data to wide-----
# Number of species: 8731 sp
spXsite	<- pt_selected %>% 
  select(species, specificEpithet, genus, family, order, descriptionYear, PROVINCE) %>%
  mutate(presence = 1) %>%
  distinct() %>%
  pivot_wider(names_from = PROVINCE, values_from = presence, values_fill = 0) %>%
  rename(OutsideAsia = `NA`)

# Endemic species------------------
# OutsideAsia==0 & sum within asia > 0
isasia	<- spXsite %>% select(Andaman:`Java Transitional`) %>% rowSums()
spXsite	<- spXsite %>% mutate(AsiaEndemic = (isasia > 0 & spXsite$OutsideAsia==0) %>% as.numeric())
spXsite	<- spXsite %>% relocate(AsiaEndemic, .after = OutsideAsia)

# Number of asia endemic
spXsite$AsiaEndemic %>% sum()
spXsite$AsiaEndemic %>% sum()/ (sp_marine %>% length())
# 833 sp. (9.5%)

# Output----------------------------
write_csv(spXsite, file.path(resDir, "obis_sp_province_matrix_sponge.csv"), na="")

# Tidy up
rm(list = ls())

