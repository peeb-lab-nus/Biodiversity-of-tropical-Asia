####################################################################################################
### Intersections between Sponges and marine ecooregions
### Checklist generated from cleaned OBIS data by Buntaro Kusumoto
### Charlie Marsh
### charliem2003@github
### 12/2024
###
### Regions map used: MEOW_BTAS
### Taxonomy used:    WoRMS (downloaded 2024-11-11)
### Checklist used:   OBIS (occurrence.parquet - official data dump 23/07/2024)
###
### CITATION:
###   OBIS (2024) Ocean Biodiversity Information System. Intergovernmental Oceanographic Commission
###     of UNESCO. https://obis.org.
###   WoRMS Editorial Board (2024). World Register of Marine Species. Available from
###     https://www.marinespecies.org at VLIZ. Accessed 2024-11-11. doi:10.14284/170 
### 
### saves csv with the presence-absence of each species within each ecooregion
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(sf)

### Locations of data, scripts and results
baseDir   <- "/mnt/Work"
regionDir <- file.path(baseDir, "NUS", "BTAS_data", "Marine_ecoregions") # dir that contains the regions data
dataDir   <- file.path(baseDir, "spatial_data", "biodiversity", "checklists", "sponges") # dir that contains the checklist data
resDir    <- file.path(baseDir, "NUS", "BTAS", "Intersections") # dir to save results to
funDir    <- file.path(baseDir, "NUS", "BTAS", "Analysis_functions")   # dir that contains the function scripts
gadmDir   <- file.path(baseDir, "spatial_data", "regions", "countries", "GADM", "GADM_4.1") # dir with gadm data
taxonDir  <- file.path(baseDir, "spatial_data", "biodiversity", "occurrences", "marine", "WoRMS") # dir that contains worms taxonomic info

### Function for cleaning dates and extracting region info from inside brackets
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Read in obis intersections
all <- read.csv(file.path(dataDir, "obis_sp_province_matrix_sponge.csv")) %>%
  as_tibble()

### The OBIS point data doesn't contain the description year, but we can pull this out of worms
taxon <- read_delim(file.path(taxonDir, "taxon.txt"), delim = "\t") %>%
  select(acceptedNameUsageID, namePublishedInYear) %>%
  right_join(all, by = "acceptedNameUsageID", relationship = "many-to-many") %>%
  mutate(Species = gsub(" ", "_", ss)) %>%
  reframe(Year = min(namePublishedInYear, na.rm = TRUE),
          .by = "Species")

### Take earliest year
taxon$Year <- clean_publication_dates(taxon$Year, max = FALSE)

### Regions data
regions <- st_read(file.path(regionDir, "MEOW_BTAS.gpkg"))

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

### The intersections have been calculated, we just need to do some tidying up
all <- all %>%
  mutate(Species = gsub(" ", "_", ss)) %>%
  rename(Genus   = genus) %>%
  rename(Family  = family) %>%
  rename_with(., ~ gsub(".", "_", .x, fixed = TRUE)) %>%
  select(-ss, -acceptedNameUsageID, -specificEpithet, -order, -class, -phylum, -kingdom) %>%
  relocate(Species, Genus, Family, AsiaEndemic) %>%
  distinct()

### 8,731 species globally
length(unique(all$Species))

### Merge with year info
all <- all %>%
  left_join(taxon, by = "Species") %>%
  relocate(Year, .after = Family)

### There are nine species with no year. Get from worms online (www.marinespecies.org)
all$Year[all$Species == "Clathrina_ceylonensis"]                   <- 1905
all$Year[all$Species == "Callyspongia_(Cladochalina)_thurstoni"]   <- 1930

### Filter out species not native to tropical Asia - 1,442 species
ecoregions <- regions$PROVINCE
all <- all[rowSums(select(all, any_of(matches(ecoregions)))) > 0, ]

### We're going to exclude "Northwest_Australian_Shelf" and "Northeast_Australian_Shelf"
### Recalculate endemism so that Australian provinces count as non-endemic
ausRegions <- c("Northwest_Australian_Shelf", "Northeast_Australian_Shelf")
nonAusregions <- ecoregions[!ecoregions %in% ausRegions]

all <- all %>%
  mutate(AsiaEndemic = case_when(Northwest_Australian_Shelf == 1 | Northeast_Australian_Shelf == 1 ~ 0,
                                 .default = AsiaEndemic))

### total number of species native to region and number of species endemic to tropical asia (excl. Aus)
all %>%
  filter(rowSums(.[, nonAusregions]) > 0) %>%
  group_by(AsiaEndemic) %>%
  reframe(n = n())

### save results
write.csv(all, file.path(resDir, "Intersections", "Intersections_bioregions_sponges.csv"),
          quote = FALSE, row.names = FALSE)

####################################################################################################
### plot richness and turnover

all <- read.csv(file.path(resDir, "Intersections", "Intersections_bioregions_sponges.csv"))
regions <- st_read(file.path(regionDir, "MEOW_BTAS.gpkg")) %>%
  st_simplify(dTolerance = 1000)

rich <- data.frame(PROVINCE = names(all[, -c(1:5)]),
                   Rich = colSums(all[, -c(1:5)]))
rich <- left_join(regions, rich, by = "PROVINCE") %>%
  filter(!PROVINCE %in% c("Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))

ggplot() + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA),
        # legend.position = c(0.1, 0.02),
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        legend.justification = c("left", "bottom"),
        legend.key = element_blank()) +
  scale_fill_viridis_c(option = "C", name = NULL, trans = "log10") +
  guides(size = guide_legend(nrow = 2)) +
  scale_radius(range = c(3, 8),
               breaks = round(10 ^ seq(log10(min(rich$Rich)), log10(max(rich$Rich)), length.out = 5), 0),
               name = "Species Richness",
               trans = "log10") +
  geom_sf(data = rich, aes(fill = Rich)) + 
  geom_sf(data = st_centroid(rich),
          col = "black", pch = 21, aes(size = Rich, fill = Rich))
