####################################################################################################
### Intersections between Trichoptera and bioregions
### Charlie Marsh
### Raw checklists processed by Leshon Lee and team
### charliem2003@github
### 06/2024
###
### Regions map used: Bioregions_checklists_noChinaJapan
### Taxonomy used:    Trichoptera World Checklist (TWC)
### Checklist used:   https://trichopt.app.clemson.edu
###
### CITATION:
###   https://trichopt.app.clemson.edu
###   Morse, J.C., 2011. The trichoptera world checklist. Zoosymposia 5, 372–380.
###      https://doi.org/10.11646/zoosymposia.5.1.29
###
### saves csv with the presence-absence of each species within each bioregion
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

# rm(list = ls())

### Libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
# library(ggplot2)
library(sf)

### Locations of data, scripts and results
baseDir   <- "/mnt/Work"
regionDir <- file.path(baseDir, "NUS", "BTAS_data", "Bioregions") # dir that contains the regions data
dataDir   <- file.path(baseDir, "spatial_data", "biodiversity", "checklists", "trichoptera") # dir that contains the checklist data
resDir    <- file.path(baseDir, "NUS", "BTAS", "Intersections") # dir to save results to
funDir    <- file.path(baseDir, "NUS", "BTAS", "Analysis_functions")   # dir that contains the function scripts
gadmDir   <- file.path(baseDir, "spatial_data", "regions", "countries", "GADM", "GADM_4.1") # dir with gadm data




### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity"                                  # project dir
dataDir   <- file.path("MilliBase", "data", "directory")  # dir that contains the checklist data
gadmDir   <- file.path("GADM", "data", "directory")       # dir that contains GADM data

### You shouldn't need to adjust these folders
resDir    <- file.path(projDir, "Intersections")          # dir to save results to
funDir    <- file.path(projDir, "Analysis_functions")     # dir that contains the function scripts
lookupDir <- file.path(resDir, "Look-up_tables")          # dir that contains look-up table for assigning bioregions

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Function for cleaning dates and extracting region info from inside brackets
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))
source(file.path(funDir, "Intersections", "extractRegionInfo.R"))
source(file.path(funDir, "Intersections", "extractParantheses.R"))

### Codes used to assign distribution info to bioregions
bioregions <- read.csv(file.path(lookupDir, "country_codes_BTAS_trichoptera.csv"))

### global list of all country and province names from GADM and bioregions
# countries <- tolower(bioregions$country[sapply(strsplit(bioregions$country, " "), length) == 1])
countries <- tolower(bioregions$country)
gadmCountries  <- st_read(file.path(gadmDir, "gadm_410-levels.gpkg"), layer = "ADM_0") %>%
  st_drop_geometry()
gadmProvinces <- st_read(file.path(gadmDir, "gadm_410-levels.gpkg"), layer = "ADM_1") %>%
  st_drop_geometry()
gadm <- as.vector(tolower(c(gadmCountries$COUNTRY, gadmProvinces$NAME_1, countries, "far east", "peninsular")))

### Cleaned data
dd <- read_csv(file.path(dataDir, "Oriental_Australasian_Trichoptera_compiled_100624.csv")) %>%
  select(`Species name`, Genus...10, Family, Year, `Distribution(with duplicates)`) %>% #, `Type locality`, Distribution) %>%
  mutate(Family  = sapply(Family, function(x) strsplit(x, " ")[[1]][2])) %>%
  dplyr::rename(Species = `Species name`, Genus = Genus...10, country = `Distribution(with duplicates)`) %>%
  mutate(Year = clean_publication_dates(Year, max = TRUE))

### Some manual corrections
dd$country[dd$Species == "Cheumatopsyche banksi"]          <- "Vietnam, Nepal, Indonesia (Sumatra), Thailand, China (Fuzhou)"
dd$country[dd$Species == "Diplectrona aurovittata"]        <- "Vietnam, Indonesia, Indonesia (Bali, Java, Lombok, Sulawesi, Sumatra), Thailand, Laos, Malaysia, Malaysia (Peninsular), China (Guangdong, Guangxi, Sichuan), Nepal, Myanmar"
dd$country[dd$Species == "Dipseudopsis robustior"]         <- "Vietnam, Thailand, Indonesia (Sumatra), Cambodia, Myanmar, Malaysia (Peninsular)"
dd$country[dd$Species == "Polyplectropus anakgugur"]       <- "Vietnam, China (Guangxi, Guangdong, Anhui, Henan, Guizhou, Jiangxi, Sichuan), Malaysia, Malaysia (Peninsular), Vietnam"
dd$country[dd$Species == "Potamyia phaidra"]               <- "Vietnam, Thailand, Laos, India, Indonesia (Java, Sumatra), Thailand"
dd$country[dd$Species == "Pseudoneureclipsis unguiculata"] <- "Indonesia (Sumatra), Philippines (Luzon, Mindanao, Basilan, Leyte, Mindoro)"

### and append to cleaned data
dd <- dd %>%
  mutate(Species = iconv(Species, from = "ISO-8859-1", to = "UTF-8")) %>%
  mutate(Species = gsub(" ", "_", Species))

#==================================================================================================#
#----------------- Get species missing from OL and AU subsets from global database ----------------#
#==================================================================================================#

### Global data - 16,891 unique species (note some species are on multiple rows)
raw <- read_csv(file.path(dataDir, "worldtrichoptera.csv")) %>%
  filter(status %in% c("extant species", "extant speciesextant speciesextant species")) %>%
  mutate(across(where(is.character), function(x) iconv(x, from = "ISO-8859-2", to = "UTF-8"))) %>%
  mutate(family = iconv(family, from = "UTF-8", to = "UTF-8", sub = "")) %>%
  mutate(taxon  = iconv(taxon,  from = "UTF-8", to = "UTF-8", sub = "")) %>%
  select("taxon", "status", "family", "genus", "distribution notes", "type country", "region(s)") %>%
  mutate(genus = gsub("^ *|(?<= ) | *$", "", genus, perl = TRUE)) %>% # remove double spaces
  mutate(taxon = gsub("^ *|(?<= ) | *$", "", taxon, perl = TRUE)) %>% # remove double spaces
  mutate(Family  = sapply(family, function(x) strsplit(x, " ")[[1]][2])) %>%
  mutate(Genus   = sapply(genus,  function(x) strsplit(x, " ")[[1]][2])) %>%
  mutate(Species = sapply(taxon,  function(x) strsplit(x, " ")[[1]][2])) %>%
  mutate(Species = gsub("[(+]", "", Species)) %>%
  mutate(Species = paste(Genus, Species, sep = "_"))

### Some manual corrections due to errors in the spreadsheet
raw$Genus[raw$Family == "Hydropsychidae" & raw$Species == "peruana"] <- "Smicridea"
raw$Genus[raw$genus == '"Limnephilus"'] <- "Limnephilus"
raw$Species[raw$Species == "NA_chereshnevi"]   <- "Limnephilus_chereshnevi"
raw$Species[raw$Species == "NA_fenestratus"]   <- "Limnephilus_fenestratus"
raw$Species[raw$Species == "NA_fumosus"]       <- "Limnephilus_fumosus"
raw$Species[raw$Species == "NA_peltus"]        <- "Limnephilus_peltus"
raw$Species[raw$Species == "NA_samoedus"]      <- "Limnephilus_samoedus"
raw$Species[raw$Species == "NA_submonilifer"]  <- "Limnephilus_submonilifer"
raw$Species[raw$Species == "NA_tricalcaratus"] <- "Limnephilus_tricalcaratus"
raw$Species[raw$Species == "NA_thor"]          <- "Montiphylax_thor"
raw$Genus[raw$Species == "Montiphylax_thor"]   <- "Montiphylax"

### Three species have messed up genus_species names, but also don't have any geog info. Remove
raw <- raw[!(raw$Family == "Hydropsychidae" & raw$Species == "NA_Plectropsyche"), ]
raw <- raw[!(raw$Family == "Hydropsychidae" & raw$Species == "NA_Streptopsyche"), ]
raw <- raw[!(raw$Family == "Limnephilidae"  & raw$Species == "NA_Nemotaulius"), ]

# ### Extract publication year. Some entries have no valid year in the field (there are others but non-Asian)
raw <- mutate(raw, Year = clean_publication_dates(taxon, max = FALSE, first = TRUE))
raw$Year[raw$Year < 1758 | raw$Year > 2024] <- NA

### Some species are missing year. Add manually
raw$Year[raw$Species == "Helicopsyche_petri"] <- 1986
raw$Year[raw$Species == "Hydropsyche_moselyi"] <- 1962
raw$Year[raw$Species == "Hydroptila_spirula"] <- 2002
raw$Year[raw$Species == "Macrostemum_pseudoneura"] <- 1865
raw$Year[raw$Species == "Macrostemum_quinquepunctatum"] <- 1920
raw$Year[raw$Species == "Oecetis_angustipennis"] <- 1936
raw$Year[raw$Species == "Oecetis_caelum"] <- 1990
raw$Year[raw$Species == "Rhyacophila_chayulpa"] <- 1970
raw$Year[raw$Species == "Triaenodes_aureus"] <- 1962
raw$Year[raw$Species == "Triaenodes_longispinus"] <- 1962

### and a few discrepancies between the oriental and global dataset. Keep the oriental years
raw$Year[raw$Species == "Agapetus_agtuuganonis"] <- dd$Year[dd$Species == "Agapetus_agtuuganonis"]
raw$Year[raw$Species == "Chimarra_bicolor"]      <- dd$Year[dd$Species == "Chimarra_bicolor"]
raw$Year[raw$Species == "Hydropsyche_closi"]     <- dd$Year[dd$Species == "Hydropsyche_closi"]
raw$Year[raw$Species == "Hydropsyche_kimminsi"]  <- dd$Year[dd$Species == "Hydropsyche_kimminsi"]
raw$Year[raw$Species == "Hydropsyche_tanung"]    <- dd$Year[dd$Species == "Hydropsyche_tanung"]

### There are some species with multiple year entries - select one based on https://eol.org/
raw$Year[raw$Species == "Cheumatopsyche_pallida"] <- 1920
raw$Year[raw$Species == "Chimarra_cornuta"]       <- 1959
raw$Year[raw$Species == "Oecetis_spinifera"]      <- 1998
raw$Year[raw$Species == "Allogamus_auricollis"]   <- 1834
raw$Year[raw$Species == "Rhyacophila_kumanskii"]  <- 1988

### Deal with duplicated species - append all information together
raw <- reframe(raw,
               Family               = paste(unique(Family), collapse = ", "),
               Genus                = paste(unique(Genus), collapse = ", "),
               Year                 = paste(unique(Year), collapse = ", "),
               `distribution notes` = paste(unique(`distribution notes`), collapse = ", "),
               `type country`       = paste(unique(`type country`), collapse = ", "),
               .by = Species)

### Extract distribution info - anything that matches names in gadm etc
raw <- raw %>%
  mutate(`distribution notes` = iconv(`distribution notes`, from = "ISO-8859-1", to = "UTF-8")) %>%
  mutate(verbatimDistribution = paste(`type country`, `distribution notes`, sep = ", ")) %>%
  mutate(country = sapply(verbatimDistribution,
                          function(x) extractRegionInfo(distribution = x, regionNames = gadm))) %>%
  filter(!is.na(country))

#==================================================================================================#
#--------------------------- Extract country info from missing species ----------------------------#
#==================================================================================================#

### get simplified set of regional names and find any matches in 'raw'
missing <- raw %>%
  mutate(Year = as.numeric(Year)) %>%
  select(Species, Genus, Family, Year, country)

### Append to original tropical asia data and separate to individual lines - 16551 species
dd <- bind_rows(dd, missing) %>%
  reframe(Species = unique(Species),
          Genus   = unique(Genus),
          Family  = unique(Family),
          Year    = unique(Year),
          country = paste(country, sep = ", "),
          .by = Species) %>%
  mutate(country = sapply(country, function(x) extractParantheses(x))) %>%
  separate_longer_delim(country, delim = ", ")

### Sometimes countries weren't separated by commas -> add commas between countries (but not provinces)
dd <- dd %>%
  mutate(country = sapply(country, function(x) {
    sep <- tolower(strsplit(x, " ")[[1]])
    countries <- which(sep %in% tolower(gadmCountries$COUNTRY))
    for(i in 1:length(sep)) {
      if((i %in% countries & (i + 1) %in% countries) | (!i %in% countries & (i + 1) %in% countries)) {
        sep[i] <- paste0(sep[i], ",")
      }
    }
    return(paste(sep, collapse = " "))
  } )) %>%
  separate_longer_delim(country, delim = ", ")

### Remove duplicated entries and trim white space
dd <- dd %>%
  mutate(country = trimws(country, "both")) %>%
  mutate(country = gsub("^ *|(?<= ) | *$", "", country, perl = TRUE)) %>% # remove double spaces
  distinct()

### Out of 1712 species with records in China and/or Japan there are 425 species in raw data with
### only 'China' or 'Japan' and no other province info
CN_JP_provinces <- c("anhui", "beijing", "chamdo", "chemo", "chokiang", "chongqing", "chusan",
                     "fujian", "fuzhou", "gansu", "guangxi", "guang-xi", "guangdong",
                     "guizhou", "hainan", "hebei", "heilongjiang", "henan", "hiroshima", "hokkaido",
                     "honanfu", "hong kong", "hubei", "hunan", "hyogo", "ishikawa", "jiangsu", "jiangxi",
                     "jiangxu", "jianxi", "jilin", "kanagawa", "kham", "korla", "kyoto", "liaoning",
                     "japan central", "manchuria", "mie", "nagano", "okayama", "okinawa", "omei shan",
                     "qinghai", "sai-chin-kette", "shaanxi", "shan", "shandong", "shanghai", "shantung",
                     "shanxi", "shenxi", "sichuan", "songpan", "tibet", "tokyo", "tsingtau",
                     "tienmuishan", "xinjiang", "xizang", "yunnan", "zhejiang")
length(unique(dd$Species[grepl(paste0("china|japan|", CN_JP_provinces, collapse = "|"), dd$country)]))
dd %>%
  filter(Species %in% dd$Species[grepl(paste0("china|japan|", CN_JP_provinces, collapse = "|"), dd$country)]) %>%
  filter(!Species %in% dd$Species[grepl(paste(CN_JP_provinces, collapse = "|"), dd$country)]) %>%
  summarise(length(unique(Species)))
#***# -> We're happy to keep yunnan, guangxi, okinawa etc in data

### If country == Indonesia or Malaysia, but more detailed province info also given, remove those rows
spMalaysia <- unique(dd$Species[dd$Species %in% dd$Species[dd$country == "malaysia"]])
for(i in 1:length(spMalaysia)) {
  sp <- filter(dd, Species == spMalaysia[i])
  if(sum(grepl(paste(c("malaysia ", "bboorrnneeoo", "borneo", "brunei", "sabah", "johor"),
                     collapse = "|"), sp$country)) >= 1) {
    dd <- filter(dd, Species != spMalaysia[i] | country != "malaysia")
  }
}

spIndonesia <- unique(dd$Species[dd$Species %in% dd$Species[dd$country == "indonesia"]])
for(i in 1:length(spMalaysia)) {
  sp <- filter(dd, Species == spIndonesia[i])
  if(sum(grepl(paste(c("indonesia ", "sulawesi", "misool island", "new britain", "papua", "rossel island",
                       "guinea", "kei island", "timor", "java", "bali", "kalimantan", "borneo", "lombok"),
                     collapse = "|"), sp$country)) >= 1) {
    dd <- filter(dd, Species != spIndonesia[i] | country != "indonesia")
  }
}

### assign to BTAS bioregions
final <- left_join(dd, bioregions, by = "country", relationship = "many-to-many") %>%
  mutate(AsiaNative = case_when(AsiaNative == 0   ~ 0,
                                AsiaNative == 1   ~ 1,
                                is.na(AsiaNative) ~ 0))

### Remove geog descriptions that are likely parsing errors that aren't def. outside trop. asia
final <- filter(final, !country %in% c("", "gulf", "lakes", "mus", "na", "new", "peninsular", "rivers"))

### which species are native to tropical asia
native <- final %>%
  group_by(Species) %>%
  reframe(native = max(AsiaNative, na.rm = TRUE)) %>%
  mutate(native = case_when(native >= 1 ~ 1,
                            native == 0 ~ 0))
final <- final %>%
  left_join(native, by = "Species") %>%
  filter(native == 1) %>%
  distinct() %>%
  mutate(Bioregion = case_when(is.na(Bioregion) & AsiaNative == 1 ~ "Unknown",
                               .default = Bioregion))

### which species are endemic to tropical asia
endemic <- final %>%
  group_by(Species) %>%
  reframe(AsiaEndemic = any(is.na(Bioregion)) | any(IncludesNonAsia == 1)) %>%
  mutate(AsiaEndemic = case_match(AsiaEndemic,
                                  FALSE ~ 1,
                                  TRUE  ~ 0))
final <- left_join(final, endemic, by = "Species")

### Filter only native intersections
final <- final %>%
  filter(AsiaNative == 1)

### convert to species x bioregion data frame
final <- final %>%
  select(Species, Genus, Family, Year, Bioregion, AsiaNative, AsiaEndemic) %>%
  filter(!duplicated(.))  %>%
  pivot_wider(id_cols = c(Species, Genus, Family, Year, AsiaEndemic),    
              names_from  = Bioregion,
              values_from = AsiaNative, 
              values_fill = 0) %>%
  relocate(Species, Genus, Family, Year, AsiaEndemic)

### total number of species native to region and number of species endemic to tropical asia
c(native = nrow(final),
  endemic = sum(final$AsiaEndemic))
colSums(final[, unique(bioregions$Bioregion[!is.na(bioregions$Bioregion)])])

### save results
write.csv(final, file.path(resDir, "Intersections", "Intersections_bioregions_trichoptera.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
