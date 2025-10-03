####################################################################################################
### Extracts description dates from the AntCat API
### I'm sure there is a much quicker way but this does the job
###
### Charlie Marsh
### charliem2003@github
### 06/2024
###
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Necessary libraries
library(httr2)
library(jsonlite)
library(dplyr)
library(sf)

### Directory structures
projDir <- "Diversity"                                    # project dir
kassDir <- file.path("Dir", "with", "Kass", "rangemaps")  # dir where we saved the Kass range maps
funDir  <- file.path(projDir, "Analysis_functions")       # dir that contains the function scripts

### Function to isolate description year from author names
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))

### list of species that we created range maps for
spList <- st_read(file.path(kassDir, "Kass_alpha_hulls.gpkg"))
spList <- st_drop_geometry(spList)$Species
spList <- gsub("_", " ", spList)

#==================================================================================================#
#--------------------------------------- Species-level info ---------------------------------------#
#==================================================================================================#

### results storage - 3619 species
all <- tibble(Species         = spList,
              level           = NA,
              id              = NA,
              status          = NA,
              name_id         = NA,
              name_cache      = NA,
              author_citation = NA)

### Loop through species and extract author info
pb <- txtProgressBar(1, nrow(all), 0, style = 3)
for(i in 1:nrow(all)) {
  setTxtProgressBar(pb, i)
  
  ### species
  sp <- all$Species[i]
  
  ### get species ID
  spID <- request(paste0("https://antcat.org/v1/taxa/search/",
                         paste(unlist(strsplit(sp, " ")), collapse = "%20"))) |>
    req_url_path(paste0("v1/taxa/search/",
                        paste(unlist(strsplit(sp, " ")), collapse = "%20"))) |>
    req_perform() |>
    resp_body_json()
  if(length(spID) > 0) {
    spID <- spID[[1]]$id
    
    ### extract author info
    spInfo <- request(paste0("https://antcat.org/v1/taxa/", spID)) |>
      req_url_path(paste0("v1/taxa/", spID)) |>
      req_perform() |>
      resp_body_json()
    
    all[all$Species == sp, ] <- data.frame(Species         = gsub(" ", "_", sp),
                                           Level           = names(spInfo),
                                           ID              = spInfo[[1]]$id,
                                           Status          = spInfo[[1]]$status,
                                           name_id         = spInfo[[1]]$name_id,
                                           name_cache      = spInfo[[1]]$name_cache,
                                           author_citation = spInfo[[1]]$author_citation)
  }
}

#==================================================================================================#
#-------------------------------------- Tidying up and save ---------------------------------------#
#==================================================================================================#

### First extract genera
all$Species <- gsub(" ", "_", all$Species)
all$Genus <- sapply(all$Species, function(x) unlist(strsplit(x, "_"))[[1]])

### All ant species are in the same family
all$Family <- "Formicidae"

### Extract description year from author information
all <- all %>%
  mutate(Year = clean_publication_dates(author_citation))

### There are some species in the Kass range maps that are no longer in antcat. Add manually
all$Year[all$Species == "Hypoponera_exoecata"]     <- 1928 # https://antwiki.org/wiki/Hypoponera_exoecata
all$Year[all$Species == "Myrmoteras_cuneonodus"]   <- 1998 # https://www.antwiki.org/wiki/Myrmoteras_cuneonodus
all$Year[all$Species == "Polyrhachis_glabrinotum"] <- 1930 # https://www.antwiki.org/wiki/Polyrhachis_glabrinotum

### Some fundamental checks
str(all)
all(is.numeric(all$Year))

### Tidying up
all <- all %>%
  rename(Level = level, ID = id, Status = status) %>%
  select(Species, Genus, Family, Level, ID, Status, Year)

### Save csv
write.csv(all, file.path(kassDir, "AntCat_taxonomy_for_BTAS_species.csv"),
                         row.names = FALSE, quote = FALSE)


