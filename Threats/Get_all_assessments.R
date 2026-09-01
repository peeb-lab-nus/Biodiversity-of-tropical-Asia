####################################################################################################
### Get IUCN asessments
### Charlie Marsh
### charliem2003@github
### 07/2026
###
### Gets IUCN RL assessments for all species using the IUCN API. To retrieve threats you will need
### an IUCN Red List API authentication token (https://api.iucnredlist.org). See the iunredlist
### package documentation for whether there are any OS-specific steps.
###
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(iucnredlist)
library(dplyr)
library(readr)
library(pbmcapply)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir  <- "Threats"                                  # project dir

### You shouldn't need to adjust these folders
interDir  <- file.path(projDir, "Intersections")       # dir with intersection results
assessDir <- file.path(projDir, "Species_assessments") # dir to save individual species IUCN assessments
resDir    <- file.path(projDir, "Results")             # dir to save collated results to

### Your IUCN api token (get from https://api.iucnredlist.org/users/edit)
IUCN_REDLIST_KEY <- "YOUR_API_TOKEN"
api <- init_api(IUCN_REDLIST_KEY)

#==================================================================================================#
#-------------------------------------- Read in summary data --------------------------------------#
#==================================================================================================#

### Taxon group information
taxonInfo  <- read_csv(file.path(projDir, "Global_totals.csv")) %>%
  filter(Taxon != "Sponges") %>%
  select(Taxon, Group, Taxon_Higher, Type, Name, Colour) %>%
  mutate(Group = gsub(" NA", "", paste(Group, Taxon_Higher, sep = " "))) %>%
  mutate(GroupForPlot = gsub(" ", "\n", Group, fixed = TRUE)) %>%
  mutate(GroupForPlot = gsub("and\n", "& ", GroupForPlot, fixed = TRUE)) %>%
  mutate(TaxonForPlot = case_when(Name == "phasmids" ~ "Stick and leaf insects",
                                  Name == "sharks"   ~ "Sharks and rays",
                                  .default = Taxon))

### Bioregion information
regionInfo <- read_csv("Bioregion_info.csv") %>%
  mutate(LabelForPlot = gsub(" ", "\n", Label, fixed = TRUE))

#==================================================================================================#
#----------------------------------- Load in BTAS intersections -----------------------------------#
#==================================================================================================#

### List of taxa
taxa <- taxonInfo$Name

### Read in intersections and determine if native to tropical Asia
intersections <- tibble()
for(taxon in taxa) {
  print(taxon)
  group <- read_csv(file.path(interDir, paste0("Intersections_bioregions_", taxon, ".csv")),
                    show_col_types = FALSE) %>%
    mutate(Name = taxon) %>%
    relocate(Name)
  intersections <- bind_rows(intersections, group)
}

### Get correct names and groups from Global_totals.csv
intersections <- intersections %>%
  left_join(select(taxonInfo, Name, Group, Name, GroupForPlot, TaxonForPlot), by = "Name") %>%
  relocate(Name, TaxonForPlot, Group, GroupForPlot, Species, Genus, Family, Year, AsiaEndemic, Range_area)

### 147,275 species in total (some fish species shared between bony fish and freshwater fish)
length(unique(intersections$Species))
reframe(intersections, n = n(), .by = TaxonForPlot) %>% t()

### Two flies share names with plants, but these flies don't feature in IUCN
intersections %>%
  filter(Group != "Marine") %>%
  filter(Species %in% Species[duplicated(Species)]) %>%
  arrange(Species)

### Alter fly names (they have no IUCN assessments) to avoid them being assigned the plant assessment
intersections <- intersections %>%
  mutate(Species = case_when(Name == "diptera" & Species == "Ormosia_formosana" ~ "Ormosia_formosana-diptera",
                             Name == "diptera" & Species == "Limnophila_glabra" ~ "Limnophila_glabra-diptera",
                             .default = Species))

####################################################################################################
### Process IUCN for name checking
####################################################################################################

## Clean up names for IUCN matching
all_intersections <- intersections %>%
  mutate(Species       = gsub(" ", "_", Species)) %>%
  mutate(genus_split   = sapply(Species, function(x) strsplit(x, "_")[[1]][1])) %>%
  mutate(GenusIUCN     = case_when(!is.na(Genus) ~ Genus,
                                   is.na(Genus) ~ genus_split)) %>%
  mutate(GenusIUCN     = gsub("[", "", GenusIUCN, fixed = TRUE)) %>%
  mutate(GenusIUCN     = gsub("]", "", GenusIUCN, fixed = TRUE)) %>%
  mutate(species_split = gsub(" ", "_", Species)) %>%
  mutate(species_split = gsub("_sensu_lato_", "_", species_split)) %>%
  mutate(species_split = sapply(species_split, function(x) {
    if(grepl("[(]", x)) {
      gen <- strsplit(x, "[_(]")[[1]][1]
      sp  <- strsplit(x, "[)_]")[[1]][4]
      return(paste(gen, sp, sep = "_"))
    } else {
      return(x)
    }
  })) %>%
  mutate(species_split = gsub(" ", "_", species_split)) %>%
  mutate(SpeciesIUCN   = sapply(species_split, function(x) strsplit(x, "_")[[1]][2])) %>%
  mutate(SpeciesIUCN   = gsub("[!]", "", SpeciesIUCN, fixed = TRUE)) %>%
  mutate(NameIUCN      = paste(GenusIUCN, SpeciesIUCN, sep = "_"))

### Where there are lumps, join together
all_intersections <- all_intersections %>%
  reframe(Year                         = min(Year),
          AsiaEndemic                  = min(AsiaEndemic),
          Range_area                   = sum(Range_area, na.rm = TRUE),
          Indian_Subcontinent          = max(Indian_Subcontinent),
          IndoChina                    = max(IndoChina),
          Philippines                  = max(Philippines),
          Malaya                       = max(Malaya),
          Sumatra                      = max(Sumatra),
          Java                         = max(Java),
          Borneo                       = max(Borneo),
          Sulawesi                     = max(Sulawesi),
          Lesser_Sundas                = max(Lesser_Sundas),
          Maluku                       = max(Maluku),
          New_Guinea                   = max(New_Guinea),
          Unknown                      = max(Unknown),
          Andaman                      = max(Andaman),                 
          Sahul_Shelf                  = max(Sahul_Shelf),
          Western_Coral_Triangle       = max(Western_Coral_Triangle),
          Eastern_Coral_Triangle       = max(Eastern_Coral_Triangle),
          Central_Indian_Ocean_Islands = max(Central_Indian_Ocean_Islands),
          Java_Transitional            = max(Java_Transitional),
          Bay_of_Bengal                = max(Bay_of_Bengal),
          Sunda_Shelf                  = max(Sunda_Shelf),
          South_China_Sea              = max(South_China_Sea),
          West_and_South_Indian_Shelf  = max(West_and_South_Indian_Shelf),
          South_Kuroshio               = max(South_Kuroshio),
          .by = c(Name, TaxonForPlot, Group, GroupForPlot, NameIUCN, GenusIUCN, SpeciesIUCN, Family)) %>%
  mutate(Unknown = case_when(Unknown == 1 & rowSums(across("Indian_Subcontinent":"New_Guinea"), na.rm = TRUE) > 0 ~ 0,
                             is.na(Unknown) ~ 0,
                             .default = Unknown))

### 147,746 sppecies from 147,253 species
nrow(all_intersections)
length(unique(all_intersections$NameIUCN))

### Save
all_intersections <- write_csv(all_intersections, file.path(resDir, "all_intersections.csv"))

####################################################################################################
### Retrieve IUCN assessments

### Function for retrieving IUCN assessment information from API
retrieveIUCN <- function(species, assessDir){
  # Retrieves IUCN RL status and threat scores
  #
  # Arguments
  #   x: name of species
  #
  # Returns
  #   data.frame of results
  # print(paste0("Retrieving IUCN species threats for ", x))
  
  outputFile <- file.path(assessDir, paste0(species, ".csv"))
  
  ### Separate genus and species name
  genName <- strsplit(species, "_")[[1]][1]
  spName  <- strsplit(species, "_")[[1]][2]
  
  if(!file.exists(outputFile)) {
    ### Get latest assessment (try 3 times)
    assessments <- NULL
    for(i in 1:3) {
      assessments <- try(iucnredlist::assessments_by_name(api     = api,
                                                          genus   = genName,
                                                          species = spName), silent = TRUE)
      if(any(class(assessments) == "data.frame")) { break }
      Sys.sleep(time = 0.25)
    }
    
    if(any(class(assessments) == "data.frame" & nrow(assessments) > 0)) {
      ### Filter just global asssements
      assessment <- assessments %>%
        filter(grepl("Global", scopes_description_en))
      
      if(nrow(assessment) > 0) {
        assessment <- assessment %>%
          mutate(year_published = as.numeric(year_published)) %>%
          arrange(-year_published) %>%
          select(taxon_scientific_name, year_published, latest, red_list_category_code, assessment_id)
      }
      
      ### Sometimes there are only local assessments
      if(nrow(assessment) == 0) {
        assessment <- data.frame(taxon_scientific_name  = species,
                                 year_published         = NA,
                                 latest                 = NA,
                                 red_list_category_code = NA,
                                 assessment_id          = NA)
      }
    }
    
    ### If species is not assessed
    if(any(class(assessments) == "try-error") |
       any(class(assessments) == "data.frame" & nrow(assessments) == 0)) {
      assessment <- data.frame(taxon_scientific_name  = species,
                               year_published         = NA,
                               latest                 = NA,
                               red_list_category_code = NA,
                               assessment_id          = NA)
    }
    
    write_csv(assessment, outputFile)
    Sys.sleep(0.5)
  }
}

### Get assessments (run in parallel if desired but be careful to avoid overloading the API)
allSp <- unique(all_intersections$NameIUCN)
assessments <- pbmclapply(allSp,
                          FUN = function(x) retrieveIUCN(species = x, assessDir = assessDir),
                          mc.cores = 1)

####################################################################################################
### Append together individual files

### Function to read in a tidy assessment data
processAssessment <- function(x) { 
  assessment <- read_csv(x, progress = FALSE, show_col_types = FALSE)
  
  ### Sometimes there are more than one assessment in a year (e.g. "Equisetum_arvense")
  if(nrow(assessment) > 1) {
    assessment <- assessment %>%
      mutate(year_published = as.numeric(year_published)) %>%
      arrange(-year_published, -assessment_id)
    
    if(any(duplicated(assessment$year_published))) {
      duplicateYears <- assessment$year_published[duplicated(assessment$year_published)]
      assessment <- assessment %>%
        filter(!year_published %in% duplicateYears |
                 year_published %in% duplicateYears & latest == TRUE)
    }
  }
  return(assessment)
}

### Process all assessments as list in parallel
allSp <- list.files(assessDir, pattern = ".csv", full.names = TRUE, recursive = FALSE)
all_assessments <- pbmclapply(allSp, processAssessment, mc.cores = 10)

### Unlist and bind together into single dataframe
all_assessments <- do.call("bind_rows", all_assessments)
all_assessments <- all_assessments %>%
  mutate(taxon_scientific_name = gsub(" ", "_", taxon_scientific_name)) %>%
  distinct()

### 147,253 species with 176,741 assessments
nrow(all_assessments)
length(unique(all_assessments$taxon_scientific_name))

### 147,253 IUCN species
length(unique(all_intersections$NameIUCN))
all_intersections$NameIUCN[!all_intersections$NameIUCN %in% all_assessments$taxon_scientific_name]
all_assessments$taxon_scientific_name[!all_assessments$taxon_scientific_name %in% all_intersections$NameIUCN]
filter(all_assessments, taxon_scientific_name %in% taxon_scientific_name[!taxon_scientific_name %in% all_intersections$NameIUCN])

### Check only one assessment per species per year
all_assessments %>% reframe(n = n(), .by = c(taxon_scientific_name, year_published)) %>% filter(n > 1)

### Save
write_csv(all_assessments, file.path(resDir, "all_assessments.csv"))
