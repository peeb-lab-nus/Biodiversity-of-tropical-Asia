####################################################################################################
### 
### Charlie Marsh
### charliem2003@github
### 12/2024
###
### Summarises cleaned occurrence data across 100 x 100 km grid cells
###
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(patchwork)
library(lubridate)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir <- "Biodiversity_gaps"                                        # project dir

### You shouldn't need to adjust these folders
dataDir <- file.path(projDir, "Data")                                 # dir that contains the function scripts
occDir  <- file.path(projDir, "Occurrences", "cleaning", "cleaned")   # dir with cleaned point data
funDir  <- file.path(projDir, "Analysis_functions")                   # dir with helper functions
resDir  <- file.path(projDir, "Results")                              # dir to save results to
figDir  <- file.path(projDir, "Figures")                              # dir to save figures to

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Function for calculating sampling coverage metrics in 100 x 100 km grid cells
source(file.path(funDir, "Point_plotting", "extractMetrics.R"))

### Various plotting helper functions
source(file.path(funDir, "Point_plotting", "samplingCoveragePlottingFunctions.R"))

### Spreadsheet with taxon information
taxonInfo <- read_csv(file.path(projDir, "Global_totals.csv")) %>%
  mutate(Taxon_Group = apply(., 1, function(x) paste(x["Group"], x["Taxon_Higher"]))) %>%
  mutate(Taxon_Group = gsub(" NA", "", Taxon_Group)) %>%
  arrange(desc(Taxon_Group))

### Groups to aggregate to
group_list <- unique(taxonInfo$Taxon_Group)

### The BTAS grid
btas <- st_read(file.path(dataDir, "BTAS_grid_100km.gpkg")) %>%
  rename(geometry = geom)

# Land for the BTAS region
land <- st_read(file.path(dataDir, "GADM_410_BTAS_land.gpkg")) %>%
  # st_simplify(dTolerance = 2500) %>%
  st_simplify(dTolerance = 5000)

### More constrained plotting region used in MS
msExt <- st_bbox(c(xmin = 7614253, xmax = 17232260, ymin = -1788714, ymax = 3279685),
                 crs = st_crs(btas))
btas <- st_crop(btas, msExt)
land <- st_crop(land, msExt)

### Convert to wgs84 for filtering occurrence records
msExt <- st_transform(msExt, 4326)

#==================================================================================================#
#--------------------------------- Summary tables of data quantity --------------------------------#
#==================================================================================================#

### Summarise all taxa
sumTable <- tibble()

taxa <- taxonInfo
taxa <- unique(taxa$Name)

for(taxon in taxa) {
  print(taxon)

  ### Obis or gbif?
  taxonSource <- tolower(taxonInfo$Occurrence_database[taxonInfo$Name == taxon])

  ### Read in cleaned point data for taxon and intersect with BTAS grid
  pts <- read_csv(file.path(occDir, paste0(taxonSource, "_filtered_", taxon, ".csv")),
                  show_col_types = FALSE, progress = FALSE)

  ### Filter for bbox used in manuscript
  pts <- pts %>%
    filter(if_all(matches("decimallongitude", ignore.case = TRUE), ~ .x >= msExt$xmin)) %>%
    filter(if_all(matches("decimallongitude", ignore.case = TRUE), ~ .x <= msExt$xmax)) %>%
    filter(if_all(matches("decimallatitude",  ignore.case = TRUE), ~ .x >= msExt$ymin)) %>%
    filter(if_all(matches("decimallatitude",  ignore.case = TRUE), ~ .x <= msExt$ymax))

  ### Is record from iNaturalist or eBIRD
  if(taxonSource == "gbif") {
    pts <- pts %>%
      mutate(citizen = case_when(!is.na(institutioncode) & institutioncode  == "iNaturalist" & yearClean >= 2008 ~ 1,
                                 !is.na(collectioncode)  & grepl("EBIRD", collectioncode)    & yearClean >= 2002 ~ 1,
                                 .default = 0))
  }
  if(taxonSource == "obis") {
    pts <- pts %>%
      mutate(citizen = case_when(!is.na(institutionCode) & institutionCode  == "iNaturalist" & yearClean >= 2008 ~ 1,
                                 !is.na(collectionCode)  & grepl("EBIRD", collectionCode)    & yearClean >= 2002 ~ 1,
                                 .default = 0))
  }

  sumTaxon <- data.frame(group = taxonInfo$Taxon_Group[taxonInfo$Name == taxon],
                         taxon = taxonInfo$Taxon[taxonInfo$Name == taxon],
                         name  = taxonInfo$Name[taxonInfo$Name == taxon],
                         nPts_clean    = nrow(pts),
                         nPts_specimen = nrow(pts) - sum(pts$citizen),
                         nPts_citizen  = sum(pts$citizen))
  sumTable <- bind_rows(sumTable, sumTaxon)

  rm(pts)
  gc()
}

### Join with spreadsheet containing pre-cleaning points
summary_df <- read_csv(file.path(resDir, "gbif_summary.csv")) %>%
  select(name, nPts_dataset, nPts_global)
sumTable <- sumTable %>%
  left_join(summary_df, by = "name") %>%
  relocate(nPts_dataset, nPts_global, .after = name)

sumTable$perc_global_pts     <- (sumTable$nPts_clean / sumTable$nPts_global) * 100
sumTable$perc_all_clean_pts  <- (sumTable$nPts_clean / sum(sumTable$nPts_clean)) * 100
sumTable$perc_taxon_specimen <- (sumTable$nPts_specimen / sumTable$nPts_clean) * 100
sumTable$perc_taxon_citizen  <- (sumTable$nPts_citizen  / sumTable$nPts_clean) * 100
print(sumTable, n = 24)

write_csv(sumTable, file.path(resDir, "cleaning_summary.csv"))

#==================================================================================================#
#--------------------------- All three measures - split up by data type ---------------------------#
#==================================================================================================#

### Append all group maps together
allMap <- tibble()

for(i in 1:length(group_list)) {
  infoGrp <- taxonInfo %>%
    filter(Taxon_Group == group_list[i])
  taxaGrp <- infoGrp$Name

  for(taxon in taxaGrp) {
    print(taxon)

    ### Obis or gbif?
    taxonSource <- tolower(infoGrp$Occurrence_database[infoGrp$Name == taxon])

    ### Read in cleaned point data for taxon
    pts <- read_csv(file.path(occDir, paste0(taxonSource, "_filtered_", taxon, ".csv")),
                    show_col_types = FALSE, progress = FALSE)

    ### Filter for bbox used in manuscript
    pts <- pts %>%
      filter(if_all(matches("decimallongitude", ignore.case = TRUE), ~ .x >= msExt$xmin)) %>%
      filter(if_all(matches("decimallongitude", ignore.case = TRUE), ~ .x <= msExt$xmax)) %>%
      filter(if_all(matches("decimallatitude",  ignore.case = TRUE), ~ .x >= msExt$ymin)) %>%
      filter(if_all(matches("decimallatitude",  ignore.case = TRUE), ~ .x <= msExt$ymax))

    if(taxonSource == "gbif") {
      pts <- pts %>%
        mutate(citizen = case_when(!is.na(institutioncode) & institutioncode  == "iNaturalist" & yearClean >= 2008 ~ 1,
                                   !is.na(collectioncode)  & grepl("EBIRD", collectioncode)    & yearClean >= 2002 ~ 1,
                                   .default = 0)) %>%
        st_as_sf(coords = c("decimallongitude", "decimallatitude"), crs =  4326)
    }
    if(taxonSource == "obis") {
      pts <- pts   %>%
        mutate(citizen = case_when(!is.na(institutionCode) & institutionCode  == "iNaturalist" & yearClean >= 2008 ~ 1,
                                   !is.na(collectionCode)  & grepl("EBIRD", collectionCode)    & yearClean >= 2002 ~ 1,
                                   .default = 0))%>%
        st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs =  4326)
    }

    ##############################################################
    ### All points

    ### Intersect with BTAS grid
    ptsAll <- pts %>%
      select(species, yearClean) %>%
      st_transform(st_crs(btas)) %>%
      st_join(btas, join = st_within)

    ### Extract metrics
    taxonMap <- extractMetrics(ptsAll) %>%
      mutate(points = "all")

    ### bind onto global dataframe
    allMap <- bind_rows(allMap, taxonMap)

    ##############################################################
    ### Specimen data

    ### Intersect with BTAS grid
    ptsSpecimen <- pts %>%
      filter(citizen == 0) %>%
      select(species, yearClean) %>%
      st_transform(st_crs(btas)) %>%
      st_join(btas, join = st_within)

    ### Extract metrics
    taxonMap <- extractMetrics(ptsSpecimen) %>%
      mutate(points = "specimen")

    ### bind onto global dataframe
    allMap <- bind_rows(allMap, taxonMap)

    ##############################################################
    ### Citizen science data

    ### Intersect with BTAS grid
    ptsCitizen <- pts %>%
      filter(citizen == 1) %>%
      select(species, yearClean) %>%
      st_transform(st_crs(btas)) %>%
      st_join(btas, join = st_within)

    ### Extract metrics
    taxonMap <- extractMetrics(ptsCitizen) %>%
      mutate(points = "citizen")

    ### bind onto global dataframe
    allMap <- bind_rows(allMap, taxonMap)

    rm(pts)
    rm(ptsAll)
    rm(ptsSpecimen)
    rm(ptsCitizen)
    gc()
  }
}

### Save
write_csv(allMap, file.path(resDir, "all_results.csv"))

#==================================================================================================#
#-------------------------------------------- Plotting --------------------------------------------#
#==================================================================================================#

### Calculate median values for each group for each metric per cell
meds <- read_csv(file.path(resDir, "all_results.csv")) %>%
  group_by(id_100km, metric, group, points) %>%
  reframe(value = median(value))

### For plotting merge back with geometry
meds <- left_join(btas, meds, by = "id_100km", relationship = "many-to-many")

### Arrangement to get right ordering of groups and plots
meds <- meds %>%
  mutate(metric = case_when(metric == "Standardised number of occurrences" ~ "A", # "A. Number of occurrence records",
                            metric == "Temporal replicates" ~ "B",                # "B. Repeated sampling through time",
                            metric == "Median record age" ~ "C")) %>%             # "C. Average record age")) %>%
  mutate(metric = factor(metric,
                         levels = c("A",   # "A. Number of occurrence records",
                                    "B",   # "B. Repeated sampling through time", 
                                    "C"))) # "C. Average record age")))

meds <- meds %>%
  mutate(group = case_when(group == "Vascular plants"                          ~ "Vascular\nplants",
                           group == "Terrestrial and freshwater vertebrates"   ~ "Terrestrial &\nfreshwater\nvertebrates",
                           group == "Terrestrial and freshwater invertebrates" ~ "Terrestrial &\nfreshwater\ninvertebrates",
                           .default = group)) %>%
  mutate(group = factor(group,
                        levels = c("Vascular\nplants", "Terrestrial &\nfreshwater\nvertebrates",
                                   "Terrestrial &\nfreshwater\ninvertebrates", "Marine")))

#=======================
### All points
#=======================

### Split into different metrics - each metric will become a column (with it's own scale bar)
medsSplit <- meds %>%
  filter(points == "all") %>%
  group_split(metric)

### Generate plot and save
pAll <- panelPlot(medsSplit, facet = "group")

ggsave(file.path(figDir, "all_records.png"), pAll,
       height = 143, width = 180, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(figDir, "all_records.pdf"), pAll,
       height = 143, width = 180, units = "mm",  dpi = 600, bg = "white", cairo_pdf)

#=======================
### Citizen science only
#=======================

### Split into different metrics - each metric will become a column (with it's own scale bar)
medsSplit <- meds %>%
  filter(points == "citizen") %>%
  group_split(metric)

### Generate plot and save
pCit <- panelPlot(medsSplit, facet = "group")

ggsave(file.path(figDir, "citizen_records.png"), pCit,
       height = 110, width = 180, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(figDir, "citizen_records.pdf"), pCit,
       height = 110, width = 180, units = "mm",  dpi = 600, bg = "white", cairo_pdf)

#=======================
### Non-citizen science records
#=======================

### Split into different metrics - each metric will become a column (with it's own scale bar)
medsSplit <- meds %>%
  filter(points == "specimen") %>%
  group_split(metric)

### Generate plot and save
pSpec <- panelPlot(medsSplit, facet = "group")

ggsave(file.path(figDir, "specimen_records.png"), pSpec,
       height = 143, width = 180, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(figDir, "specimen_records.pdf"), pSpec,
       height = 143, width = 180, units = "mm",  dpi = 600, bg = "white", cairo_pdf)

#==================================================================================================#
#---------------------------- Plotting each taxon plotted individually ----------------------------#
#==================================================================================================#

### Read in data and just keep the all points data
meds <- read_csv(file.path(resDir, "all_results.csv")) %>%
  filter(points == "all")

### For plotting merge back with geometry
meds <- left_join(btas, meds, by = "id_100km", relationship = "many-to-many")

### Arrangement to get right ordering of groups and plots
meds <- meds %>%
  mutate(metric = case_when(metric == "Standardised number of occurrences" ~ "A", # "A. Number of occurrence records",
                            metric == "Temporal replicates" ~ "B",                # "B. Repeated sampling through time",
                            metric == "Median record age" ~ "C")) %>%             # "C. Average record age")) %>%
  mutate(metric = factor(metric,
                         levels = c("A",   # "A. Number of occurrence records",
                                    "B",   # "B. Repeated sampling through time", 
                                    "C"))) # "C. Average record age")))

meds <- meds %>%
  mutate(group = case_when(group == "Vascular plants"                          ~ "Vascular\nplants",
                           group == "Terrestrial and freshwater vertebrates"   ~ "Terrestrial &\nfreshwater\nvertebrates",
                           group == "Terrestrial and freshwater invertebrates" ~ "Terrestrial &\nfreshwater\ninvertebrates",
                           .default = group)) %>%
  mutate(group = factor(group,
                        levels = c("Vascular\nplants", "Terrestrial &\nfreshwater\nvertebrates",
                                   "Terrestrial &\nfreshwater\ninvertebrates", "Marine")))

for(taxa in c("Vascular\nplants", "Terrestrial &\nfreshwater\nvertebrates",
              "Terrestrial &\nfreshwater\ninvertebrates", "Marine")) {
  
  ### Split into different metrics - each metric will become a column (with it's own scale bar)
  medsSplit <- meds %>%
    filter(group == taxa) %>%
    group_split(metric, group)
  
  pTax <- panelPlot(medsSplit, facet = "taxon")
  
  ### plotting params
  if(length(unique(medsSplit[[1]]$taxon)) == 4)  { plotHeight <- 143 }
  if(length(unique(medsSplit[[1]]$taxon)) == 5)  { plotHeight <- 173 }
  if(length(unique(medsSplit[[1]]$taxon)) == 11) { plotHeight <- 323 }
  
  taxa <- gsub("\n", "_", taxa)
  taxa <- gsub(" " , "_", taxa)
  
  ggsave(file.path(figDir, paste0(taxa, "_records.png")), pTax,
         height = plotHeight, width = 180, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(figDir,  paste0(taxa, "_records.pdf")), pTax,
         height = plotHeight, width = 180, units = "mm",  dpi = 600, bg = "white", cairo_pdf)
}
