# CALCULATE ENDEMICITY ============
# Generate

rm(list = ls())

## PACKAGES ====================
library(plyr)
library(tidyverse)
library(reshape2)

## DIRECTORIES ====================

# Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
proj.dir <- "Diversity_and_endemism"

# You shouldn't need to adjust these folders
data.dir <- file.path(proj.dir, "Intersections", "Intersections")
fun.dir  <- file.path(proj.dir, "Analysis_functions")
res.dir  <- file.path(proj.dir, "Endemicity", "Results")
if(!dir.exists(res.dir)) { dir.create(res.dir, recursive = TRUE) }

# ==================================================================================================
## IMPORT DATA ====================
# ==================================================================================================

# Function for reading in and processing intersection spreadsheet
source(file.path(fun.dir, "Intersections", "processIntersection.R"))

### Read in bioregion info which has info on labelling
bioregion_info <- read.csv(file.path(main.dir, "Bioregion_info.csv"))

# Spreadsheet with summary info on each taxon
taxonInfo <- read.csv(file.path(main.dir, "Global_totals.csv")) %>%
  mutate(Realm = apply(., 1, function(x) trimws(paste(x["Group"], x["Taxon_Higher"])))) %>%
  arrange(desc(Group), desc(Taxon_Higher), Taxon) %>%
  mutate(Order = 1:nrow(.)) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon[order(Order)]))

# List of taxa
taxa <- taxonInfo$Taxon

# Loop through and read in data [uses Global_totals spreadsheet to control taxa read in]
combined_intersect <- tibble()
for(taxon in taxa) {
  taxonPA <- processIntersection(taxon = taxon, taxonInfo = taxonInfo, dataDir = data.dir)
  combined_intersect <- bind_rows(combined_intersect,
                                  taxonPA$mean)
}

# ==================================================================================================
## DIVERSITY FOR TROPICAL ASIA ====================
# ==================================================================================================

# Pivot to match formatting of previous code
combined_distr_df <- combined_intersect %>%
  pivot_longer(cols = -any_of(c("Taxon", "Group", "Taxon_Higher", "Global_total",
                                "Species", "Genus", "Family", "Year", "AsiaEndemic")),
               names_to  = "Region",
               values_to = "Present") %>%
  filter(Present == 1)

# Get correct labelling from bioregion_info
combined_distr_df <- left_join(combined_distr_df,
                               select(bioregion_info, Bioregion, Label),
                               join_by(Region == Bioregion)) %>%
  relocate(Label, .after = Region)

write.csv(combined_distr_df,
          file.path(res.dir, "taxon_full_list.csv"),
          row.names = FALSE)

# Calculate diversity at taxon level
combined_distr_all_df <- ddply(.data = combined_distr_df,
                               .variables = .(Group, Taxon_Higher, Taxon),
                               .fun = function(x){
                                 unique(x[,c("Species", "AsiaEndemic")])
                               })

taxon_div_summary <- ddply(.data = combined_distr_all_df,
                           .variables = .(Group, Taxon_Higher, Taxon),
                           .fun = summarise,
                           n_tropasia_species  = length(Species),
                           n_tropasia_endemics = sum(AsiaEndemic >= 1)) %>%
  mutate(Taxon = factor(Taxon, levels = taxa)) %>%
  arrange(Taxon)

# Save
write.csv(taxon_div_summary,
          file.path(res.dir, "taxon_div_summary.csv"),
          row.names = FALSE)

# ==================================================================================================
## DIVERSITY FOR BIOREGIONS ====================
# ==================================================================================================

### Loop through terrestrial/freshwater and marine
realms <- unique(taxonInfo$Realm)

### loop through mean, lower and upper estimates
measures <- c("mean", "lower", "upper")

for(measure in measures) {
  for(realm in realms) {
    print(paste(measure, realm))
    ### Only taxa relevant to realm
    taxonInfoGroup <- filter(taxonInfo, Realm == realm)
    taxa <- taxonInfoGroup$Taxon
    
    combined_intersect <- tibble()
    for(taxon in taxa) {
      ### Only proceeed with lower/upper estimate if taxon used range maps
      taxonType <- taxonInfoGroup$Type[taxonInfoGroup$Taxon == taxon]
      if(measure == "mean" | (measure %in% c("lower", "upper") & taxonType == "Range maps")) {
        taxonPA <- processIntersection(taxon = taxon, taxonInfo = taxonInfoGroup, dataDir = data.dir)
        combined_intersect <- bind_rows(combined_intersect,
                                        taxonPA[[measure]])
      }
    }
    
    if(nrow(combined_intersect) > 0) {
      # Pivot to match formatting of previous code
      combined_distr_df <- combined_intersect %>%
        pivot_longer(cols = -any_of(c("Taxon", "Group", "Taxon_Higher", "Global_total",
                                      "Species", "Genus", "Family", "Year", "AsiaEndemic")),
                     names_to  = "Region",
                     values_to = "Present") %>%
        filter(Present == 1)
      
      # Remove 'Unknown' region cases (but should be included in the tropical Asia species count above)
      combined_distr_df <- filter(combined_distr_df, Region != "Unknown")
      
      # Calculate how many bioregions each species occurs in 
      species_nareas_df <- combined_distr_df %>%
        ddply(.variables = .(Taxon, Species),# .(Taxon, Group, Taxon_Higher, Species, AsiaEndemic),
              .fun = summarise, 
              nareas = length(unique(Region)),
              .progress = "text")
      
      # Calculate endemicity in each region
      taxa_endemic_df <- combined_distr_df %>%
        left_join(species_nareas_df, by = join_by(Species, Taxon)) %>%
        ddply(.variables = .(Taxon, Group, Taxon_Higher, Region),
              .fun = summarise, 
              n_species  = length(Species),
              n_endemics = length(Species[nareas == 1 & AsiaEndemic == 1]),
              endemicity = n_endemics / n_species ) %>%
        mutate(label = paste0(round(endemicity, 2)," (", n_endemics, "/",n_species,")"))
      
      # Get correct labelling from bioregion_info
      bioregion_info$Bioregion[bioregion_info$Realm == "Marine"] <- gsub(" ", "_",  bioregion_info$Bioregion[bioregion_info$Realm == "Marine"])
      taxa_endemic_df <- left_join(taxa_endemic_df,
                                   select(bioregion_info, Bioregion, Label),
                                   join_by(Region == Bioregion)) %>%
        relocate(Label, .after = Region)
      
      # Prepare region by taxon endemicity table
      taxa_endemic_mat <- acast(Taxon ~ Label, fill = "", value.var = "label",
                                data = taxa_endemic_df)

      # Save
      write.csv(taxa_endemic_df,
                file.path(res.dir, paste0("taxon_region_endemicity_long_df", "_",
                                          tolower(gsub(" ", "_", realm)), "_", measure, ".csv")),
                row.names = FALSE)
      
      write.csv(taxa_endemic_mat,
                file.path(res.dir, paste0("taxon_region_endemicity_wide_df", "_",
                                          tolower(gsub(" ", "_", realm)), "_", measure, ".csv")),
                row.names = TRUE)
    }
  }
}
