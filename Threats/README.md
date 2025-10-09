## Analysis and results scripts for the Biogeography of Tropical Asia project

## - Threats analyses

Scripts to carry out analyses for the 'Threats' section of Lim et al. 'The rich, unique, and threatened biodiversity of tropical Asia'. First, we generate checklists of species present in each subregion, then these are summarised to create the tables and figures in the manuscript.

NOTES:

1)  For groups with IUCN range maps available, we need to run the species intersections with the subregions in the same manner as for the 'Diversity and endemism' section. The results will differ slightly with the results from that section when we have used a different range map source (e.g. we used MDD for mammals for the diversity and endemism analyses) or we have used a different version of the IUCN range maps.

2)  For groups without range maps (ie. checklist data, we will directly use the files outputted during the diversity and endemism analyses. Therefore, you must run those files (in the 'Diversity_and_endemism/Intersections' folder) first.

3)  The 'Rates of species discovery’ analyses are run as part of the previous section, 'Diversity_and_endemism'.

### 1. Species intersections with tropical Asia and subregions

Scripts inside the folder 'Intersections'. Scripts process the biodiversity data (IUCN rangemaps) to generate a species list for each subregion. Each taxon has a separate script (note, some groups may take a week or more to run depending on how many threads you have).

As we don't own the underlying data, in each case you will need to access the data sources from the data provider outlined in the Supplementary Material. This may also be the case for the taxonomic information if the spatial data does not have all the relevant information (e.g. the description year). You will also need some other spatial data sources, such as GADM 4.1 for coastline data or hydroBASINs data from the IUCN, outlined below. Other spatial data required are in the folder 'Data'.

*IMPORTANT: to fully replicate the exact results you will need the same versions as the ones outlined at the top of each script. Using more recent versions may result in differing results.*

Input data required (depends on taxon):

-   various spatial data (e.g. shapefiles with the subregions) stored in folder 'Data'.

-   the biodiversity data source (outlined at top of each script).

-   if necessary, a separate taxonomy information source for description dates and higher-order taxonomic information (outlined at top of each script).

-   You will need to create a version of GADM called 'GADM_410_land_Equal_Area.gpkg'. This is simply the global GADM that has been unioned (i.e. all land borders removed) and reprojected to that of the rest of the project: World Cylindrical Equal Area (+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs). This is to act as a mask on top of the range maps to remove parts of range that fall in ocean areas. You may need to do some manual cleaning within New Guinea and a few other locations where the internal borders within GADM have not been aligned properly.

-   HydroBASINS data (2023/08 version) available on the IUCN website. You will need levels 08, 10 and 12 - see <https://www.iucnredlist.org/resources/spatial-data-download>

-   The GBIF taxonomic backbone - available at <https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c>

### 2. Process intersections (Final_processing.R)

Once all the IUCN range map intersections have been run the script 'Intersections/Final_processing.R' will read in the intersections spreadsheets for all taxa and do final processing. For taxa without IUCN range maps it will read in the intersections spreadsheets generated in the 'Diversity_and_endemism' folder. Each taxon's spreadsheet is saved in the 'Intersections' folder.

There are a few final processing steps:

-   Merges Andamans into IndoChina.

-   Merges Palawan into Philippines.

-   Removes Northern Australia and species only found there within region.

-   Removes Northwest and Northeast Australian Shelves and species only found within those regions.

### 3. Generate look-up table to harmonise WCVP and IUCN names (get_threats_vascularplants.R)

This script will generate a spreadsheet to harmonise WCVP names (that used in our intersections) with the IUCN names (that used for the IUCN assessments).

### 4. Align the IUCN assessments with intersections (Process_intersections_and_iucn.R)

This script reads in the IUCN assessments and joins them with the intersections data. Will generate a file 'All_iucn_and_threats.csv' that has the IUCN Red List status and threat score (or NA if absent) for all species intersecting with tropical Asia.

It is important that you use the same assessment and threat year as those used for generating the assessments. Some care will be needed to make sure you have the same naming structure as used in this script, and adjust as necessary. For this project we used:

***IUCN taxonomies:***

IUCN 2025: plants, bony_fish, freshwater_fish, corals, all non-chordates

IUCN 2024: amphibians, freshwater_crabs, mammals, reptiles

BirdLife 2024: birds

WCVP 2025: plants

***IUCN assessments - save in folder with structure 'taxon_name/IUCN_v2024' etc:***

IUCN 2025: plants, bony_fish, freshwater_fish, corals, all non-chordates

IUCN 2024: amphibians, freshwater_crabs, mammals, reptiles

BirdLife 2024: birds

***IUCN threat scores:***

IUCN 2025: data for all species (see top of script

IUCN 2024: amphibians, birds, freshwater_crabs, mammals, reptiles

### 5. Summarise across tropical Asia and subregions and generate figures (New_threats.R)

Summarises the assessment and threat information across the region, for subregions and across taxonomic groups. Examines three facets:

1)  The proportion of species in each Red List status (or non-evaluated) for each taxon across tropical Asia, and for each taxonomic group in each subregion.

2)  The proportion of species in each taxon listed under different threats (threats are assigned to six broad categories).

3)  The year of the latest assessment (if the species is assessed) for each taxon and taxonomic group.

Figures are written to the folder 'Figures'.
