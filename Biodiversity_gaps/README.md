---
editor_options: 
  markdown: 
    wrap: 100
---

## Analysis and results scripts for the Biogeography of Tropical Asia project

## - Biodiversity knowledge gaps analyses

Scripts to carry out analyses for the 'Gaps in biodiversity knowledge' section of Lim et al. 'The
rich, unique, and threatened biodiversity of tropical Asia'. Takes the occurrence data from GBIF and
OBIS, removes points outside the region, with incomplete data, are spatio-temporal duplicates and
other issues, and summarises the remaining data for each taxonomic group in 100 x 100 km grid cells.

*NOTE: for running the 'Rates of species discovery' analyses see the 'Diversity_and_endemism'
folder*

### 1. Clean occurrence data for each taxonomic group

Filters out the points for a each taxonomic group from the global GBIF and OBIS .parquet data dumps,
and applies the cleaning steps outlined in Table S4. Scripts for each taxonomic group are available
in the folder 'Occurrences/cleaning'. It will output the cleaned data to the folder
'Occurrences/cleaning/cleaned'.

NOTE: the global GBIF data dump is \>250GB. For some taxonomic groups you will need considerable RAM
(128GB+) to run the cleaning steps. Even then, the birds need to be processed in smaller chunks. If
you are struggling to process smaller groups, see the birds cleaning script for how to adapt the
script to process chunks instead.

Input data required (depends on taxon):

-   for terrestrial groups, the .parquet file of the entire GBIF data dump
    (<https://www.gbif.org/occurrence-snapshots>). For this study we use the 01/10/2024 snapshot.

-   for marine groups, the .parquet file of the entire OBIS data dump
    (<https://obis.org/data/access/>). For this study we use the 23/07/2024 snapshot.

-   The cleaning steps to be applied read in from the csv files stored in 'Occurrences/cleaning' -
    'cleaningSteps_BTAS_10_2024.csv' for terrestrial groups and 'cleaningSteps_OBIS_12_2024.csv' for
    marine groups (the steps are the same, but for terrestrial groups points falling outside of land
    areas are retained and non-land areas removed, and vice-versa for marine groups). *The path to
    the GADM and the BTAS grid shapefiles may need to be changed for your directory structure*.

-   The equal-area grid for the BTAS region (BTAS_grid_100km.gpkg), available in the 'Data' folder.

-   A (slightly cleaned) version of GADM 4.1 for the BTAS region '(GADM_410_land_Equal_Area.gpkg)'
    which defines land and non-land areas, available in the 'Data' folder.

### 2. Summarise occurrence data for 100 x 100 km grid cells (Spatial_coverage.R)

Generate the various metrics described in the Supplementary Material. Figures are saved to the
folder 'Figures'. Results tables are saved to the folder 'Results'.
