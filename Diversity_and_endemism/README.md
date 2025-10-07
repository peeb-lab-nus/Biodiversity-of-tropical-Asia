## Analysis and results scripts for the Biogeography of Tropical Asia project
## - Diversity and endemism analyses

Scripts to carry out analyses for the 'Diversity and endemism' section of Lim et al. 'The rich, 
unique, and threatened biodiversity of tropical Asia'. First, we generate checklists of species 
present in each subregion, then these are summarised to create the tables and figures in the 
manuscript.

### 1. Species intersections with tropical Asia and subregions

Processes the biodiversity data (rangemaps, checklists or occurrence data) to generate a species 
list for each subregion. Each taxon has a separate script (note, some groups may take a week or 
more to run depending on how many threads you have).

As we don't own the underlying data, in each case you will need to access the data sources from the 
data provider outlined in the Supplementary Material. This may also be the case for the taxonomic 
information if the spatial data does not have all the relevant information (e.g. the description 
year). You will also need some other spatial data sources, such as GADM 4.1 for coastline data or 
hydroBASINs data from the IUCN.

Important: to fully replicate the exact results you will need the same versions as the ones outlined 
at the top of each script. Using more recent versions may result in differing results.

Input data required (depends on taxon):

- various spatial data (e.g. shapefiles with the subregions) stored in folder 'Data'

- the biodiversity data source (outlined at top of each script)

- if necessary, a separate taxonomy information source (for description dates and higher-order taxonomic information)

- You will need to create a version of GADM called 'GADM_410_land_Equal_Area.gpkg'. This is simply GADM that has
been unioned (i.e. all land borders removed) and preprojected to that of the rest of the project:
World Cylindrical Equal Area (+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs)
This is to act as a mask on top of the range maps to remove parts of range that fall in ocean areas.
You may need to do some manual cleaning within New Guinea and a few other locations where the internal 
borders within GADM have not been aligned properly

- HydroBASINS data (2023/08 version) available on the IUCN website. You will need levels 08, 10 and 
12 - see https://www.iucnredlist.org/resources/spatial-data-download

- The GBIF taxonomic backbone - available at https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c

- WoRMS taxonomic data - request here https://www.marinespecies.org/usersrequest.php

### 2. Summarise diversity and endemism in tropical Asia and subregions



