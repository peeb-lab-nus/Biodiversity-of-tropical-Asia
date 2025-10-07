####################################################################################################
### Code to regenerate the alpha-hull range maps used in Kass et al. 2022
### Charlie Marsh
### charliem2003@github
### 06/2024
###
### For each species the code produces a shapefile of the alpha-hull. Where a taxon is at a
### subspecies-level in Kass et al. we merge to species-level using the points from all subspecies
### under the same parent species.
###
### We follow the methods outlined in Kass et al 2022 using the cleaned occurrence points provided
### in the supplementary material of that paper, namely:
### 
### For species with ≥3 occurrence localities (9156 species)
###   - Alpha hulls were generated using the R package alphahull with alpha value 15.
###   - For 4 species we generate alpha shapes because of issues fitting alpha hulls
###   - In both cases we then buffered the polygon by 30 km
### For species with <3 occurrence localities (5168 species)
###   - We generated buffered points (30 km)
###
### Alpha hulls and points were buffered to account for spatial uncertainty. Points and alpha-hulls
### were buffered after projecting the data to the World Behrmann equal-area.
###
### In all cases, we masked areas corresponding to marine areas using GADM v4.1 (https://gadm.org/)
### and to inland water bodies sourced from natural earth (https://naturalearthdata.com).
###
### Citation:
### Kass, J.M., Guénard, B., Dudley, K.L., Jenkins, C.N., Azuma, F., Fisher, B.L., Parr, C.L.,
###   Gibb, H., Longino, J.T., Ward, P.S., Chao, A., Lubertazzi, D., Weiser, M., Jetz, W.,
###   Guralnick, R., Blatrix, R., Lauriers, J.D., Donoso, D.A., Georgiadis, C., Gomez, K.,
###   Hawkes, P.G., Johnson, R.A., Lattke, J.E., MacGown, J.A., Mackay, W., Robson, S.,
###   Sanders, N.J., Dunn, R.R., Economo, E.P., 2022. The global distribution of known and
###   undiscovered ant biodiversity. Science Advances 8, eabp9908.
###   https://doi.org/10.1126/sciadv.abp9908
### Kass, Jamie et al. (2022). The global distribution of known and undiscovered ant biodiversity
###   [Dataset]. Dryad. https://doi.org/10.5061/dryad.wstqjq2pp
###
####################################################################################################

### Load libraries
library(dplyr)
library(sf)
library(alphahull)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity"                                           # project dir
gabiDir   <- file.path("Dir", "with", "GABI", "occurrence_points") # dir containing clean GABI data
lakesDir  <- file.path("Dir", "with", "Natural_Earth", "Lakes")    # dir with Natural Earth lakes data
kassDir   <- file.path("Dir", "to", "save", "final", "rangemaps")  # dir to save final Kass range maps to

### You shouldn't need to adjust these folders
regionDir <- file.path(projDir, "Data")                            # dir that contains the subregions data
funDir    <- file.path(projDir, "Analysis_functions")              # dir that contains the function scripts
resDir    <- file.path(baseDir, "NUS", "BTAS", "Intersections")    # dir to save results to

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Intersections function
source(file.path(funDir, "Taxon_specific", "Ants", "Alpha_hull_functions.R"))

### regions data to check species intersects with tropical asia
regions <- st_read(file.path(regionDir, "Bioregions", "Bioregions_rangemaps.gpkg"))

### bounding box for tropical asia for quick exclusion of non-tropical asia species
regionsExtent <- st_transform(st_as_sfc(st_bbox(regions)), crs = 4326)

### global GADM for masking out marine areas
gadm <- st_read(file.path(gadmDir, "GADM_410_land_Equal_Area.gpkg"))

### lakes for masking out freshwater
lakes <- st_read(file.path(lakesDir, "ne_10m_lakes.shp")) %>%
  select(name) %>%
  st_transform(st_crs(gadm)) %>%
  st_make_valid()

### cleaned occurrence data
spPts <- read.csv(file.path(gabiDir, "processing_data", "3A_forAnalysis_processed_database_SPECIES.csv")) %>%
  as_tibble %>%
  mutate(valid_species_name = gsub("[.]", "_", valid_species_name)) %>%
  mutate(valid_species_name = gsub(" ", "_", valid_species_name))

### Merge subspecies into single species
spPts <- spPts %>%
  mutate(valid_species_name = sapply(valid_species_name,
                                     function(x) paste(strsplit(x, "_")[[1]][1:2], collapse = "_")))

### generate list of species to run
spList <- unique(spPts$valid_species_name)

done <- c(list.files(file.path(kassDir, "species_rangemaps"), pattern = ".gpkg", recursive = FALSE),
          list.files(file.path(kassDir, "species_rangemaps", "Non-native"), pattern = ".csv", recursive = FALSE))
done <- gsub(".csv", "", done)
done <- gsub(".gpkg", "", done)
spList <- spList[!spList %in% done]

#==================================================================================================#
#----------------------------------------- Generate MCPs ------------------------------------------#
#==================================================================================================#

### loop through species and append to overall shapefile
pb <- txtProgressBar(1, length(spList), style = 3)
for(i in 1:length(spList)) {
  setTxtProgressBar(pb, i)
  sp <- spList[i]
  
  ### extract species points, remove duplicates and convert to shapefile
  spPt <- spPts %>%
    filter(valid_species_name == sp) %>%
    filter(!(duplicated(lon_opt) & duplicated(lat_opt))) %>%
    st_as_sf(coords = c("lon_opt", "lat_opt"), crs = 4326) %>%
    st_make_valid() %>%
    select(valid_species_name) %>%
    mutate(valid_species_name = gsub(" ", "_", valid_species_name))
  sp <- gsub(" ", "_", sp)
  
  native <- TRUE
  if(nrow(spPt) > 1) {
    ### the buffer is there because some points are too close together to make a bounding box
    native <- st_intersects(st_as_sfc(st_bbox(st_buffer(spPt, 100))),
                            regionsExtent, sparse = FALSE)
  }
  if(native == FALSE) {
    write.csv(data.frame(Species = sp),
              file.path(kassDir, "species_rangemaps", "Non-native", paste0(sp, ".csv")),
              quote = FALSE, row.names = FALSE)
  } else {
    
    ### if at least 3 points then calculate alpha hull
    if(nrow(spPt) >= 3) {
      ### MCP around points and convert to sf object using functions modified from functions developed
      ### by Cecina Babich Morrow (https://babichmorrowc.github.io/post/2019-03-18-alpha-hull/
      spMCP <- NULL
      catchError <- try({
        spMCP <- ahull(st_coordinates(spPt), alpha = 15) %>%
          ahull2lines() %>%
          st_set_crs(4326)
      }, silent = TRUE)
      
      ### for some species it is not possible to construct the alpha hull. Construct alpha shape instead
      ### method to convert ashape to sf from https://stackoverflow.com/questions/72449104/converting-alpha-hulls-to-spatial-polygon
      if(is.null(spMCP)) {
        spMCP <- ashape(st_coordinates(spPt), alpha = 30)
        apices <- data.frame(spMCP$edges)[, c( 'x1', 'y1', 'x2', 'y2')]
        edges <- st_linestring(matrix(as.numeric(apices[1, ]), ncol = 2, byrow = TRUE))
        for(j in 2:nrow(apices)){
          edges <- c(edges, st_linestring(matrix(as.numeric(apices[j, ]), ncol = 2, byrow = TRUE)))
        }
        spMCP <- st_sf(geom = st_sfc(edges), crs = 4326) %>%
          st_polygonize() %>%
          st_collection_extract()
      }
      
      ### reproject and buffer MCP by 30 km
      spMCP <- spMCP %>%
        st_transform(st_crs(regions)) %>%
        st_make_valid() %>%
        st_buffer(30000)
    }
    
    ### if < 3 points then use a 30km buffer
    if(nrow(spPt) < 3) {
      spMCP <- spPt %>%
        st_transform(st_crs(regions)) %>%
        st_buffer(30000) %>%
        st_union() %>%
        st_sf() %>%
        mutate(Species = sp)
    }
    
    ### check that the species overlaps with tropical asia
    native <- st_intersects(spMCP, regions, sparse = FALSE)
    if(all(native == FALSE)) {
      write.csv(data.frame(Species = sp),
                file.path(kassDir, "species_rangemaps", "Non-native", paste0(sp, ".csv")),
                quote = FALSE, row.names = FALSE)
    } else {
      ### clip to gadm land
      spMCP <- spMCP %>%
        st_intersection(gadm) %>%
        st_make_valid() %>%
        mutate(Species = sp) %>%
        select(-Land)
      
      ### clip out lakes
      lakesSp <- st_crop(lakes, spMCP) %>%
        st_union() %>%
        st_make_valid()
      if(length(lakesSp) > 0) {
        spMCP <- spMCP %>%
          st_difference(lakesSp) %>%
          st_make_valid()
      }
      
      ### save
      st_write(spMCP, file.path(kassDir, "species_rangemaps", paste0(sp, ".gpkg")),
               quiet = TRUE, append = FALSE)
      
      ### optional plotting to check matches up with original raster
      # library(ggplot2)
      # library(terra)
      # r <- rast(paste0("/mnt/Work/spatial_data/biodiversity/rangemaps/ants/Kass/sdm_species/", gsub("_", ".", sp), ".tif"))
      # names(r) <- "Prob"
      # ggplot() +
      #   geom_tile(data = as.data.frame(r, xy = TRUE), aes(x = x, y = y, fill = Prob)) +
      #   geom_sf(data =  st_transform(st_crop(regions, spMCP), st_crs(r)), fill = NA) +
      #   geom_sf(data =  st_transform(st_crop(lakes, spMCP), st_crs(r)), col = "orange", fill = NA) +
      #   geom_sf(data = st_transform(spMCP, st_crs(r)), col = "red", fill = NA) +
      #   geom_sf(data = st_transform(spPt, st_crs(r)))
      
      # ### append to global results
      # if(is.null(all)) {
      #   all <- spMCP
      # } else {
      #   all <- bind_rows(all, spMCP)
      # }
    }
  }
}

#==================================================================================================#
#------------------------------------- Merge in to single gpkg ------------------------------------#
#==================================================================================================#

### list of all species with alpha hulls - 3619 species
spList <- list.files(file.path(kassDir, "species_rangemaps"),
                     full.names = TRUE, pattern = ".gpkg", recursive = FALSE)

### Append all together
all <- do.call(bind_rows, sapply(spList, function(x) { st_read(x) }))
all <- NULL
pb <- txtProgressBar(1, length(spList), style = 3)
for(i in 1:length(spList)) {
  setTxtProgressBar(pb, i)
  spMap <- st_read(spList[i], quiet = TRUE)
  all <- bind_rows(all, spMap)
}

### Save
st_write(all, file.path(kassDir, "Kass_alpha_hulls.gpkg"))


