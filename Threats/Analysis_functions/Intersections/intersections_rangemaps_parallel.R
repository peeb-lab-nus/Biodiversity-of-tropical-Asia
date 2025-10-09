####################################################################################################
### Parallel version of function for generating intersections with countries, regions or grids
### Charlie Marsh
### charliem2003@github
### 05/2024
###
### returns data frame with the proportion of range area occurring within each region/grid cell and
### also saves species intersections as individual shapefiles
###
### NOTE: the function used for parallelisation (mclapply) only works on linux and mac. For windows
###       machines it will only run the job on a single thread
###
####################################################################################################

intersections_rangemaps_parallel <- function(
    rangemaps,                     # species rangemaps
    spNameCol        = "sci_name", # column name where species names are kept
    regions,                       # shapefile with regions to summarise by (grid or regions)
    regionNameCol    = "ISO",      # column name where unique region names are kept
    gadm,                          # global gadm for initial masking out where coastlines don't overlap
    marine           = FALSE,      # if TRUE masks out land areas. If FALSE masks out marine areas
    crop_to_region   = TRUE,       # exclude rangemaps not overlapping with region?
    areaThresh       = 10000000,   # area threshold in km2 for when to split the shapefile in to multiple parts
    spList,                        # species list to subset rangemaps with
    spDir,                         # directory to save individual species intersections
    overwriteSpFiles = FALSE,      # overwrite existing species intersections, otherwise skips
    threads          = 10          # for parallelisation (for windows only works where threads = 1)
) {
  ### Libraries
  require(sf)
  require(dplyr)
  require(tidyr)
  require(units)
  require(parallel)
  require(pbmcapply)
  
  ##################################################################################################
  ### check projections - st_overlaps assumes maps are projected (not epsg 4326)
  crsRangemap <- st_crs(rangemaps)
  crsRegions  <- st_crs(regions)
  if(crsRangemap != crsRegions) {
    stop("Range maps and regions map are in different projections. Reproject")
  }
  if(is.na(crsRangemap)) {
    stop("Range maps are not assigned a projection.")
  }
  if(is.na(crsRegions)) {
    stop("Regions map has not been assigned a projection.")
  }
  if(!is.na(crsRangemap$epsg) & crsRangemap$epsg == 4326) {
    warning("Range maps are not in a planar projection (epsg = 4326). ",
            "Intersections will probably be calculated incorrectly. ",
            "Reproject or set sf_use_s2(FALSE) (this is much much slower!).")
  }
  if(!is.na(crsRegions$epsg) & crsRegions$epsg == 4326) {
    warning("Regions map is not in a planar projection (epsg = 4326). ",
            "Intersections will probably be calculated incorrectly. ",
            "Reproject or set sf_use_s2(FALSE) (this is much much slower!).")
  }
  
  ##################################################################################################
  ### calculate intersections for each species
  
  ### if cropRegion then use simple bbox to check if any part of rangemap overlaps region
  if(crop_to_region == TRUE) {
    regionsExtent <- st_as_sfc(st_bbox(regions))
  }
  
  ### loop through species list
  pbmclapply(1:length(spList), function(i) {
    sp <- spList[i]
    spFile <- file.path(spDir, paste0(gsub(" ", "_", sp), ".csv"))
    
    ### file to show which species are running
    runningFile <- file.path(spDir, "Running", paste0(gsub(" ", "_", sp), ".txt"))
    write.table(data.frame(Species = sp), runningFile, quote = FALSE, row.names = FALSE)
    
    catchError <- try({
      if(overwriteSpFiles == TRUE | file.exists(spFile) == FALSE) {
        ### extract species map
        spMap <- filter(rangemaps, get(spNameCol) == sp)
        
        ### to speed things up, check whether extent of species rangemaps overlaps with asia
        extOverlap <- FALSE
        if(crop_to_region == TRUE) {
          extOverlap <- st_intersects(regionsExtent,
                                      st_make_valid(st_as_sfc(st_bbox(spMap))),  sparse = FALSE)
        }
        
        if(crop_to_region == TRUE & extOverlap == FALSE) {
          write.csv(data.frame(Species = sp,
                               Range_area = round(drop_units(sum(set_units(st_area(spMap), km^2))), 4)),
                    file.path(spDir, "Outside_extent", paste0(gsub(" ", "_", sp), ".csv")),
                    quote = FALSE, row.names = FALSE)
        }
        
        ### if there is overlap calculate area overlap within each region
        if(crop_to_region == FALSE | extOverlap == TRUE) {
          
          ### Some BirdLife International/IUCN range maps have invalid MULTISURFACE geometries
          # repair if necessary (code adapted from https://github.com/r-spatial/sf/issues/2289)
          badGeometry <- try(st_cast(spMap, "MULTIPOLYGON"), silent = TRUE)
          if(any(class(badGeometry) == "try-error")) {
            
            ### use gdal instead
            if(file.exists(file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), ".gpkg")))) { 
              file.remove(file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), ".gpkg")))
            }
            if(file.exists(file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), "_fix.gpkg")))) {
              file.remove(file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), "_fix.gpkg")))
            }
            st_write(spMap, file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), ".gpkg")))
            
            ### gdal changed from -append to -upsert from v3.6
            vGDAL <- system("gdal-config --version", intern = TRUE)
            vGDAL <- as.numeric(paste(strsplit(vGDAL, "[.]")[[1]][1:2], collapse = "."))
            system(paste("ogr2ogr",
                         file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), "_fix.gpkg")),
                         file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), ".gpkg")),
                         ifelse(vGDAL < 3.6, "-append", "-upsert"),
                         "-skipfailures -explodecollections -nlt CONVERT_TO_LINEAR"))
            spMap <- st_read(file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), "_fix.gpkg"))) %>%
              st_make_valid()
            file.remove(file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), "_fix.gpkg")))
            file.remove(file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), ".gpkg")))
          }
          
          ### Merge multiple polygons if necessary
          spMap <- spMap  %>%
            st_cast("MULTIPOLYGON") %>%
            st_make_valid() %>%
            select(all_of(spNameCol)) %>%
            st_union() %>%
            st_sf() %>%
            mutate(Species =  gsub(" ", "_", sp))
          
          ### check that range map has valid spherical geometry (some older iucn maps aren't)
          if(any(st_is_valid(spMap) == FALSE)) {
            spMap <- st_make_valid(spMap)
          }
          
          ### First, mask out global range based on GADM 
          ###  - many IUCN maps use very inaccurate coastlines which messes up the overlaps with our
          ###    bioregions, so we overlap the global range map with GADM and remove any non-overlaps
          
          ### If range size is very large then split into smaller chunks first
          if(shpBBoxArea(spMap) > areaThresh) {
            spMapQrtrs <- splitShp(speciesShp = spMap, areaThresh = areaThresh, verbose = FALSE)
            
            ### if you want to plot the splits
            # par(mfrow = c(floor(sqrt(length(speciesShp) + 1)), ceiling(sqrt(length(speciesShp) + 1))))
            # plot(spMap[1], border = T, reset = F, main = "Original rangemap")
            # lapply(speciesShp, function(x) plot(x[1], border = F, reset = F, main = NULL))
            # par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))
            
            ### mask each split to the gadm
            if(marine == FALSE) {
              spMapQrtrs <- lapply(spMapQrtrs, function(x) { st_make_valid(st_intersection(x, gadm)) })
            }
            if(marine == TRUE) {
              spMapQrtrs <- lapply(spMapQrtrs, function(x) { st_make_valid(st_difference(x, gadm)) })
            }
            
            ### In Phalaropus fulicarius there is one quarter that ends up with no data.
            ### Remove this special case
            spMapQrtrs <- spMapQrtrs[lapply(spMapQrtrs, nrow) > 0]
            
            ### merge back together
            spMap <- bind_rows(spMapQrtrs) %>%
              st_union() %>%
              st_make_valid()
            
            ### in some cases this makes a geometry collection rather than a single multipolygon
            ### because we are left with points (rather than polygons) that need to be removed
            if(st_geometry_type(spMap) == "GEOMETRYCOLLECTION") {
              geoms <- lapply(spMap[1], `[`)[[1]]
              polys <- lapply(1:length(geoms), function(x) { st_geometry_type(geoms[[x]]) })
              polys <- which(unlist(polys) %in% c("POLYGON", "MULTIPOLYGON"))
              geoms <- geoms[polys]
              mp <- lapply(geoms, function(x) st_polygon(x = x))
              sfc_mp <- st_sfc(mp) %>%
                st_cast(to = "MULTIPOLYGON") %>%
                st_union() %>%
                st_make_valid()
              
              ### and make final multipolygon
              spMap <- st_sf(Species = sp,
                             geom =  sfc_mp,
                             crs = st_crs(spMap))
            } else {
              spMap <- st_as_sf(spMap, Species = sp)
            }
          }
          
          ### Otherwise mask the whole map with the gadm
          if(shpBBoxArea(spMap) <= areaThresh) {
            if(marine == FALSE) {
              spMap <- spMap %>%
                st_intersection(gadm) %>%
                st_make_valid()
            }
            if(marine == TRUE) {
              spMap <- spMap %>%
                st_difference(gadm) %>%
                st_make_valid()
            }
          }
          
          ### Total range area
          spMap <- spMap %>%
            mutate(Range_area = round(drop_units(set_units(st_area(.), km^2)), 4))
          
          ### calculate area of overlap between range map and each region. Some cases we may need to lower precision
          canIntersect <- try({
            overlap <- spMap %>%
              st_set_precision(0) %>%
              st_intersection(regions)
          }, silent = TRUE)
          if(inherits(canIntersect, "try-error")) {
            canIntersect <- try({
              overlap <- spMap %>%
                st_set_precision(1e5) %>%
                st_make_valid() %>%
                st_intersection(regions)
            }, silent = TRUE)
          }
          if(inherits(canIntersect, "try-error")) {
            canIntersect <- overlap <- spMap %>%
                st_set_precision(1e4) %>%
                st_make_valid() %>%
                st_intersection(regions)
          }
          overlap <- overlap %>%
            st_make_valid() %>%
            rename(Region = all_of(regionNameCol)) %>%
            mutate(Intersect_area = round(drop_units(set_units(st_area(.), km^2)), 4))
          
          ### save csv depending on whether there is overlap with bioregions
          if(nrow(overlap) == 0) {
            write.csv(data.frame(Species = sp,
                                 Range_area = round(drop_units(sum(set_units(st_area(spMap), km^2))), 4)),
                      file.path(spDir, "Non-native", paste0(gsub(" ", "_", sp), ".csv")),
                      quote = FALSE, row.names = FALSE)
          }
          
          if(nrow(overlap) > 0) {
            ### tidying up and save
            overlap <- st_drop_geometry(data.frame(overlap)) %>%
              select(Species, Region, Range_area, Intersect_area) %>%
              relocate(Species)
            
            write.csv(overlap, spFile, quote = FALSE, row.names = FALSE)
          }
        }
      }
      file.remove(runningFile)
      
      return(NULL)
    }, silent = TRUE)
    
    ### save error to file
    if(inherits(catchError, "try-error")) {
      write.table(catchError, file.path(spDir, "Errors", paste0(gsub(" ", "_", sp), ".txt")),
                  quote = FALSE, row.names = FALSE)
    }
    
    ### remove running file
    file.remove(runningFile)
  }, mc.cores = threads)
  
  return(NULL)
}
