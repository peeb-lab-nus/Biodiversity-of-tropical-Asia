####################################################################################################
### Run cleaning steps for GBIF data.
###
### The filters to apply, which order to run them in, and any necessary parameters are read in from
### the filterSteps spreadsheet.
### 
### Returns:
###   Original GBIF data frame with columns for whether records were flagged ('Flag'), and if so
###    at which stage ('Flag_step')
###   Optionally saves full GBIF data frame with new flag columns (pathFlaggedData)
###   Optionally saves filtered GBIF data frame with flagged records removed (pathFilteredData)
### 
### Charlie Marsh
### charliem2003@github
### 10/2024
####################################################################################################

runCleaningSteps <- function(gbif,                    # prepared gbif data frame, with cols Flag and Flag_step added
                             filterSteps,             # data frame with details of flags to apply
                             id_col,                  # Column name containing unique record ID
                             lon_colname,             # Column name containing longitude coords
                             lat_colname,             # Column name containing latitude coords
                             pathFlaggedData  = NULL, # path and file name of full gbif data with column indicated if flag applied. If NULL saves nothing
                             pathFilteredData = NULL, # path and file name of gbif data with flagged data filtered out. If NULL saves nothing
                             plot             = FALSE # optional plotting of points and flags
) {
  ### Get steps to be taken from filterSteps csv
  steps <- sort(filterSteps$Order, decreasing = FALSE)
  
  for(i in steps) {
    ### Retrieve function for flagging step
    stepFun <- get(filterSteps$Function[filterSteps$Order == i])
    
    ### Run
    gbif <- do.call(stepFun,
                    args = list(gbif        = gbif,
                                filterSteps = filterSteps,
                                id_col      = colID,
                                lon_colname = colLon,
                                lat_colname = colLat))
  }
  
  ### Write full data frame with column indicating flagged data to file
  if(!is.null(pathFlaggedData)) {
    print(paste0("Saving flagged gbif data to ", pathFlaggedData))
    write_csv(gbif, pathFlaggedData)
  }
  
  ### Write filtered data frame to file
  if(!is.null(pathFilteredData)) {
    print(paste0("Saving flagged gbif data to ", pathFilteredData))
    gbifFiltered <- gbif %>%
      filter(is.na(Flag))
    write_csv(gbifFiltered, pathFilteredData)
  }
  
  ### Optional plotting
  if(plot == TRUE) {
    print("Plotting results ...")
    require(ggplot2)
    require(sf)
    require(rnaturalearth)
    
    ### Convert data frame to points
    gbifPts <- gbif %>%
      select(all_of(c(colID, colLon, colLat, "Flag", "Flag_step"))) %>%
      filter(!is.na(get(colLon)) & !is.na(get(colLat))) %>%
      mutate(Flag_step = formatC(Flag_step, flag = "0", width = 2)) %>%
      mutate(Flag = apply(., 1, function(x) paste(x["Flag_step"], x["Flag"], sep = "_"))) %>%
      mutate(Flag = case_when(Flag == "NA_NA" ~ "Passed",
                              .default = Flag)) %>%
      st_as_sf(coords = c(colLon, colLat), crs = 4326)
    
    ### Background land layer to plot
    # land <- rnaturalearth::ne_countries(scale = 50, type = "countries") %>%
    #   st_transform(st_crs(gbifPts)) %>%
    #   st_crop(st_bbox(gbifPts)) %>%
    #   st_make_valid()
    
    land <- rnaturalearth::ne_coastline(scale = 50) %>%
      st_transform(st_crs(gbifPts)) %>%
      # st_crop(st_bbox(gbifPts)) %>%
      st_make_valid()
    gbifPts <- st_transform(gbifPts, st_crs(land))
    
    p1 <- ggplot() +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = NA, colour = "black")) +
      coord_sf(xlim = st_bbox(gbifPts)[c("xmin", "xmax")],
               ylim = st_bbox(gbifPts)[c("ymin", "ymax")],
               expand = FALSE) +
      geom_sf(data = land) +# , fill = "grey70") +
      geom_sf(data = filter(gbifPts, Flag != "Passed"), aes(col = Flag)) +
      geom_sf(data = filter(gbifPts, Flag == "Passed"), aes(shape = Flag), col = "black")
    print(p1)
  }
  
  ### Return full flagged data frame
  print(paste0(prettyNum(sum(!is.na(gbif$Flag)), big.mark = ","), " records flagged out of ",
               prettyNum(nrow(gbif), big.mark = ","), " records. ",
               prettyNum(sum( is.na(gbif$Flag)), big.mark = ","), " records remaining"))
  return(gbif)
}
