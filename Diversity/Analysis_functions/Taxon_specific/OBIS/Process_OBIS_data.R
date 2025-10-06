####################################################################################################
### Reads in the OBIS data dump in chunks, does some initial basic filtering, writes them to file
### and finally appends them all together
###
### Charlie Marsh
### charliem2003@github
### 10/2024
####################################################################################################

### Load libraries
library(dplyr)
library(readr)

### Directory where the data is, and also where we will save results - ADJUST TO YOUR STRUCTURE
baseDir <- file.path("Directory", "to", "write", "results") # Dir where all results will be stored
obisDir <- file.path("OBIS", "data dump", "directory")      # Dir where the OBIS data dump is located

### Directory to save subset files to
subsetDir <- file.path(baseDir, "subsets")
if(!file.exists(subsetDir)) { dir.create(subsetDir) }

### Determine number of rows of data frame - 108,504,418
### (linux only, windows needs to be something like pipe("C:\\Rtools\\bin\\wc -l myfile.csv")))
csvRows <- system(paste0("wc -l ", file.path(baseDir, "obis_20230208.csv")), intern = TRUE)
csvRows <- as.numeric(strsplit(csvRows, " ")[[1]][1]) - 1 # 1st row is column headers

### Get column names
csvCols <- read_csv(file.path(baseDir, "obis_20230208.csv"), n_max = 1)
csvCols <- colnames(csvCols)

### Which columns are likely to be unnecessary?
### Determine from first million rows which columns mostly contain NAs
# obisSub <- read_csv("/mnt/Work/NUS/obis_20230208.csv", n_max = 1000000)
# csvColsSub <- apply(obisSub, 2, function(x) sum(is.na(x))) / 1000000
# csvColsSub <- names(csvColsSub[csvColsSub < 1])

csvColsSub <- c("id", "decimallongitude", "decimallatitude", "scientificname", "originalscientificname",
                "depth", "coordinateuncertaintyinmeters", "flags", "marine", "brackish", "freshwater",
                "terrestrial", "taxonrank", "redlist_category", "superdomain", "kingdom", "phylum",
                "class", "order", "family", "genus", "species", "basisofrecord","individualcount",
                "occurrencestatus", "eventdate", "year", "month", "day", "typestatus", "scientificnameauthorship")

### We'll go in 1,000,000 row chunks
rowSeq <- seq(1, csvRows, by = 1000000)
pb <- txtProgressBar(1, length(rowSeq), style = 3)
for(i in 1:length(rowSeq)) {
  setTxtProgressBar(pb, i)
  obisSub <- read_csv(file.path(baseDir, "obis_20230208.csv"),
                      skip = rowSeq[i], n_max = 1000000,
                      col_names = csvCols, col_select = csvColsSub,
                      progress = FALSE, show_col_types = FALSE)
  
  ### Apply initial basic filtering
  obisSub <- obisSub %>%
    filter(taxonrank == "Species") %>%
    filter(marine == 1 | brackish == 1) %>%
    filter(!is.na(decimallatitude) & !is.na(decimallongitude))
  
  ### write to csv
  if(nrow(obisSub) > 0) {
    write_csv(obisSub, file.path(subsetDir, paste0("obis_20230208_subset", i, ".csv")),
              progress = FALSE, append = FALSE)
  }
}

### Bind everything together
all <- tibble()
pb <- txtProgressBar(1, length(rowSeq), style = 3)
for(i in 1:length(rowSeq)) {
  setTxtProgressBar(pb, i)
  
  ### Get column types from first file
  if(i == 1) {
    colTypes <- attributes(obisSub)$spec
  }
  
  obisSub <- read_csv(file.path(subsetDir, paste0("obis_20230208_subset", i, ".csv")),
                      progress = FALSE, show_col_types = FALSE, col_types = colTypes)
  all <- bind_rows(all, obisSub)
}

### There shouldn't be any duplicated entries
any(duplicated(all$id))

### And save
write_csv(all, file.path(baseDir, "obis_20230208_filtered.csv"), progress = TRUE, append = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
gc()
