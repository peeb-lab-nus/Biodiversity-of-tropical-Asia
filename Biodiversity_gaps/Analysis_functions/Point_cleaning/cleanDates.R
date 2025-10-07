####################################################################################################
### cleanDates.R
### Charlie Marsh
### charliem2003@github
### 07/2025
###
### Function for cleaning and extracting dates in GBIF and OBIS. Dates are messy (e.g. 2002/2011
### and 2018-07-31T06:30:00Z etc)
###
### 1. Now functional when dates are non-numeric or ambiguous (e.g 09-Sep-09)
### 2. Defaults to eventDate over year/month/day
### 3. If there is a month and/or day, but not year/month, resets them to NA
### 4. Cleans out impossible numbers
### 5. If a range of dates is given (e.g 2007-01/2007-02), keep earliest date
### 
### Returns:
###   The original GBIF data frame with new columns dateClean, yearClean, monthClean, dayClean but
###   keeps original data as well
###
####################################################################################################

cleanDates <- function(gbif) {
  require(lubridate)
  require(dplyr)
  require(stringi)
  
  # ### First standardise names
  # gbif <- gbif %>%
  #   mutate(eventName = select(matches("eventdate",   ignore.case = TRUE))) %>%
  #   mutate(dayName   = if_all(matches("day",   ignore.case = TRUE))) %>%
  #   mutate(dayName   = as.numeric(dayName)) %>%
  #   mutate(monthName = if_all(matches("month", ignore.case = TRUE))) %>%
  #   mutate(monthName = as.numeric(monthName)) %>%
  #   mutate(yearName  = if_all(matches("year",  ignore.case = TRUE))) %>%
  #   mutate(yearName  = as.numeric(yearName))
  # ### First standardise names
  gbif$eventName <- unlist(as.vector(gbif[, which(tolower(colnames(gbif)) == "eventdate")]))
  gbif$dayName   <- unlist(as.vector(gbif[, which(tolower(colnames(gbif)) == "day")]))
  gbif$dayName   <- as.numeric(gbif$dayName)
  gbif$monthName <- unlist(as.vector(gbif[, which(tolower(colnames(gbif)) == "month")]))
  gbif$monthName <- as.numeric(gbif$monthName)
  gbif$yearName  <- unlist(as.vector(gbif[, which(tolower(colnames(gbif)) == "year")]))
  gbif$yearName  <- as.numeric(gbif$yearName)
  
  ### Then extract eventDate to account for variations in entries (such as "01/01/2001-02/01/2001" vs "01-01-2001/02-01-2001")
  gbif$nSlash  <- stri_count_fixed(gbif$eventName, "/")
  gbif$nHyphen <- stri_count_fixed(gbif$eventName, "-")
  gbif <- gbif %>%
    mutate(dateClean = case_when(nSlash == 0 & nHyphen >= 0 ~ eventName,
                                 nSlash == 1 & nHyphen == 0 ~ gsub("/", "-", eventName),
                                 nSlash == 1 & nHyphen == 1 ~ max(stri_split_regex(eventName, "[/-]")[[1]]),
                                 nSlash == 1 & nHyphen  > 1 ~ eventName,
                                 nSlash  > 1 & nHyphen == 0 ~ gsub("/", "-", eventName),
                                 nSlash  > 1 & nHyphen == 1 ~ gsub("/", "-", stri_split_fixed(eventName, "-")[[1]][1]),
                                 # nSlash  > 1 & nHyphen  > 1 ~ ,
                                 .default = NA))
  
  ### Where we have data range ("/") or time ("T" and " ") take first split
  gbif$dateClean <- stri_split_fixed(gbif$dateClean, "/", simplify = TRUE)[, 1]
  gbif$dateClean <- stri_split_fixed(gbif$dateClean, "T", simplify = TRUE)[, 1]
  gbif$dateClean <- stri_split_fixed(gbif$dateClean, " ", simplify = TRUE)[, 1]
  
  ### Deals with cases of eg "2016-2017"
  gbif$dateClean <- sapply(gbif$dateClean, function(x) {
    splits <- as.numeric(stri_split_fixed(x, "-")[[1]])
    return(ifelse(length(splits) > 1 & min(splits, na.rm = TRUE) > 1000 & max(splits, na.rm = TRUE) > 1000,
                  min(splits),
                  x))
  })
  
  ### Is date in ymd, ym or just y - this vector saves that info because 'year' doesn't work with only year
  gbif$type <- stri_count_fixed(gbif$dateClean, "-")
  
  ### Now assign final values based on dateClean, or if that is missing, based on the specific columns
  gbif <- gbif %>%
    mutate(yearClean = case_when(type == 0 & !is.na(dateClean) ~ as.numeric(dateClean),
                                 type == 1 ~ as.numeric(year(parse_date_time(dateClean, c("my",  "ym")))),
                                 type == 2 ~ as.numeric(year(parse_date_time(dateClean, c("dmy", "ymd")))),
                                 .default = as.numeric(yearName))) %>%
    
    mutate(monthClean = case_when(type == 1 ~ as.numeric(month(parse_date_time(dateClean, c("my", "ym")))),
                                  type == 2 ~ as.numeric(month(parse_date_time(dateClean, c("dmy", "ymd")))),
                                  .default =  as.numeric(monthName))) %>%
    
    mutate(dayClean = case_when(type == 2 ~ as.numeric(day(parse_date_time(dateClean, c("dmy", "ymd")))),
                                .default = as.numeric(dayName)))
  
  ### If there is no valid year, but there is month and/or day, set everything to NA
  gbif <- gbif %>%
    mutate(monthClean = case_when(is.na(yearClean) ~ NA,
                                  .default = monthClean)) %>%
    mutate(dayClean   = case_when(is.na(monthClean) ~ NA,
                                  .default = dayClean))
  
  ### OBIS has some impossible dates, in which case set to NAs
  gbif <- gbif %>%
    mutate(yearClean = case_when(yearClean < 1000 ~ NA,
                                 yearClean > 2024 ~ NA,
                                 .default = yearClean)) %>%
    mutate(monthClean = case_when(monthClean < 1  ~ NA,
                                  monthClean > 12 ~ NA,
                                  .default = monthClean)) %>%
    mutate(dayClean = case_when(dayClean < 1  ~ NA,
                                dayClean > 31 ~ NA,
                                .default = dayClean))
  
  ### Put together final clean data
  gbif <- gbif %>%
    mutate(dateClean = paste(yearClean,
                             formatC(as.numeric(monthClean), width = 2, format = "d", flag = "0"),
                             formatC(as.numeric(dayClean),   width = 2, format = "d", flag = "0"),
                             sep = "-")) %>%
    mutate(dateClean = gsub("-NA", "", dateClean)) %>%
    mutate(dateClean = case_when(dateClean == "NA" ~ NA, .default = dateClean)) %>%
    
    ### Remove unncessary columns
    select(-eventName, -dayName, -monthName, -yearName, -nSlash, -nHyphen, -type)
  
  return(gbif) 
}


# (gbifClean <- tibble(eventdate = c("2004-04-24/2004-04-28", "2004/04/24-2004/04/28", "2016-03-25", "03-2016", "2016-03", "2016",   NA,   NA, "15-Sep-84",    NA, NA, "2016-03-25 00:00", "2016/2017", "2016"),
#                      year      = c(                   2004,                    2004,           NA,        NA,        NA,   2016, 2016, 2016,          NA,   900, NA,                 NA,          NA,     NA),
#                      Month     = c(                     04,                      04,           NA,        NA,        03,     03,   03,   NA,          NA, "Sep", 10,                 NA,          NA,     NA),
#                      day       = c(                     24,                      24,           NA,        NA,        NA,     24,   24,   NA,          NA,    35, 15,                 NA,          NA,     NA)))
# cleanDates(gbif = gbifClean)

