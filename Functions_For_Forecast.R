################################################################################
#                                                                              #
# Functions for generating seasonal forecasts for the red-billed quelea (Quelea#
# quelea)                                                                      # 
#                                                                              #
# Code for: Seasonal forecasting of mobile species distributions               #
# for dynamic management under extreme weather events.                         #
#                                                                              #
# Rachel Dobson, Stephen G. Willis, Stewart Jennings, Robert A. Cheke,         #
# Andrew J. Challinor and Martin Dallimer                                      #
#                                                                              #
################################################################################
#
# This script does not need to be opened or edited. It contains custom functions
# for generating seasonal forecasts of red-billed quelea distribution. This file
# will be used to load functions into the main script
# "Seasonal_Forecast_Quelea.R"

#' get_evi_characterstics Extracts EVI characterstics for vegetation growth
#' stage classification using 16-day MODIS EVI values from past 52-weeks.
#' @param directory path to 16-day EVI rasters. 
#' @param land_cover raster of land cover cells from MODIS Land Cover Yearly. 
#' @param spatial_ext spatial polygon, the spatial extent to crop EVI data to. 

get_evi_characterstics <- function(directory,
                                   land_cover,
                                   spatial_ext) {
  
  library(dplyr)
  
  # List and stack the EVI rasters 
  evi_files <- list.files(directory, full.names = T)
  
  stack <-  evi_files %>% # list files to read in
    lapply(terra::rast)
  
  stack <- terra::rast(stack)
  
  
  # Extract EVI characteristics for stage classification
  
  peakinlastyear <- max(stack, na.rm = T)
  
  troughinlastyear <- min(stack, na.rm = T)
  
  nolay <- terra::nlyr(stack)
  
  timesincepeak <- (nolay - which.max(stack))
  
  currentEVI <- stack[[nolay]]
  
  changeinpastfourweek <- currentEVI - stack[[nolay - 2]]
  
  changeinpastsixweek <- currentEVI - stack[[nolay - 3]]
  
  changeinpasteightweek <- currentEVI - stack[[nolay - 5]]
  
  changeinsixteenweek <- currentEVI - stack[[nolay - 8]]
  
  changeintwentysixweek <- currentEVI - stack[[nolay - 13]]
  
  changeinthirtytwoweek <- currentEVI - stack[[nolay - 16]]
  
  changeinfourtyweek <- currentEVI - stack[[nolay - 20]]
  
  percentofamplitude <- currentEVI / peakinlastyear - troughinlastyear
  
  # Resample land cover data to match the EVI data
  land_cover2<- resample(land_cover,changeinsixteenweek, method = "near")

  # Combine all rasters into one stack
  combined_rast <- terra::rast(list(land_cover2,
                                    peakinlastyear,
                                    troughinlastyear,
                                    timesincepeak,
                                    currentEVI,
                                    changeinpastfourweek,
                                    changeinpastsixweek,
                                    changeinpasteightweek,
                                    changeinsixteenweek,
                                    changeintwentysixweek,
                                    changeinthirtytwoweek,
                                    changeinfourtyweek,    
                                    percentofamplitude))

  # Rename the layers to match that of classification model 
  names(combined_rast) <- c("LC",
                            "peakinlastyear",
                            "troughinlastyear",
                            "timesincepeak",
                            "currentEVI",
                            "changeinpastfourweek",
                            "changeinpastsixweek",
                            "changeinpasteightweek",
                            "changeinsixteenweek",
                            "changeintwentysixweek",
                            "changeinthirtytwoweek",
                            "changeinfourtyweek",
                            "percentofamplitude")
  
  # Crop and mask rasters to specified spatial extent
  combined_rast <- terra::crop(combined_rast, spatial_ext)
  combined_rast <- terra::mask(combined_rast, spatial_ext)
  
  # Transform into data frame for Random Forest model projections
  combined_rast_df <- terra::as.data.frame(combined_rast, xy = TRUE)
  
  return(combined_rast_df)
  
}

#' project_seed_availability Projects seed availability for each forecast
#' interval using EVI characteristics and average vegetation growth stage
#' lengths. 
#' @param EVI_data_frame a data frame, the output of `get_evi_characterstics()`.
#' @param mean_lengths raster stack, the average lengths of each vegetation growth stage. 
#' @param type one of; `cereal` or `grass`, the land cover cell type to project
#'   seed availability for.
#' @param model Random Forest model for classifying vegetation growth stages based upon EVI characteristics. 
#' @param save_dir path to save seed availability projections to. 
#' @param forecast_intitiation a character, the date that the seasonal forecast was initiated.  
#' @param forecast_intervals a character vector, the date of each seasonal forecast interval. 

project_seed_availability <- function(EVI_data_frame ,
                                      mean_lengths,
                                      type,
                                      model,
                                      save_dir,
                                      forecast_intitiation,
                                      forecast_intervals) {
  
  
  # Set value used for vegetation type in MODIS Land Cover Yearly dataset
  if (type == "cereal") {
    LC_TYPE <- 7
  }
  
  if (type == "grass") {
    LC_TYPE <- 6
  }
  
  # Filter EVI to just cells of specified land cover type. 
  EVI_data_frame<-EVI_data_frame[EVI_data_frame$LC==LC_TYPE,]
  
  # Get month that EVI characteristics were extracted for. 
  EVI_data_frame$month<-rep(lubridate::month(forecast_intitiation-14),nrow(EVI_data_frame))
  
  # Project Random Forest classification model onto EVI characterstics
  EVI_data_frame$preds<-predict(model,EVI_data_frame)
  
  # If NA value (no EVI characteristics available), set cell as vegetation dormancy
  EVI_data_frame$preds[is.na(EVI_data_frame$preds)]<-7
  
  vegetation_stages <- data.frame(unique(cbind(c("_Greenup",
                                                 "MidGreenup",
                                                 "Senescence",
                                                 "MidGreendown",
                                                 "dormancy"), 
                                               c(1, 2, 5, 6, 7))))
  
  stage_track<-NULL
  
  # Iterates through each vegetation stage in chronological order,
  # adds the average length to initial stage. 
  
  for (i in 1:nrow(vegetation_stages)){
    
    name <- vegetation_stages[i, 1]
    
    number <- vegetation_stages[i, 2]
   
    croped <-  dplyr::filter(EVI_data_frame, preds == number)
    
    if(i ==1){
      
      phenology<-terra::extract(mean_lengths,y=as.matrix(croped[,c("x","y")])) 
      phenology[,1] <- phenology[,1]/2
      phenology[,2]<-phenology[,1]+phenology[,2]
      phenology[,3]<-phenology[,2]+phenology[,3]
      phenology[,4]<-phenology[,3]+phenology[,4]
      
      stage_track<-rbind(stage_track,cbind(croped[,c("x","y")],phenology))
      
    }
    
    if(i ==2){
      
      phenology<-terra::extract(mean_lengths,y=as.matrix(croped[,c("x","y")])) 
      phenology[,1]<-rep(0,nrow(phenology))
      phenology[,2]<-phenology[,1]+phenology[,2]
      phenology[,2] <- phenology[,2]/2
      phenology[,3]<-phenology[,2]+phenology[,3]
      phenology[,4]<-phenology[,3]+phenology[,4]
      stage_track<-rbind(stage_track,cbind(croped[,c("x","y")],phenology))}
    
    
    if(i ==3){
      
      phenology<-terra::extract(mean_lengths,y=as.matrix(croped[,c("x","y")])) 
      phenology[,1]<-rep(0,nrow(phenology))
      phenology[,2]<-rep(0,nrow(phenology))
      phenology[,3]<-phenology[,2]+phenology[,3]
      phenology[,3] <- phenology[,3]/2
      phenology[,4]<-phenology[,3]+phenology[,4]
      stage_track<-rbind(stage_track,cbind(croped[,c("x","y")],phenology))}
    
    
    
    if(i ==4){
      
      phenology<-terra::extract(mean_lengths,y=as.matrix(croped[,c("x","y")])) 
      phenology[,1]<-rep(0,nrow(phenology))
      phenology[,2]<-rep(0,nrow(phenology))
      phenology[,3]<-rep(0,nrow(phenology))
      phenology[,4]<-phenology[,3]+phenology[,4]
      phenology[,4] <- phenology[,4]/2
      stage_track<-rbind(stage_track,cbind(croped[,c("x","y")],phenology))}
    
    
    if(i ==5){
      
      phenology<-terra::extract(mean_lengths,y=as.matrix(croped[,c("x","y")])) 
      phenology[,1]<-rep(0,nrow(phenology))
      phenology[,2]<-rep(0,nrow(phenology))
      phenology[,3]<-rep(0,nrow(phenology))
      phenology[,4]<-rep(0,nrow(phenology))
      
      stage_track<-rbind(stage_track,cbind(croped[,c("x","y")],phenology))}
    
    
  }
  
  colnames(stage_track)<-c("x","y","midgreenup","peak","midgreen","dorm")
  
  forecast_doys<-as.numeric(forecast_intervals - (forecast_intitiation-14))
  

  # Iterate through each forecast interval, and if vegetation growth stage
  # between peak and dormancy, then seed inferred to be available.
  for(int in 1:length(forecast_intervals)){
    
    date_1 <- forecast_intervals[int]
    d1 <- forecast_doys[int]
    
    seed <- matrix(0, nrow = nrow(stage_track), ncol = 1)
    seed <- dplyr::between(rep(d1, nrow(stage_track)), stage_track[, 4] - 7, stage_track[, 6] + 7)
    seed[stage_track[, 6] == 0] <- 0
    rastero <- terra::rast(as.matrix(cbind(stage_track[, c("x", "y")],
                                           as.numeric(seed))),
                           type = "xyz")
    
    # Aggregate to sum precipitation at coarse resolution (e.g. 12 for 0.05 degree)
    rastero <- terra::aggregate(rastero, agg.factor, fun = "sum", na.rm = TRUE)
    
    # If data are categorical then moving.window.matrix with weights = 1
    moving.window.matrix[1:nrow(moving.window.matrix),
                         1:ncol(moving.window.matrix)] <- 1
    
    # Calculate sum across moving.window.matrix
    # This sums total available seed in surrounding radius. 
    rastero <- terra::focal(rastero,
                            moving.window.matrix,
                            fun = "sum",
                            na.rm = TRUE)
    
    terra::writeRaster(
      rastero,
      file = paste0(save_dir, "/", date_1, "_", type, "_seed.tif"),
      overwrite = T
    )
    
  }
  
  
}




# Function to generate next 7 dates, one month apart
generate_dates <- function(start_date) {
  
  # Convert the input to a Date object
  start_date <- as.Date(start_date)
  
  # Create a vector of 7 subsequent months
  subsequent_dates <- seq(from = start_date, by = "month", length.out = 8)
  
  # Format the dates and return them
  return(as.Date(subsequent_dates)) # Exclude the starting date
}

# Function to generate next 7 dates, one month apart
generate_dates_past <- function(start_date) {
  
  # Convert the input to a Date object
  start_date <- as.Date(start_date)
  year_date <- start_date -years(1)
  # Create a vector of 7 subsequent months
  subsequent_dates <- seq(from = year_date, by = "month", length.out = 12)
  
  # Format the dates and return them
  
  return(as.Date(subsequent_dates)) # Exclude the starting date
}


# Function to get indices for a specific month and year in a raster stack
get_monthly_indices <- function(year, month, start_year = 2002, end_year = 2016) {
  
  duration <- (end_year - start_year)+2
  
  daysin <- as.numeric(lubridate::days_in_month(month))
  
  total_layers <- daysin * duration
  
  yearsin <-(year - start_year) +1
  
  layers <- seq(from = 1, to = total_layers, by = daysin)
  
  end <- layers[yearsin+1]-1
  
  start<- layers[yearsin]
  
  return(c(start))
}



# Function to get indices for a specific month and year in a raster stack
get_monthly_indices_end <- function(year, month, start_year = 2002, end_year = 2016) {
  
  duration <- (end_year - start_year)+2
  
  daysin <- as.numeric(lubridate::days_in_month(month))
  
  total_layers <- daysin * duration
  
  yearsin <-(year - start_year) +1
  
  layers <- seq(from = 1, to = total_layers, by = daysin)
  
  end <- layers[yearsin+1]-1
  
  start<- layers[yearsin]
  
  return(c(end))
}



# Function to get indices for a specific month and year in a raster stack
get_monthly_indices_past <- function(year, month, start_year = 2002, end_year = 2016) {
  
  duration <- (end_year - start_year)+2
  
  daysin <- as.numeric(lubridate::days_in_month(month))
  
  total_layers <- daysin * duration
  
  yearsin <-(year - start_year) +1
  
  layers <- seq(from = 1, to = total_layers, by = daysin)
  
  end <- layers[yearsin+1]-1
  
  start<- layers[yearsin]
  
  return(c(start))
}

# Function to get indices for a specific month and year in a raster stack
get_monthly_indices_end_past <- function(year, month, start_year = 2002, end_year = 2016) {
  
  duration <- (end_year - start_year)+2
  
  daysin <- as.numeric(lubridate::days_in_month(month))
  
  total_layers <- daysin * duration
  
  yearsin <-(year - start_year) +1
  
  layers <- seq(from = 1, to = total_layers, by = daysin)
  
  end <- layers[yearsin+1]-1
  
  start<- layers[yearsin]
  
  return(c(end))
}




get_details <- function(forecast_date, ee, bias_corrected_SEAS5){
  
  next_dates <- generate_dates(forecast_date)
  
  dataframe <- data.frame(date = next_dates,
                          month = lubridate::month(next_dates),
                          year = lubridate::year(next_dates)  )
  dataframe <- dataframe[1:7,]
  dataframe$LT <- rep(1:7)
  dataframe$EM <- rep(ee, nrow(dataframe))
  
  dataframe$filename_prec <- paste0(bias_corrected_SEAS5,"/total_precipitation_LT",dataframe$LT,"_EM",dataframe$EM,"_MN",sprintf("%002d",dataframe$month), "_.nc")
  
  dataframe$filename_temp <- paste0(bias_corrected_SEAS5,"/2m_temperature_LT",dataframe$LT,"_EM",dataframe$EM,"_MN",sprintf("%002d",dataframe$month), "_.nc")
  
  # Applying the function to each row
  data <- dataframe %>%
    mutate(
      start_indices = mapply(get_monthly_indices, year, month)
    )
  data <- data %>%
    mutate(
      end_indices = mapply(get_monthly_indices_end, year, month)
    )
  
  data[data$year =="2017","start_indices"]<-1
  
  data[data$year =="2017","end_indices"]<-1+lubridate::days_in_month(data[data$year =="2017","month"])
  
  return(data)
  
}


# Function to get the number of days in the month
get_days_in_month <- function(year, month) {
  ymd(paste(year, month, "01", sep = "-")) %>% days_in_month()
}

library(lubridate)
get_details_past <- function(forecast_date, bias_corrected_ERA5){
  
  next_dates <- generate_dates_past(forecast_date)
  
  dataframe <- data.frame(date = next_dates,
                          month = lubridate::month(next_dates),
                          year = lubridate::year(next_dates)  )
  
  dataframe$filename_prec <-  paste0(bias_corrected_ERA5,"/total_precipitation_LT",0,"_EM",0,"_MN",sprintf("%002d",dataframe$month),"_.nc")
  
  dataframe$filename_temp <- paste0(bias_corrected_ERA5, "/2m_temperature_LT",0,"_EM",0,"_MN",sprintf("%002d",dataframe$month),"_.nc")
  
  # Applying the function to each row
  data <- dataframe %>%
    mutate(
      start_indices = mapply(get_monthly_indices_past, year, month)
    )
  data <- data %>%
    mutate(
      end_indices = mapply(get_monthly_indices_end_past, year, month)
    )
  
  # Add a new column for the number of days in the month
  data$days_in_month <- mapply(get_days_in_month, data$year, data$month)
  
  
  return(data)
  
}




library(terra)


# Define a function to process each row of the DataFrame
process_nc_file <- function(filename, start_idx, end_idx, start_date,obs,var) {
  # Read the NetCDF file
  filename2<- filename
  print(filename)
  rast <- rast(filename2)
  # terra::ext(rast) <- terra::ext(obs)
  # Extract the required layers based on start and end indices
  rast_subset <- rast[[as.numeric(start_idx):as.numeric(end_idx)]]
  print(as.character(start_date))
  # Parse the start date
  start_date <- as.Date(as.character(start_date))
  
  # Set the time for each layer based on the start date
  dates <- seq(start_date, by = "day", length.out = length(as.numeric(start_idx):as.numeric(end_idx)))
  
  # Save each layer as a separate TIFF file
  for (i in seq_along(dates)) {
    # Create a filename for each date
    date_str <- format(dates[i], "%Y%m%d")
    
    out_filename <- paste0(dates[i],"_SEAS5_", gsub(".nc", "",basename(filename2)), "_.tif")
    out_filenamenc <- paste0(dates[i],"_SEAS5_", gsub(".nc", "",basename(filename2)), "_.nc")
    # Write the single layer raster to a file
    writeRaster(rast_subset[[i]], filename = out_filename,  overwrite=TRUE)
    
    #terra::writeCDF(rast_subset[[i]], 
    #            zname = "time",
    #               varname = short_v,
    #           file =out_filenamenc,
    #           overwrite=T)
    
  }
}

process_nc_file_past <- function(filename, start_idx, end_idx, start_date,obs,var) {
  # Read the NetCDF file
  filename2<- as.character(filename)
  print(filename)
  rast <- rast(filename2)
  # terra::ext(rast) <- terra::ext(obs)
  # Extract the required layers based on start and end indices
  rast_subset <- rast[[as.numeric(start_idx):as.numeric(end_idx)]]
  print(start_date)
  print(as.character(start_date))
  # Parse the start date
  start_date <- as.Date(as.character(start_date))
  
  # Set the time for each layer based on the start date
  dates <- seq(start_date, by = "day", length.out = length(as.numeric(start_idx):as.numeric(end_idx)))
  
  # Save each layer as a separate TIFF file
  for (i in seq_along(dates)) {
    # Create a filename for each date
    
    date_str <- format(dates[i], "%Y%m%d")
    
    out_filename <- paste0(dates[i],"_ERA5_", gsub(".nc", "",basename(filename2)), "_.tif")
    
    out_filenamenc <- paste0(dates[i],"_ERA5_", gsub(".nc", "",basename(filename2)), "_.nc")
    
    print(out_filenamenc)
    
    writeRaster(rast_subset[[i]], filename = out_filename,  overwrite=TRUE)
    
  }
}





extract_chelsa <- function(dataset, precipitation_dir, temperature_dir, output_dir ){
  
  precipitation_files <- list.files(precipitation_dir, pattern = "nc", full.names=T)
  
  temperature_files <- list.files(temperature_dir, pattern = "nc", full.names=T)
  
  dataset <- dataset[dataset$year <2017,] # CHELSA temporal extent
  
  month_year_combos <- paste0(dataset$year, sprintf("%02d", dataset$month)) #Get unique dates
  
  dataset$month_year_combo <- month_year_combos
  
  all_together_combined <- NULL # Object to add extracted data too.
  
  for(u in 1:length(unique(month_year_combos))){ # For every unique month/year combination
    
    print(paste0(u, " of ", length(unique(month_year_combos))))
    
    combination <- unique(month_year_combos)[u]
    
    split_dataset <- dataset[dataset$month_year_combo == combination,]
    
    n_2 <- which(grepl(combination, precipitation_files))
    
    n_1 <- n_2 - 12 # Each nc file contains one month if data, we need one year 
    
    all_prec <- terra::rast(precipitation_files[n_1:n_2])
    
    n_2 <- which(grepl(combination,temperature_files))
    
    n_1 <- n_2 - 12
    
    all_temp <- terra::rast(temperature_files[n_1:n_2])
    
    days <- unique(split_dataset$day)
    
    
    for (d in 1:length(days)) {    # for every unique day in month/year combination
      
      split_dataset_further <- split_dataset[split_dataset$day == days[d],]
      
      dt <- paste0(split_dataset_further$year[1], "-",sprintf("%02d",split_dataset_further$month[1]),"-",sprintf("%02d",days[d]))
      
      n_3 <- which(grepl(dt, terra::time(all_temp)))
      
      n_4 <- n_3 - 365
      
      cropped_rast <- all_temp[[n_4:n_3]]
      
      split_dataset_further<-rbind(split_dataset_further,split_dataset_further)
      
      values <- terra::extract(cropped_rast,
                               y = as.matrix(split_dataset_further[, c("x", "y")]),
                               method = "simple")
      
      ann_mean_temp <- rowMeans(values, na.rm=T)   
      ann_sd_temp <- apply(values,1,sd, na.rm = T)   
      eight <- values[(length(values) - (7*8)):length(values)]   
      eight_mean_temp <- rowMeans(eight, na.rm=T)   
      eight_sd_temp <- apply(eight,1,sd, na.rm = T)
      
      
      n_3 <- which(grepl(dt, terra::time(all_prec)))
      
      n_4 <- n_3 - 365
      
      cropped_rast <- all_prec[[n_4:n_3]]
      
      values <- terra::extract(cropped_rast,
                               y = as.matrix(split_dataset_further[, c("x", "y")]),
                               method = "simple")
      
      
      ann_sum_prec <- rowSums(values, na.rm=T)   
      ann_sd_prec <- apply(values,1,sd, na.rm = T)    
      eight <- values[(length(values) - (7*8)):length(values)]   
      eight_sum_prec <- rowSums(eight, na.rm=T)   
      eight_sd_prec <- apply(eight,1,sd, na.rm = T)   
      
      
      extracted_data<- data.frame(mean_annual_temperature = ann_mean_temp,
                                  sd_annual_temperature = ann_sd_temp,
                                  mean_eight_temperature = eight_mean_temp,
                                  sd_eight_temperature =  eight_sd_temp,                               
                                  sum_annual_precipitation = ann_sum_prec,
                                  sum_eight_precipitation = eight_sum_prec)
      
      
      all_together <- cbind(split_dataset_further, extracted_data)
      
      all_together_combined <- rbind(all_together_combined, all_together)
      
    }
    
    
  }
  
  return(all_together_combined)
  
}


