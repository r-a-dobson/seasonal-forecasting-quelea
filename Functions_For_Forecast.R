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
# This script does not need to be opened. It contains custom functions for
# generating seasonal forecasts of red-billed quelea distribution. This file
# will be used to load functions in the main script "Seasonal_Forecast_Quelea.R"
# 
#' extract_daily_temperature Extracts daily temperature from SEAS5 .nc file
#' @param filepath path to nc file.
#' @param forecast_start forecast initiation date.
#' @param save_dir path to directory to save daily tifs to.
#' @param one_month logical, whether to extract one or seven months of data. 
 
extract_daily_temperature <- function(filepath, 
                                      forecast_start,
                                      save_dir,
                                      one_month = FALSE){
  
  

  # Get forecast end date 
  forecast_end <- forecast_start %m+% months(7)
  
  if (one_month) {
    forecast_end <- forecast_start %m+% months(1)
  }
  
  # Work out which layers to extract - every four intervals = one day for temp
  daysintheforecast <- as.numeric(forecast_end - forecast_start)
  
  endlayer <- as.numeric(daysintheforecast * 4)
  
  layer_numbers <- c(1)
  
  for (layer in 1:(endlayer / 4)) {
    layer_numbers <- c(layer_numbers, 4 * layer)
  }
  
  
  # Iterate through each layer and calculate the mean for day
  
  for (day in 1:((length(layer_numbers)-1))){
    
    print(paste(day,"/",daysintheforecast))
    
    # Get layers of nc file to read in for this day
    startlayer <- layer_numbers[day]
    endlayer <- layer_numbers[day + 1]
    
    # Empty list to bind ensemble members rasts too 
    ensembles <- vector("list", length = 51)
    
    # Iterate through each ensemble member to take the median for this day
    
    for (ensemble_member in 1:51) {
      
      r <- subset(brick(filepath, level = ensemble_member), startlayer:endlayer)
      
      r <- calc(r, mean, na.rm = T) # Mean across day to get average daily temperature
      
      ensembles[[ensemble_member]] <- r
      
    }
    
    # Calculate ensemble median value for this day
    r <- calc(raster::stack(ensembles), median, na.rm = T) 
    
    # Get the date for saving file
    daily_date <- forecast_start + day
    
    writeRaster(r,
                file = paste0(save_dir,
                              "/SEAS5_Dailymean_2m_air_temperature_",
                              daily_date,
                              "_.tif"),
                overwrite = T)
    
    listoffiles <- list.files(save_dir, full.names = T)
    
    file.remove(listoffiles[grep(".aux.xml",listoffiles)])
    
  }
  
}



#' extract_daily_precipitation Extracts daily precipitation from SEAS5 .nc file
#' @param filepath path to nc file.
#' @param forecast_start forecast initiation date.
#' @param save_dir path to directory to save daily tifs to.
#' @param one_month logical, whether to extract one or seven months of data. 

extract_daily_precipitation <- function(filepath, 
                                        forecast_start,
                                        save_dir,
                                        one_month = FALSE){
  
  # Get forecast end date 
  forecast_end <- forecast_start %m+% months(7)
  
  if(one_month) {
    forecast_end <- forecast_start %m+% months(1)
    
    forecast_end <- forecast_end +1
  }
  
  # Work out which layers to extract - every intervals = one day for temp
  daysintheforecast <- as.numeric(forecast_end - forecast_start)
  
  endlayer <- as.numeric(daysintheforecast)
  
  # Create temporary directory to save intermediate precipitation files to.
  # Precipitation is aggregate from forecast initiation, so must be subtracted 
  # from day after to get daily value. 
  
  processing_dir <- paste0(tempdir(),"/processing_prec")
  
  dir.create(processing_dir)
  
  for (day in 1:endlayer){
    
    print(paste("Initial processing: ",day,"/",daysintheforecast))
    
    ensembles <- vector("list", length = 51)
    
    # Extract daily value for each ensemble member
    for (ensemble_member in 1:51){
      
      r<-subset(brick(filepath,level=ensemble_member),day)
      
      ensembles[[ensemble_member]] <- r
    }
    
    # Take the ensemble median value for this day
    r <- calc(stack(ensembles), median, na.rm = T)
    
    daily_date <- forecast_start + day
    
    writeRaster(r,
                file = paste0(processing_dir,
                              "/SEAS5_Dailytotal_precipitation_",
                              daily_date,
                              "_.tif"),
                overwrite = T)
    
    listoffiles <- list.files(processing_dir, full.names = T)
    
    file.remove(listoffiles[grep(".aux.xml",listoffiles)])
    
  }
  
  # Get directory to save processed precipitation files to
  if(one_month) {
    processing_dir2 <- paste0(tempdir(), "/processing_prec2")
    dir.create(processing_dir2)
    file.remove(list.files(processing_dir2, full.names = T))
  }
  
  if (!one_month) {
    processing_dir2 <- save_dir
  }
  
  
  # List of daily precipitation summed from forecast initiation
  list <- list.files(processing_dir, full.names= T, pattern = paste0("_.tif"))
  
  
  # Working backwards from final forecast date...
  
  for (date in length(list):1){
    
    print(paste("Final processing: ",(length(list)-date)+1,"/",daysintheforecast))
    
    name <- list[date]
    
    name2 <- strsplit(name, "/")[[1]]
    name2 <- name2[(length(name2))]
    
    # For all but the final date, read in the raster and the raster of previous day.
    # Subtract one from the other to get the daily precipitation value
    
    if (!date == 1) {
      name_prev <- list[date - 1]
      
      r1 <- terra::rast(name)
      r2 <-  terra::rast(name_prev)
      
      r3 <- r1 - r2
    }
    
    # The first day has aggregated precipitation from the forecast initiation
    # (i.e. one day), so this does not need to be subtracted
    if (date == 1) {
      r3 <- terra::rast(name)
    }
    
    # Write the daily precipitation file to save directory 
    writeRaster(r3,file=paste0(processing_dir2,"/",name2),overwrite=T)
    
    listoffiles <- list.files(processing_dir2, full.names = T)
    
    file.remove(listoffiles[grep(".aux.xml",listoffiles)])
    
  }
  
  # Only copy across one-month tif files to the save directory if one_month = TRUE
  
  if(one_month){
    
    forecast_end <- forecast_start %m+% months(1)
    
    # Get one-month dates needed
    dates_to_keep <- dynamicSDM::dynamic_proj_dates(as.character(forecast_start),
                                                    as.character(forecast_end),
                                                    interval.level = "day",
                                                    interval = 1)
    
    list_of_files <- list.files(processing_dir2, full.names=T)
    
    selected_files <- list_of_files[grep(paste(dates_to_keep, collapse = "|"), list_of_files)]
    
    file.copy(selected_files,save_dir)
  }
  
  
  file.remove(list.files(processing_dir2, full.names = T))
  file.remove(list.files(processing_dir,full.names=T))
  
  
}  


#' resample_daily_weather Resample daily weather data to save spatial resolution
#' and extent.
#' @param daily_dir path to directory containing daily values for forecast and
#'   historical datasets.
#' @param spatial_resolution spatial resolution in degrees to resample weather
#'   data to.
#' @param spatial_extent sf polygon, the spatial extent to crop weather data to. 


resample_daily_weather <- function(daily_dir,
                                   spatial_resolution,
                                   spatial_extent){
  
 
  # Create rast for specified resolution and extent to resample weather data to.
  r <- terra::rast(spatial_extent)
  terra::res(r) <- spatial_resolution
  r <-  terra::setValues(r, 1:terra::ncell(r))
  r <- terra::crop(r, southernafrica)
  r <- terra::mask(r, southernafrica)

  
  # Iterate through each daily tif and resample/crop/mask so that they can
  # be stacked together. 
  
  list<-list.files(daily_dir, full.names=T)
  
  for (file in 1:length(list)){
    
    print(paste("Processing: ",file," / ",length(list)))
    r2<-terra::rast(list[file])
    
    r2<- terra::resample(r2,r)
    
    r2<-terra::crop(r2,r)
    
    terra::ext(r2)<-terra::ext(r)
    
    r2<-terra::mask(r2,r)
    
    writeRaster(r2,file=list[file],overwrite=T)
    
    listoffiles <- list.files(daily_dir, full.names = T)
    
    file.remove(listoffiles[grep(".aux.xml",listoffiles)])
    
  }
}


#' extract_proc_var Extracts processed 8- and 52- week weather variables for
#' forecast intervals
#' @param forecast_intervals character vector, the seasonal forecast interval
#'   dates to extract processed weather variables for
#' @param save_dir path to directory to save processed weather variables to. 
#' @param daily_dir path to directory containing forecast and historical daily values. 
#' @param period one of `eight` or `annual`, the temporal period to calculate
#'   variable across.
#' @param variable one of `precipitation` or `temperature` to specify which
#'   weather data to process.

extract_proc_var <- function(forecast_intervals,
                             save_dir,
                             daily_dir,
                             period,
                             variable ){
  
  # Specify number of days to calculate variable across
  if (period == "eight") {
    days <- 56
  }
  
  if (period == "annual") {
    days <- 365
  }

  # For every forecast interval, calculate the variable across period prior to
  # this date.
  
  for (interval in 1:length(forecast_intervals)){
    
    print(paste("Processing ", interval, "/", length(forecast_intervals)))
    
    date_1 <- forecast_intervals[interval]
    
    # List names of all daily files for this variable
    list <- list.files(daily_dir, full.names = T)
    list <- list[grep(variable, list)]
    
    # Get forecast file list
    list1 <- list[grepl("SEAS5", list)]

    # Get historical file list
    list2 <- list[grepl("ERA5", list)]
    
    list <- c(list2, list1)
    
    # Get location of first date tif within the list. e.g. 52-weeks prior
    initial <- list2[grepl(as.character(date_1 - days), list2)]
    
    # If not in the historical list, try seasonal list. 
    if (length(initial) == 0) {
      initial <- list1[grepl(as.character(date_1 - days), list1)]
    }
    
    # Get location of final date dir within the list. i.e. interval date.
    # This will definitely be in the seasonal dataset
    final <- list1[grepl(as.character(date_1), list1)]
    
    start <- as.numeric(match(initial, list))
    
    end <- as.numeric(match(final, list))
    
    # Get names of all daily files within this period
    files <- list[start:end]
    
    # Stack the rasters for the period
    all_daily <- terra::rast(files)
      
    # If temperature, then we want the mean and sd across period
    if (variable == "temperature") {
      mean_annual_temp <- terra::mean(all_daily,  na.rm = T)
      sd_annual_temp <- terra::stdev(all_daily,  na.rm = T)
      
      writeRaster(
        mean_annual_temp,
        file = paste0(save_dir, "/", date_1, "_mean_", period, "_temperature_.tif"),
        overwrite = T
      )
      writeRaster(
        sd_annual_temp,
        file = paste0(save_dir, "/", date_1, "_sd_", period, "_temperature_.tif"),
        overwrite = T
      )
    }
    # If precipitation, then we want the sum across period
    if (variable == "precipitation") {
      sum_annual_prec <- terra::app(all_daily, sum,  na.rm = T)
      
      writeRaster(
        sum_annual_prec,
        file = paste0(save_dir, "/", date_1, "_sum_", period, "_precipitation_.tif"),
        overwrite = T
      )
    }
    
    listtoremove <- list.files()
    file.remove(listtoremove[grep(".aux.xml", listtoremove)])
    
    
  }
  
}


#' get_evi_characterstics Extracts evi characterstics for vegetation growth
#' stage classification using 16-day MODIS Evi values from past 52-weeks.
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
#' @param model Random Forest classification model for vegetation growth stage based upon EVI characteristics. 
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
  
  # Now we have for each cell, the days since the EVI characteristics were
  # recorded, in which each vegetation growth stage likely began
  
  colnames(stage_track)<-c("x","y","midgreenup","peak","midgreen","dorm")
  
  forecast_doys<-as.numeric(forecast_intervals - (forecast_intitiation-14))
  

  # Iterate through each forecast interval, and if vegeation growth stage
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
    
    # Calculate sum across moving.window.matrix for the raster
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

