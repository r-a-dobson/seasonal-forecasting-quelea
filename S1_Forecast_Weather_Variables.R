#-------------------------------------------------------------------------------
# (1) Code for processing near-term hindcasts of weather variables
#------------------------------------------------------------------------------

# Please ensure that you have downloaded the GitHub repository that contains the
# associated data sets and custom functions.
# Found at: https://github.com/r-a-dobson/seasonal-forecasting-quelea
#
#------------------------------------------------------------------------------
# Step 1: Set up working directory and load packages
#------------------------------------------------------------------------------

# Provide path to "seasonal_forecasting_quelea" directory downloaded from GitHub
directory <- "C:/Users/XXXXX/Downloads/seasonal_forecasting_quelea/"

# Set this as your working directory
setwd(directory)

# Install packages required for this analysis
packages_to_install <- readLines("packages_to_install.txt")
install.packages(packages_to_install)

# If any packages installs fail, try installing directly from GitHub 
# install_github("SantanderMetGroup/downscaleR") # https://github.com/SantanderMetGroup/downscaleR

# Load required packages for this R script
library(dynamicSDM)
library(lubridate)
library(rnaturalearth)
library(terra)

# Read in the custom functions for near-term hindcasting
source(paste0(directory, "/", "Functions_For_Forecast.R"))

# Before running this script, please ensure that all required packages have been
# installed and that you have registered for free accounts on:
# > Climate Data Store (https://cds.climate.copernicus.eu/)
# > Google Earth Engine (https://earthengine.google.com/)
# > Google Drive  (https://workspace.google.com/intl/en_uk/products/drive/). 
# These are used to download hindcast and historical datasets.


#------------------------------------------------------------------------------
# Step 2: Set spatial extent and resolution for near-term hindcasts
#------------------------------------------------------------------------------

# Set the spatial extent and resolution for near-term distribution forecasting
# -	Spatial extent: country or countries to forecast across (or your own custom
# `sf` polygon). Please note all SDMs are trained on Quelea quelea lathamii
# records from southern Africa countries.
# 
# -	Spatial resolution: given in degrees, this value specifies the spatial
# resolution for the output seasonal forecasts.

countries_of_interest <- c('South Africa', # Can be replaced with any country 
                           'Botswana',
                           'Lesotho',
                           'Eswatini',
                           'Mozambique',
                           'Namibia',
                           'Zimbabwe',
                           'Angola',
                           'Zambia',
                           'Malawi')

spatial_extent <- ne_countries(country = countries_of_interest,
                               scale = "small",
                               returnclass = "sf")

spatial_extent_sp <-ne_countries(country = countries_of_interest,
                                 scale = "small",
                                 returnclass = "sp")

# Spatial resolution for forecasts in degrees
spatial_resolution <- 0.05

#------------------------------------------------------------------------------
# Step 3: Download SEAS5 Seasonal Forecast precipitation and temperature data
#------------------------------------------------------------------------------

#See details on seasonal forecast data at:
#https://www.ecmwf.int/en/newsletter/154/meteorology/ecmwfs-new-long-range-forecasting-system-seas5

library(ecmwfr)
library(keyring)
library(raster)

# Input your Climate Data Store account and log-in details
user_number <- "XXXXX" # Replace with your details
user_key <- "XXXXX" 
wf_set_key(user = user_number, key = user_key, service = "cds")

# Spatial extent to extract (CDS requires non-polygon format)
forecast_extent <- c(terra::ymax(terra::vect(spatial_extent)),
                     terra::xmin(terra::vect(spatial_extent)),
                     terra::ymin(terra::vect(spatial_extent)),
                     terra::xmax(terra::vect(spatial_extent)))

# Temporal extent to extract (training and testing period for bias correction) 

forecast_years <- 1982:2016

forecast_months <- 1:12

all_to_download <- expand.grid(forecast_years, forecast_months)

# Create directories to store SEAS5 hindcasts in
save_dir_tas_SEAS5 <- paste0(directory, "Data/SEAS5_TAS")
dir.create(save_dir_tas_SEAS5, recursive = T)

save_dir_pr_SEAS5 <- paste0(directory, "Data/SEAS5_PR")
dir.create(save_dir_pr_SEAS5, recursive = T)

# Iterate through each file required and download from CDS.

for(x in 1:nrow(all_to_download)){
  
  forecast_year <- all_to_download$Var1[x]
  
  forecast_months <- all_to_download$Var2[x]
  
#----------------------------
# Total precipitation
#----------------------------

input = 1:5160 # Number of hourly intervals in seven months
multiple_of_8 = (input %% 24) == 0 # Get the lead-time hour for each day
leadtimefordaily <- input[multiple_of_8]# Final lead-times to extract for daily 

# Parameters for CDS request - this shouldn't need changing
request <-
  list(
    "dataset_short_name"  = "seasonal-original-single-levels",
    "originating_centre"  = "ecmwf",
    "system"              = "5",
    "variable"            = "total_precipitation" ,
    "year"                = forecast_year,
    "month"               = forecast_month,
    "day"                 = "01",
    "leadtime_hour"       = leadtimefordaily,
    "format"              = "netcdf",
    "target"              = paste0("total_precipitation" , "_", forecast_year, "_",
                                          forecast_month, ".nc"), 
    "area"                = forecast_extent)

setwd(save_dir_pr_SEAS5)

# This sends the request to CDS. 
ncfile <- wf_request(user = user_number,
                     request = request,
                     transfer = TRUE,
                     time_out = 100000,
                     path = paste0(save_dir_pr_SEAS5),
                     verbose = T)

#----------------------------
# 2m temperature 
#----------------------------
input = 1:5160
multiple_of_8 = (input %% 6) == 0
leadtimefor6hour <- input[multiple_of_8]

# Parameters for CDS request - this shouldn't need changing
request <-
  list(
    "dataset_short_name" = "seasonal-original-single-levels",
    "originating_centre" = "ecmwf",
    "system"             = "5",
    "variable"           = "2m_temperature",
    "year"               = forecast_year,
    "month"              = forecast_month,
    "day"                = "01",
    "leadtime_hour"      = leadtimefor6hour,
    "format"             = "netcdf",
    "target"             = paste0("2m_temperature", "_", forecast_year, "_",
                                  forecast_month, ".nc"),
    "area"               = forecast_extent
  )

setwd(save_dir_tas_SEAS5)

# Send request to CDS and wait for download
ncfile <- wf_request(user = user_number,
                     request = request,   
                     transfer = TRUE, 
                     time_out= 100000,
                     path = paste0(save_dir_tas_SEAS5),
                     verbose = T)


}


#------------------------------------------------------------------------------
# Step 4: Extract daily data from near-term seasonal hindcasts
#------------------------------------------------------------------------------

#----------------------------------
# Process hindcast temperature data
#---------------------------------

seasonal_forecasts <- save_dir_tas_SEAS5

forecast_files <- list.files(seasonal_forecasts, full.name =T, pattern=".nc")

forecast_files <- forecast_files[grepl("temperature", forecast_files)]

# Create directory to store processed hindcast data
save_dir_tas_SEAS5_processed <- paste0(save_dir_tas_SEAS5, "/SEAS5_PROCCESSED")
dir.create(save_dir_tas_SEAS5_processed)

for(ens in 1:25){ # For each ensemble member (N = 25)
  
  ensemble_member <- ens
  
  for(x in 1:length(forecast_files)){ # For every SEAS5 hindcast downloaded
    
    print(x)
    
    forecast_data <- terra::rast(forecast_files[x])
    
    sequence <- seq(1, nlyr(forecast_data), by = 25) + (ensemble_member -1)
    
    forecast_data<-forecast_data[[sequence]]
    
    start_date <- as.Date(time(forecast_data[[1]]))
    
    end_date <- as.Date(time(forecast_data[[nlyr(forecast_data)]]))
    
    duration <- as.numeric(end_date - start_date)
    
    boundaries <- as.Date(c(start_date,
                            start_date %m+% months(1),
                            start_date %m+%  months(2),
                            start_date %m+%  months(3),
                            start_date %m+%  months(4),
                            start_date %m+%  months(5),
                            start_date %m+%  months(6),
                            end_date))
    
    all_dates <- time(forecast_data, format="days")
    
    for (day in 0:duration){
      
      print(day)
      curr_date <- start_date + day
      
      # Determine which dates in date_list the target_date falls between
      lead_time <- max(which(boundaries-1 < curr_date))
      
      if (day == 0) {
        lead_time <- 1
      }
      
      one_day <- forecast_data[[which(all_dates %in% (start_date + day))]]
      
      month <- sprintf("%002d",lubridate::month(curr_date))
      
      seas_fore <-save_dir_tas_SEAS5_processed
      
      dir.create(paste0(seas_fore, "/2m_temperature/LT_",lead_time,"/EM_",ensemble_member,"/MN_",month), recursive =T)
      
      if(!file.exists(paste0(seas_fore, "/2m_temperature/LT_",lead_time,"/EM_",ensemble_member,"/MN_",month,"/",
                             curr_date,"_tas_seas5hist_LT",lead_time,"_EM",ensemble_member,"_MT",month,"_.nc"))){
        
        setwd(paste0(seas_fore, "/2m_temperature/","LT_",lead_time,"/EM_",ensemble_member,"/MN_",month))
        
        r <-terra::mean(one_day,na.rm=T)
        
        terra::time(r)<-rep(curr_date, 1)
        
        names(r)<- rep("tas", length(names(r)))
        
        terra::writeCDF(r,
                        varname = "tas",
                        longname = "2 metre temperature",
                        unit = "K",
                        zname = "time",
                        file = paste0(curr_date,"_tas_seas5hist_LT",lead_time,"_EM",ensemble_member,"_MT",month,"_.nc"), overwrite=T)
        
      }
    }
  }
}

#--------------------------------------
# Process hindcast  precipitation data
#--------------------------------------

seasonal_forecasts <- save_dir_pr_SEAS5

forecast_files <- list.files(seasonal_forecasts, full.name = T, pattern = ".nc")

forecast_files <-forecast_files[grepl("precipitation",forecast_files)]

# Create directory to store processed hindcast data
save_dir_pr_SEAS5_processed <- paste0(save_dir_pr_SEAS5, "/SEAS5_PROCCESSED")
dir.create(save_dir_pr_SEAS5_processed)

for(ens in 1:25){ # For each ensemble member (N = 25)
 
  ensemble_member <- ens
  
  for(x in 1:length(forecast_files)){ # For every SEAS5 hindcast downlaoded
    
    forecast_data <- terra::rast(forecast_files[x])
    uni_dates<-unique(terra::time(forecast_data, format = "days"))
    
    sequence <- seq(1, nlyr(forecast_data), by = 25) + (ensemble_member -1)
    
    forecast_data<-forecast_data[[sequence]]
    
    forecast_data <- app(forecast_data, fun = function(x) { c(x[1], diff(x))  }) # pr is accumulated, so calculate the difference
    
    time(forecast_data) <- uni_dates
    
    start_date <- as.Date(time(forecast_data[[1]]))
    
    end_date <- as.Date(time(forecast_data[[nlyr(forecast_data)]]))
    
    duration <- as.numeric(end_date - start_date)
    
    boundaries <- as.Date(c(start_date,
                            start_date %m+% months(1),
                            start_date %m+%  months(2),
                            start_date %m+%  months(3),
                            start_date %m+%  months(4),
                            start_date %m+%  months(5),
                            start_date %m+%  months(6),
                            end_date))
    
    all_dates <- time(forecast_data, format="days")
    
    for (day in 0:duration){
      
      print(day)
      
      curr_date <- start_date + day
      
      # Determine which dates in date_list the target_date falls between
      lead_time <- max(which(boundaries-1 < curr_date))
      
      if(day == 0){
        lead_time <- 1}
      
      
      one_day <- forecast_data[[(day + 1)]]
      
      seas_fore <- save_dir_pr_SEAS5_processed
      
      month <- sprintf("%002d",lubridate::month(curr_date))
      
      if(!file.exists(paste0(seas_fore, "/total_precipitation/LT_",lead_time,"/EM_",ensemble_member,"/MN_",month,"/",
                             curr_date,"_pr_seas5hist_LT",lead_time,"_EM",ensemble_member,"_MT",month,"_.nc"))){
        
        dir.create(paste0(seas_fore, "/total_precipitation/LT_",lead_time,"/EM_",ensemble_member,"/MN_",month), recursive =T)
        
        setwd(paste0(seas_fore, "/total_precipitation/","LT_",lead_time,"/EM_",ensemble_member,"/MN_",month))
        
        names(one_day)<- rep("pr", length(names(one_day)))
        
        print(paste0("writing", curr_date))
        
        terra::writeCDF(one_day,
                        varname = "pr",
                        longname = "2 metre temperature",
                        unit = "K",
                        zname = "time",
                        file = paste0(curr_date,"_pr_seas5hist_LT",lead_time,"_EM",ensemble_member,"_MT",month,"_.nc"), overwrite=T)
      }
    }
  }
}

#------------------------------------------------------------------------------
# Step 5: Bias correction of seasonal hindcasts  
#------------------------------------------------------------------------------
# Bias correction is split by variable, lead-time, month, and ensemble member.

# First, you need to download CHELSA-W5E5 daily temperature and precipitation
# data for 1982-2016 from:
# https://chelsa-climate.org/chelsa-w5e5-v1-0-daily-climate-data-at-1km-resolution/

precipitation_dir <- "XXX" # Replace with path to CHELSA-W5E5 precipitation files 

temperature_dir <- "XXX" # Replace with path to CHELSA-W5E5 temperature files

# Load packages for bias correction using the DownscaleR package 
options(java.parameters = "-Xmx10g")
library(reticulate)
library(terra)
library(parallel)
library(ncdf4)
library(loadeR)
library(loadeR.java)
library(climate4R.UDG)
library(climate4R.datasets)
library(transformeR)
library(downscaleR)
library(visualizeR)
library(rJava)
.jinit() 
library(reticulate)
library(loadeR.2nc)
library(ncdf4)

variables <- c("total_precipitation","2m_temperature")

LTS <- c(1:7)

MNS <- c(1:12)

EMS <- c(1:25)

save_directory_bc <- paste0(directory,"/Data/SF_BC/")

setwd(save_directory_bc)

train_years <- c(1982:2001)

test_years <- c(2002:2016)

combinations <- expand.grid(variables, LTS, MNS, EMS)

# Iterate through each variable, lead-time, month and ensemble member
for (x in 1:nrow(combinations)) {
  
  print(paste(x, " of ", nrow(combinations)))
  
  VAR1 <- as.character(combinations$Var1[x])
  
  LTS1 <- as.character(combinations$Var2[x])
  
  MNS1 <- sprintf('%02d', combinations$Var3[x])
  
  EMS1 <- as.character(combinations$Var4[x])
  
  if(!file.exists(paste0(save_directory_bc, "/ ",VAR1,"_LT",LTS1,"_EM",EMS1,"_MN",MNS1,"_.nc"))){
    
    ###############################
    # Read in the forecast data 
    ################################
   
    if(VAR1 == "total_precipitation"){
      
      director <- paste0(save_dir_pr_SEAS5_processed, "/", VAR1, "/LT_", LTS1, "/EM_", EMS1, "/MN_", MNS1)
      
      makeAggregatedDataset(source.dir = director,
                            ncml.file = paste0(x,"combined_data.ncml"))
      short_v <- "pr"
      short_vera <- "tp"
      UNIT1 <- "MM"
    }
    
    if(VAR1 == "2m_temperature"){
      
      director <- paste0(save_dir_tas_SEAS5_processed, "/", VAR1, "/LT_", LTS1, "/EM_", EMS1, "/MN_", MNS1)
      
      makeAggregatedDataset(source.dir = director,
                            ncml.file = paste0(x,"combined_data.ncml"))
      short_v <- "tas"
      short_vera <- "t2m"
      UNIT1 <- "K"
    }
    
    print("reading in the training forecast")
    training_grid <- loadGridData(paste0(x,"combined_data.ncml"), var = short_v, years = train_years, season = combinations$Var3[x])
    
    print("reading in the testing forecast")
    testing_grid <- loadGridData(paste0(x,"combined_data.ncml"), var = short_v, years = test_years, season = combinations$Var3[x])
    
    
    ###############################
    # Load the CHELSA-W5E5 data 
    ################################
    
    print("Loading CHELSA-W5E5")
    
    if(VAR1 == "total_precipitation"){
      
      makeAggregatedDataset(source.dir = precipitation_dir,
                            ncml.file = paste0("historical_CHELSA_pr.ncml"))
      
      historical_grid <- loadGridData("historical_CHELSA_pr.ncml", var = short_v, years = train_years, season = combinations$Var3[x])
      tp_daily <- historical_grid
      tp_daily$Variable$varName <- "pr"
      
    }
    
    if(VAR1 == "2m_temperature"){
      
      makeAggregatedDataset(source.dir = temperature_dir,
                            ncml.file = paste0("historical_CHELSA_tas.ncml"))
      
      historical_grid <- loadGridData("historical_CHELSA_tas.ncml", var = short_v, years = train_years, season = combinations$Var3[x])
      tp_daily <- historical_grid
      tp_daily$Variable$varName <- "tas"
      
    }
    
    print("Bias correction")
    
    # Regrid CHELSA to SEAS5 
    new.coordinates <- getGrid(training_grid)
    
    tp_daily <- interpGrid(tp_daily, new.coordinates = new.coordinates,
                 method = "bilinear")
    
    # Run bias correction
    if(VAR1 == "total_precipitation"){
      bc_data <- biasCorrection(y = tp_daily, x = training_grid, newdata = testing_grid, method = "eqm", precip=T,cross.val = "none")
    }
    
    if(!VAR1 == "total_precipitation"){
      bc_data <- biasCorrection(y = tp_daily, x = training_grid, newdata = testing_grid, method = "eqm", precip=F,cross.val = "none")
    }
    
    
    bc_data_proc <- lapply(1:nrow(bc_data$Data),  function(x, bc = bc_data) {
      matrix <- bc$Data[x, ,] 
      terra::rast(matrix[nrow(matrix):1, ])})
    
    bc_data_proc <- terra::rast(bc_data_proc)
    
    # Set extent of BC data
    get_grid <- terra::rast(list.files(director, full.names=T)[1])
    
    terra::ext(bc_data_proc) <-  terra::ext(get_grid)
    
    bc_data_proc <-  terra::resample(bc_data_proc,get_grid )
    
    print(length(as.Date(bc_data$Dates$start)))
    
    terra::time(bc_data_proc) <- as.Date(bc_data$Dates$start)
    
    # Create save directory for bias corrected SEAS5
    dir.create(paste0(save_directory_bc,"/"))
    
    names(bc_data_proc) <- rep(short_v, length(names(bc_data_proc)))
    
    terra::writeCDF(bc_data_proc,
                    varname = short_v,
                    longname = VAR1,
                    unit = UNIT1,
                    zname = "time",
                    file = paste0(save_directory_bc, "/ ",VAR1, "_LT", LTS1, "_EM", EMS1, "_MN", MNS1, "_.nc"), overwrite=T)
    
    file.remove(paste0(x, "combined_data.ncml"))
    file.remove(list.files(tempdir(), recursive = T, full.names = T))
  }
}

#------------------------------------------------------------------------------
# Step 6: Download ERA5 precipitation and temperature data 
#------------------------------------------------------------------------------

# Download ERA5 re-analysis data from the CDS. For information on the dataset see:
# https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5

library(ecmwfr)
library(keyring)
library(raster)

years<-c(1979:2016)

variable <- c("total_precipitation","2m_temperature") 

combinations <- as.data.frame(expand.grid(years, variable))

save_directory <- paste0(directory, "Data/ERA5")

dir.create(save_directory, recursive = T)

for (c in 1:nrow(combinations)){
  
  print(c)
  
  setwd(save_directory)
  
  yer <- combinations$Var1[c]
  
  var <- as.character(combinations$Var2[c])
  
  dir.create(var)
  
  variable <- var
  
  if (var == "total_precipitation") {
    variable <- "tp"
  }
  
  path <- paste0(save_directory, variable)
  
  setwd(path)
  
  
  if(!file.exists(paste0(var, "_", yer, "_.nc"))){
    
    request <- list("dataset_short_name" = 'reanalysis-era5-single-levels',
                    'product_type'       = 'reanalysis',
                    "variable"           = var,
                    "year"               = yer,
                    "month"              = c('01',
                                             '02',
                                             '03',
                                             '04',
                                             '05',
                                             '06',
                                             '07',
                                             '08',
                                             '09',
                                             '10',
                                             '11',
                                             '12'), 
                    "day"                = c('01', '02', '03',
                                             '04', '05', '06',
                                             '07', '08', '09',
                                             '10', '11', '12',
                                             '13', '14', '15',
                                             '16', '17', '18',
                                             '19', '20', '21',
                                             '22', '23', '24',
                                             '25', '26', '27',
                                             '28', '29', '30',
                                             '31'),
                    'time'               = c('00:00', '01:00', '02:00',
                                             '03:00', '04:00', '05:00',
                                             '06:00', '07:00', '08:00',
                                             '09:00', '10:00', '11:00',
                                             '12:00', '13:00', '14:00',
                                             '15:00', '16:00', '17:00',
                                             '18:00', '19:00', '20:00',
                                             '21:00', '22:00', '23:00'),
                    "grid"               = "1.0/1.0",
                    "format"             = "netcdf",
                    "target"             = paste0(var, "_", yer, "_.nc"), 
                    'area'               = forecast_extent
    )
    
    
    ncfile <- wf_request(
      user = "127243",
      request = request,
      transfer = TRUE,
      time_out = 100000000,
      path = path,
      verbose = T
    )
    
  }
}

#------------------------------------------------------------------------------
# Step 7: Bias correction of ERA5 data 
#------------------------------------------------------------------------------

output_directory <- paste0(directory, "/Data/ERA5_BC/")

variables <- c("total_precipitation", "2m_temperature")

LTS <- c(0) 

MNS <- c(1:12)

EMS <- c(0)

precipitation_ERA5 <- ""#  Directory containing ERA5 data, must be stored in folder named "tp"  

temperature_ERA5 <- "" #  Directory containing ERA5 data, must be stored in folder named "t2m" 

train_years <- c(1982:2001)

test_years <- c(2002:2016)

combinations <- expand.grid(variables, LTS, MNS, EMS)


for (x in 1:nrow(combinations)) {
  
  print(paste(x, " of ", nrow(combinations)))
  
  VAR1 <- as.character(combinations$Var1[x])
  
  LTS1 <- as.character(combinations$Var2[x])
  
  MNS1 <- sprintf('%02d', combinations$Var3[x])
  
  EMS1 <- 0
  
  if(!file.exists(paste0(output_directory, VAR1, "_LT", LTS1, "_EM", EMS1, "_MN", MNS1, "_.nc"))) {
    
    if (VAR1 == "total_precipitation") {
      short_v <- "pr"
      short_vera <- "tp"
      UNIT1 <- "MM"
    }
    
    if (VAR1 == "2m_temperature") {
      short_v <- "tas"
      short_vera <- "t2m"
      UNIT1 <- "K"
    }
    
    
    if (VAR1 == "total_precipitation") {
      
      reanalysis_directory <- precipitation_ERA5
      
      print("reading in the training ERA")
      
      if (!file.exists("historical_agg_prec.ncml")) {
        makeAggregatedDataset(source.dir = reanalysis_directory,
                              ncml.file = "historical_agg_prec.ncml")
      }
      
      # Sum hourly to daily resolution
      training_grid <- loadGridData("historical_agg_prec.ncml", var = short_vera, years = train_years, season = combinations$Var3[x])
      training_grid <- aggregateGrid(training_grid, aggr.d = list(FUN = "sum", na.rm = TRUE))
      training_grid$Variable$varName <- "pr"
   
      testing_grid <- loadGridData("historical_agg_prec.ncml", var = short_vera, years = test_years, season = combinations$Var3[x])
      testing_grid <- aggregateGrid(testing_grid , aggr.d = list(FUN = "sum", na.rm = TRUE))
      testing_grid$Variable$varName <- "pr"
      
      
    }
    
    if(VAR1 == "2m_temperature"){
      
      reanalysis_directory <- temperature_ERA5
      print("reading in the training ERA")
      
      if (!file.exists("historical_agg_temp.ncml")) {
        makeAggregatedDataset(source.dir = reanalysis_directory,
                              ncml.file = "historical_agg_temp.ncml")
      }
      
      # Average hourly to daily resolution
      training_grid <- loadGridData("historical_agg_temp.ncml", var = short_vera, years = train_years, season = combinations$Var3[x])
      training_grid <- aggregateGrid(training_grid, aggr.d = list(FUN = "mean", na.rm = TRUE))
      training_grid$Variable$varName<- "tas"
      
      testing_grid <- loadGridData("historical_agg_temp.ncml", var = short_vera, years = test_years, season = combinations$Var3[x])
      testing_grid<- aggregateGrid(testing_grid, aggr.d = list(FUN = "mean", na.rm = TRUE))
      testing_grid$Variable$varName <- "tas"
      
    }
    
    ###############################
    # Load CHELSA-W5E5 data 
    ################################
    
    print("Loading CHELSA-W5E5 data")
    
    if(VAR1 == "total_precipitation"){
      historical_grid <- loadGridData("historical_CHELSA_pr.ncml", var = short_v, years = train_years, season = combinations$Var3[x])
      tp_daily <- historical_grid
      tp_daily$Variable$varName <- "pr"
    }
    
    if(VAR1 == "2m_temperature"){
      historical_grid <- loadGridData("historical_CHELSA_tas.ncml", var = short_v, years = train_years, season = combinations$Var3[x])
      tp_daily <- historical_grid
      tp_daily$Variable$varName <- "tas"
    }
    
    # Regrid CHELSA to ERA5 
    new.coordinates <- getGrid(training_grid)
    tp_daily <- interpGrid(tp_daily, new.coordinates = new.coordinates,
                           method = "bilinear")
    
    # Correct bias
    if(VAR1 == "total_precipitation"){
      bc_data <- biasCorrection(y = tp_daily, x = training_grid, newdata = testing_grid, method = "eqm", precip = T, cross.val = "none")
    }
    
    if(!VAR1 == "total_precipitation"){
      bc_data <- biasCorrection(y = tp_daily, x = training_grid, newdata = testing_grid, method = "eqm", precip = F, cross.val = "none")
    }
    
    
    bc_data_proc <-lapply(1:nrow(bc_data$Data), function(x, bc = bc_data) {
      matrix <- bc$Data[x, ,]  # Extract the layer's data
      terra::rast(matrix[nrow(matrix):1, ])
    })
    
    bc_data_proc <- terra::rast(bc_data_proc)
    
    # Extent sort
    
    get_grid <- terra::rast(list.files(reanalysis_directory, full.names = T)[1])
    
    terra::ext(bc_data_proc) <-  terra::ext(get_grid)
    
    bc_data_proc <-  terra::resample(bc_data_proc,get_grid )
    
    print(length(as.Date(bc_data$Dates$start)))
    
    terra::time(bc_data_proc) <- as.Date(bc_data$Dates$start)
    
    dir.create(paste0(output_directory))
    
    names(bc_data_proc) <- rep(short_v, length(names(bc_data_proc)))
    
    terra::writeCDF(bc_data_proc,
                    varname = short_v,
                    longname = VAR1,
                    unit = UNIT1,
                    zname = "time",
                    file =paste0(output_directory, VAR1, "_LT", LTS1, "_EM", EMS1, "_MN", MNS1, "_.nc"), overwrite = T)
    
    file.remove(paste0(x, "combined_data.ncml"))
    file.remove(list.files(tempdir(), recursive = T, full.names = T))
    
  }
}

#------------------------------------------------------------------------------
# Step 8: Combine ERA5 and SEAS5 to generate weather variable hindcasts
#------------------------------------------------------------------------------

# Load required packages
library(readr)
library(terra)
library(pROC)
library(dplyr)
library(magrittr)
library(rnaturalearth)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(lubridate)

working_dir <- directory

years <- c(2004:2016)

months <- c(1:12)

days <- c(1)

dataframe <- as.data.frame(expand.grid(years, months, days))

colnames(dataframe) <- c("year", "month", "day")

dataframe <- dataframe[order(dataframe$year), ]

dataframe$date <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d")

#---------------------------------------------------------------
# Create data frame of relevant dates for near-term hindcasting
#---------------------------------------------------------------
dataframe$forecastdate_1 <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m+% months(1)
dataframe$forecastdate_2 <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m+% months(2)
dataframe$forecastdate_3 <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m+% months(3)
dataframe$forecastdate_4 <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m+% months(4)
dataframe$forecastdate_5 <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m+% months(5)
dataframe$forecastdate_6 <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m+% months(6)
dataframe$forecastdate_7 <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m+% months(7)

dataframe$minus_month_one   <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m-% months(1)
dataframe$minus_month_two   <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m-% months(2)
dataframe$minus_month_three <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m-% months(3)

# Historical data period
dataframe$historical_time_start <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") - 365
dataframe$historical_time_end <- as.Date(with(dataframe, paste(year, month, day, sep = "-")), "%Y-%m-%d") %m-% months(3) -  1

Forecasting_Dates <- as.data.frame(dataframe)

# Specify directories containing bias corrected SEAS5 and ERA5 data
bias_corrected_ERA5 <- paste0(directory, "/Data/ERA5_BC/")

bias_corrected_SEAS5 <- paste0(directory, "/Data/SF_BC/")

for (f in 1:nrow(Forecasting_Dates)) {
  
  row <- Forecasting_Dates[f, ]
  
  output_dir_this <- paste0(working_dir, row$date, "/combined_weather/")
  
  backupfilename <- row$date
  
  print(row)
  
  working_dir_curr <- paste0(working_dir, row$date, "/daily_raster", recursive = T)
  
  output_dir_current <- paste0(working_dir, row$date, "/ensemble_medians/")
 
  file.remove(list.files(paste0(working_dir, row$date), full.names = T, recursive = T))
      
  working_dir_curr <- paste0(working_dir, row$date, "/daily_raster", recursive = T)
  
  dir.create(working_dir_curr, recursive = T)
      
      for(ens in 1:25){
        
        output_dir_this <- paste0(working_dir, row$date, "/combined_weather/")
        
        dir.create(output_dir_this, recursive = T)
        
        output_dir_curr <- paste0(working_dir,row$date,"/processing_ensemble/")
        
        dir.create(output_dir_curr, recursive = T)
        
        if(!file.exists(paste0(output_dir_this,"/", row$forecastdate_7, "_", ens, "_eight_sum_precipitation_", "14", "_.tif"))){
      
          if(ens==1){
            
            data <- get_details_past(row$date, bias_corrected_ERA5)
            
            setwd(working_dir_curr)
            
            for (x in 1:nrow(data)) {
              
              process_nc_file_past(data$filename_prec[x],
                                   data$start_indices[x],
                                   data$end_indices[x],
                                   data$date[x],
                                   obs)
              
            }
          
            for (x in 1:nrow(data)) {
              process_nc_file_past(data$filename_temp[x],
                                   data$start_indices[x],
                                   data$end_indices[x],
                                   data$date[x],
                                   obs)
            }
          }
          
          # Get resolution and extent
          r <- raster(spatial_extent_sp)
          res(r) <- spatial_resolution
          r <-  raster::setValues(r, 1:ncell(r))
          r <- crop(r, spatial_extent_sp)
          r <- mask(r, spatial_extent_sp)
          
          data <- get_details(row$date, ens = ens, bias_corrected_SEAS5)
          
          setwd(working_dir_curr)
          
          # Apply function to process files for precipitation
          apply(data, 1, function(x)
            process_nc_file(x["filename_prec"], x["start_indices"], x["end_indices"], x["date"], obs))
          
     
          # Apply function to process files for temperature
          apply(data, 1, function(x)
            process_nc_file(x["filename_temp"], x["start_indices"], x["end_indices"], x["date"], obs))
          
          
          setwd(working_dir_curr)
          
          list <- list.files(pattern = ".tif$")
          
          list <- list[grep("temperature", list)]
          
          listx <- list.files(pattern = ".tif$")
          
          listx <- listx[grep("precipitation", listx)]
          
          if(ens ==1){
            lapply(1:length(list), function(x) {
              terra::writeRaster(
                terra::resample(terra::rast(list[x]), y = terra::rast(list[length(list)])),
                file = list[x],
                overwrite = T
              )
            })
            
            lapply(1:length(listx), function(x) {
              terra::writeRaster(
                terra::resample(terra::rast(listx[x]), y = terra::rast(listx[length(listx)])),
                file = listx[x],
                overwrite = T
              )
            })
          }
          
          dates<-c(as.Date(paste0(lubridate::year(row$date),"-",lubridate::month(row$date),"-14")),
                   row$forecastdate_1,
                   as.Date(paste0(lubridate::year(row$forecastdate_1),"-",lubridate::month(row$forecastdate_1),"-14")),
                   row$forecastdate_2,
                   as.Date(paste0(lubridate::year(row$forecastdate_2),"-",lubridate::month(row$forecastdate_2),"-14")),
                   row$forecastdate_3,
                   as.Date(paste0(lubridate::year(row$forecastdate_3),"-",lubridate::month(row$forecastdate_3),"-14")),
                   row$forecastdate_4,
                   as.Date(paste0(lubridate::year(row$forecastdate_4),"-",lubridate::month(row$forecastdate_4),"-14")),
                   row$forecastdate_5,
                   as.Date(paste0(lubridate::year(row$forecastdate_5),"-",lubridate::month(row$forecastdate_5),"-14")),
                   row$forecastdate_6,
                   as.Date(paste0(lubridate::year(row$forecastdate_6),"-",lubridate::month(row$forecastdate_6),"-14")),
                   row$forecastdate_7)
          
          Dates <- dates
        
          for (dt in 1:length(dates)) {
            
            print(dt)
            
            setwd(working_dir_curr)
            
            savename<- dates[dt]
            
            curr_date<-dates[dt]
            
            if (dt == 14) {
              curr_date <- curr_date - 1
            }
            
            current_date <- which(grepl(curr_date, list))
            
            if (length(current_date) == 0) {
              current_date <- which(grepl(curr_date - 1, list))
            }
            
            annual_indices <- current_date - 365
            
            eight_indices <- current_date - 56
            
            annual_temp_rasts <- terra::rast(list[annual_indices:(current_date-1)])
            
            eight_temp_rasts <- terra::rast(list[eight_indices:(current_date-1)])
            
            # Calculate 8-week and 52-week weather variables from data
            
            ann_mean_temp <- terra::mean(annual_temp_rasts, na.rm = T)
            ann_sd_temp <- terra::app(annual_temp_rasts, "sd", na.rm = T)
            
            eight_mean_temp <- mean(eight_temp_rasts, na.rm = T)
            eight_sd_temp <- terra::app(eight_temp_rasts, "sd", na.rm = T)
            
            current_date <- which(grepl(curr_date, listx))
            
            if (length(current_date) == 0) {
              current_date <- which(grepl(curr_date - 1, list))
            }
            
            annual_indices <- current_date - 365
            
            eight_indices <- current_date - 56
         
            annual_prec_rasts <- terra::rast(listx[annual_indices:(current_date-1)])
            
            eight_prec_rasts <- terra::rast(listx[eight_indices:(current_date-1)])
            
            eight_sum_prec <- terra::app(eight_prec_rasts, "sum", na.rm = T)
            
            annual_sum_prec <- terra::app(annual_prec_rasts, "sum", na.rm = T)
            
            setwd(output_dir_this)
            
            writeRaster(ann_mean_temp, paste0(savename,"_", ens, "_annual_mean_temperature_", dt, "_.tif"), overwrite = T)
            writeRaster(ann_sd_temp, paste0(savename, "_", ens, "_annual_sd_temperature_", dt, "_.tif"), overwrite = T)
            writeRaster(eight_mean_temp, paste0(savename, "_", ens, "_eight_mean_temperature_", dt, "_.tif"), overwrite = T)
            writeRaster(eight_sd_temp, paste0(savename,"_", ens, "_eight_sd_temperature_", dt, "_.tif"), overwrite = T)
            writeRaster(annual_sum_prec, paste0(savename, "_", ens, "_annual_sum_precipitation_", dt, "_.tif"), overwrite = T)
            writeRaster(eight_sum_prec, paste0(savename, "_", ens, "_eight_sum_precipitation_", dt, "_.tif"), overwrite = T)
         
          }
          
          file.remove(list.files(output_dir_curr, full.names = T))
          file.remove(list.files(working_dir_curr, full.names = T, pattern = "SEAS5"))
          
        }
        
        # Now processed for all ensemble members, need to calculate the median
        
        dates<-c(as.Date(paste0(lubridate::year(row$date), "-", lubridate::month(row$date), "-14")), 
                 row$forecastdate_1,
                 as.Date(paste0(lubridate::year(row$forecastdate_1), "-", lubridate::month(row$forecastdate_1), "-14")), 
                 row$forecastdate_2,
                 as.Date(paste0(lubridate::year(row$forecastdate_2), "-", lubridate::month(row$forecastdate_2), "-14")), 
                 row$forecastdate_3,
                 as.Date(paste0(lubridate::year(row$forecastdate_3), "-", lubridate::month(row$forecastdate_3), "-14")), 
                 row$forecastdate_4,
                 as.Date(paste0(lubridate::year(row$forecastdate_4), "-", lubridate::month(row$forecastdate_4), "-14")), 
                 row$forecastdate_5,
                 as.Date(paste0(lubridate::year(row$forecastdate_5), "-", lubridate::month(row$forecastdate_5), "-14")), 
                 row$forecastdate_6,
                 as.Date(paste0(lubridate::year(row$forecastdate_6), "-", lubridate::month(row$forecastdate_6), "-14")), 
                 row$forecastdate_7)
        
        Dates <- dates
        vars <-
          c("_annual_mean_temperature_",
            "_annual_sd_temperature_",
            "_eight_mean_temperature_",
            "_eight_sd_temperature_",
            "_annual_sum_precipitation_",
            "_eight_sum_precipitation_")
        
        output_dir_current <- paste0(working_dir, row$date, "/ensemble_medians/")
        
        dir.create(output_dir_current, recursive = T)
        
        dir.create(paste0(tempdir(), "/proc/"), recursive = T)
        
        for (date_x in 1:length(dates)){
          
          print(date_x)
          
          setwd(output_dir_this)
          
          savename<- dates[date_x]
          
          curr_date <- dates[date_x]
          
          for(v in 1:length(vars)){
            
            filelist <- list.files(full.names=T)
            
            filelist <- filelist[grepl(vars[v], filelist)]
            
            filelist <- filelist[grepl(curr_date, filelist)]
            
            ensemble_median <- terra::app(terra::rast(filelist), median, na.rm = T)
            
            terra::time(ensemble_median) <- as.Date("2000-01-02") # Assign generic time as not relevant
            
            writeRaster(ensemble_median, file = paste0(output_dir_current, curr_date, "_", vars[v], "_.tif"), overwrite = T)
            
          }
        }
      }
  }