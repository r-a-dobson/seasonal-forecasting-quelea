################################################################################
#                                                                              #
# Generate seasonal forecasts for the red-billed quelea (Quelea quelea)        #
#                                                                              #    
# Code for: Seasonal forecasting of mobile species geographical distributions  #
# for dynamic management under extreme weather events.                         #
#                                                                              #
# Rachel Dobson, Stephen G. Willis, Stewart Jennings, Robert A. Cheke,         #
# Andrew J. Challinor and Martin Dallimer                                      #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------
# Step 1: Set up and seasonal forecast specification
#------------------------------------------------------------------------------
# This script includes the code for every step in generating seasonal forecasts
# of red-billed quelea distribution suitability. Please ensure that you have
# downloaded the Zenodo folder that contains associated data sets and custom
# functions written for this script.
#
# Our code is highly flexible. In this first step, you can specify the: 
#
# -	Spatial extent: country or countries to forecast across (or your own custom
# `sf` polygon). Please note all SDMs are trained on Quelea quelea lathamii
# records from southern Africa countries.
# 
# -	Spatial resolution: given in degrees, this value specifies the spatial
# resolution for the output seasonal forecasts.
#
# -	Temporal extent: the month and year to initiate the seasonal forecast. ECMWF
# SEAS5 forecasts have an initial date of the 1st of each month, and run for 7
# months. These data are released on the 5th of each month. We recommend
# choosing the closest available month.
#
# -	Temporal resolution: the intervals to generate seasonal forecasts at between
# zero and seven months ahead of the initiation date. This can range from daily
# to weekly to monthly intervals depending on your needs.
#
# Once you have specified these details in the first step, the rest of the code
# does not need to be edited aside from giving your Climate Data Store (CDS)
# username and key (Step 2), and your Google Earth Engine/Google Drive email
# (Step 4).
#
# Final output of this script are quelea distribution suitability forecasts
# for each interval within the seasonal forecast. This is exported in both tif
# format and png image of plotted suitability maps in ggplot2.

# Provide path to "seasonal_forecasting_quelea" directory downloaded from Zenodo
directory <- "C:/Users/xxxxx/Downloads/seasonal_forecasting_quelea/"

# Set this directory as your working directory for analyses
setwd(directory)

# Before beginning the script, please ensure that all necessary packages have
# been installed and that you have registered for free Climate Data Store,
# Google Earth Engine and Google Drive accounts. These are required for download
# of forecast and historical data sets throughout the seasonal forecasting
# workflow.

packages_to_install <- readLines("packages_to_install.txt")

install.packages(packages_to_install)


# Temporal extent to extract 
forecast_year<- "2024"
forecast_month<- "1"

# Create forecast initiation date
forecast_intitiation <- as.Date(paste0(forecast_year, "-",forecast_month,"-01"))

# Forecast intervals - here, fortnightly between 0 and 7 months
library(dynamicSDM)
library(lubridate)

forecast_intervals <- dynamic_proj_dates(startdate = as.character(forecast_intitiation),
                                         enddate =    as.character(forecast_intitiation %m+% months(7)),  
                                         interval = 2, # Represents two of each interval i.e. fortnightly
                                         interval.level = "week") # Change to day or month if desired

# Example for monthly forecasts

#forecast_intervals <- dynamic_proj_dates(startdate = as.character(forecast_intitiation),
#                                         enddate =    as.character(forecast_intitiation %m+% months(7)),  
#                                         interval = 1, 
#                                         interval.level = "month") 


# Spatial extent 
library(rnaturalearth)
# Alternatively provide an `sf` polygon for specific region
spatial_extent <- ne_countries(country = c('South Africa', # Any country globally could be given here
                                           'Botswana',
                                           'Lesotho',
                                           'Swaziland',
                                           'Mozambique',
                                           'Namibia',
                                           'Zimbabwe',
                                           'Angola',
                                           'Zambia',
                                           'Malawi'),
                               scale = "small",
                               returnclass = "sf")

# Spatial resolution for forecasts in degrees
spatial_resolution <- 0.05


# Read in the custom functions written for seasonal forecasting
source(paste0(directory, "/", "Functions_For_Forecast.R"))

#------------------------------------------------------------------------------
# Step 2: Download SEAS5 Seasonal Forecast precipitation and temperature data
#------------------------------------------------------------------------------
library(ecmwfr)
library(keyring)
library(raster)

#----------------------------
# Set key extraction variables 
#----------------------------

# Get Climate Data Store account and set log-in details
user_number <- "" # Replace with your details
user_key <- "" 
wf_set_key(user = user_number, key = user_key, service = "cds")

# Spatial extent to extract (CDS requires non-polygon format)
forecast_extent <- c(terra::ymax(terra::vect(spatial_extent)),
                     terra::xmin(terra::vect(spatial_extent)),
                     terra::ymin(terra::vect(spatial_extent)),
                     terra::xmax(terra::vect(spatial_extent)))


#----------------------------
# Extract total precipitation
#----------------------------

input = 1:5160 # Number of hourly intervals in seven months
multiple_of_8 = (input %% 24) == 0 # Get the lead-time hour for each day
leadtimefordaily <- input[multiple_of_8]# Final lead-times to extract for daily 

# Parameters for CDS request - this shouldn't need changing
request <- list("dataset_short_name" = "seasonal-original-single-levels",
                "originating_centre"= "ecmwf",
                "system" = "5", 
                "variable"       = "total_precipitation" , 
                "year"           = forecast_year,
                "month"          = forecast_month,
                "day"            = "01",
                "leadtime_hour" = leadtimefordaily, 
                "format"         = "netcdf",
                "target"         = paste0("total_precipitation" , "_",
                                          forecast_year, "_",
                                          forecast_month, ".nc"), 
                "area" = forecast_extent)
    
setwd(paste0(directory, "/Data/"))

# This sends the request to CDS. 
ncfile <- wf_request(user = user_number,
                     request = request,
                     transfer = TRUE,
                     time_out = 100000,
                     path = paste0(directory, "/Data/"),
                     verbose = T)

#----------------------------
# Extract temperature 
#----------------------------

input = 1:5160
multiple_of_8 = (input %% 6) == 0
leadtimefor6hour<-input[multiple_of_8]

# Parameters for CDS request - this shouldn't need changing
request <- list("dataset_short_name" = "seasonal-original-single-levels",
                "originating_centre" = "ecmwf",
                "system" = "5",
                "variable"       = "2m_temperature",
                "year"           = forecast_year,
                "month"          = forecast_month,
                "day"            = "01",
                "leadtime_hour" = leadtimefor6hour, 
                "format"         = "netcdf",
                "target"         = paste0("2m_temperature", "_", 
                                          forecast_year, "_",
                                          forecast_month, ".nc"),
                "area" =          forecast_extent)

setwd(paste0(directory, "/Data/"))

# This sends the request to CDS. 
ncfile <- wf_request(user = user_number,
                     request = request,   
                     transfer = TRUE, 
                     time_out= 100000,
                     path = paste0(directory, "/Data/"),
                     verbose = T)


#------------------------------------------------------------------------------
# Step 3: Process SEAS5 Seasonal Forecast data into daily ensemble-median "tif"
#------------------------------------------------------------------------------
library(ncdf4)
library(terra)
library(lubridate)
library(raster)

setwd(directory)

# File name for seven-month temperature SEAS5 data
temperature_forecast <- paste0("Data/2m_temperature_",
                               forecast_year, "_",
                               forecast_month, ".nc")

# File name for seven-month precipitation SEAS5 data
precipitation_forecast <- paste0("Data/total_precipitation_", 
                                 forecast_year, "_", 
                                 forecast_month, ".nc")

# Create directory to store daily data in 
dir.create(paste0(directory,"/Data/daily/"))


#---------------------
# Temperature
#---------------------

extract_daily_temperature(filepath = paste0(directory, "/", temperature_forecast), 
                          forecast_start = forecast_intitiation, 
                          save_dir = paste0(directory, "/Data/daily/"))

#---------------------
# Precipitation
#---------------------

extract_daily_precipitation(filepath = paste0(directory,"/",precipitation_forecast),
                                forecast_start = forecast_intitiation,
                                save_dir = paste0(directory,"/Data/daily/"))


#------------------------------------------------------------------------------
# Step 4: Download historical ERA5 for available dates
#------------------------------------------------------------------------------
# For this data download, you will need Google Earth Engine and Google Drive
# log-ins. Please register for these. 
# Then install `rgee` package and run `ee_install()` and `ee_Authenticate()`


GEE_email <- "" # Set the email address registered with your GEE account

library(dynamicSDM)
library(googledrive)
library(rgee)
library(stars)

rgee::ee_Authenticate(user = GEE_email, drive=T)
rgee::ee_check()
rgee::ee_check_credentials()
rgee::ee_Initialize()

library(googledrive)
googledrive::drive_auth()
googledrive::drive_user()

# Get the dates for historical data extraction 
# 52-week variables so at least 12 months before initiation
historical_start <- forecast_intitiation %m-% months(12)

# However, 3-month lag in data release so only to three months before initiation
historical_end <- forecast_intitiation %m-% months(3)

# Get daily dates for historical data extraction
dates <- dynamicSDM::dynamic_proj_dates(as.character(historical_start),
                                        as.character(historical_end),
                                        interval = 1,
                                        interval.level = "day")

#---------------------
# Temperature
#---------------------

# See Google Earth Engine catalog for more dataset details
# https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY

dynamicSDM::extract_dynamic_raster(dates = dates,
                                   spatial.ext = spatial_extent,
                                   datasetname = "ECMWF/ERA5_LAND/HOURLY",
                                   bandname = "temperature_2m",
                                   spatial.res.metres = 111320,
                                   GEE.math.fun = "mean",
                                   user.email = GEE_email,
                                   varname = "ERA5_historical_temperature",
                                   temporal.res = 1 ,
                                   temporal.direction = "post",
                                   save.directory = paste0(directory, "/Data/daily/"),
                                   resume = T) # Resume continues download from
                                               # previous point if cut off. 


#---------------------
# Precipitation
#---------------------
dynamicSDM::extract_dynamic_raster(dates = dates,
                                   spatial.ext = spatial_extent,
                                   datasetname = "ECMWF/ERA5_LAND/HOURLY",
                                   bandname = "total_precipitation_hourly",
                                   spatial.res.metres = 111320,
                                   GEE.math.fun = "sum",
                                   user.email = GEE_email,
                                   varname = "ERA5_historical_precipitation",
                                   temporal.res = 1 ,
                                   temporal.direction = "post",
                                   save.directory = paste0(directory, "/Data/daily/"),
                                   resume = T)

#------------------------------------------------------------------------------
# Step 5: Download daily from three months previous to in-fill gap
#------------------------------------------------------------------------------

# Now extract the data for the three month gap from previous three months before 
# forecast initiation. Using each month at one-month lead-time. 

# Iterates through each month, extracting SEAS5 data for this interval

for (month in 1:3){
  
  month_to_extract <- forecast_intitiation %m-% months(month)
  
  first_year <- lubridate::year(month_to_extract)
  
  first_month <- lubridate::month(month_to_extract)

#----------------------------
# Extract total precipitation
#----------------------------

input = 1:5160 # Number of hourly intervals in seven months
multiple_of_8 = (input %% 24) == 0 # Get the lead-time hour for each day
leadtimefordaily <- input[multiple_of_8]# Lead-times to extract for daily data

request <- list("dataset_short_name" = "seasonal-original-single-levels", 
                "originating_centre" = "ecmwf", 
                "system" = "5", 
                "variable"       = "total_precipitation" ,
                "year"           = first_year,
                "month"          = first_month,
                "day"            = "01",
                "leadtime_hour" = leadtimefordaily, 
                "format"         = "netcdf",
                "target"         = paste0("Data/total_precipitation", "_",
                                          first_year, "_",
                                          first_month, ".nc"), 
                "area" = forecast_extent)

setwd(directory)


ncfile <- wf_request(
  user = user_number,
  request = request,
  transfer = TRUE,
  time_out = 100000,
  path = directory,
  verbose = T
)


#----------------------------
# Extract temperature 
#----------------------------


input = 1:5160
multiple_of_8 = (input %% 6) == 0
leadtimefor6hour <- input[multiple_of_8]

request <- list("dataset_short_name" = "seasonal-original-single-levels",
                "originating_centre" = "ecmwf",
                "system" = "5", 
                "variable"       = "2m_temperature",
                "year"           = first_year,
                "month"          = first_month,
                "day"            = "01",
                "leadtime_hour" = leadtimefor6hour,
                "format"         = "netcdf",
                "target"         = paste0("Data/2m_temperature", "_",
                                          first_year, "_",
                                          first_month, ".nc"), 
                "area" =          forecast_extent)

setwd(directory)


ncfile <- wf_request(
  user = user_number,
  request = request,
  transfer = TRUE,
  time_out = 100000,
  path = directory,
  verbose = T
)

}

#------------------------------------------------------------------------------
# Step 6: Extract daily ensemble-median for three months previous to in-fill gap
#------------------------------------------------------------------------------

# Again iterates through one to three months from initiation and extracts daily

for (month in 1:3) {
  
  month_to_extract <- forecast_intitiation %m-% months(month)
  first_year <- lubridate::year(month_to_extract)
  first_month <- lubridate::month(month_to_extract)
  
  month_to_end <- month_to_extract %m+% months(1)
  
#---------------------
# Temperature
#---------------------

  extract_daily_temperature(one_month = TRUE, 
                            filepath = paste0(directory, "/Data/", 
                                              "2m_temperature_",
                                              first_year, "_",
                                              first_month, ".nc"), 
                            forecast_start = month_to_extract, 
                            save_dir = paste0(directory, "/Data/daily/"))

#---------------------
# Precipitation
#---------------------
extract_daily_precipitation(one_month = TRUE,
                            filepath = paste0(directory, "/Data/",
                                              "total_precipitation_",
                                              first_year, "_",
                                              first_month, ".nc"), 
                            forecast_start = month_to_extract, 
                            save_dir = paste0(directory, "/Data/daily/"))

}

#------------------------------------------------------------------------------
# Step 7: Combine historical + forecast to calculate 8- and 52-week variables
#------------------------------------------------------------------------------

# Run function to resample all weather data to spatial resolution and extent specified
# This ensures that they will all stack for calculation

resample_daily_weather(daily_dir = paste0(directory, "/Data/daily/"),
                       spatial_resolution = spatial_resolution,
                       spatial_extent = spatial_extent)

# Create directory to store processed weather variables too
dir.create(paste0(directory, "/Data/processed/"))

# 52-week temperature mean and sd
extract_proc_var(
  forecast_intervals = forecast_intervals,
  save_dir = paste0(directory, "/Data/processed/"),
  daily_dir = paste0(directory, "/Data/daily/"),
  period = "annual",
  variable = "temperature"
)

# 8-week temperature mean and sd
extract_proc_var(
  forecast_intervals = forecast_intervals,
  save_dir = paste0(directory, "/Data/processed/"),
  daily_dir = paste0(directory, "/Data/daily/"),
  period = "eight",
  variable = "temperature"
)

# 52-week precipitation sum
extract_proc_var(
  forecast_intervals = forecast_intervals,
  save_dir = paste0(directory, "/Data/processed/"),
  daily_dir = paste0(directory, "/Data/daily/"),
  period = "annual",
  variable = "precipitation"
)

# 8-week precipitation sum
extract_proc_var(
  forecast_intervals = forecast_intervals,
  save_dir = paste0(directory, "/Data/processed/"),
  daily_dir = paste0(directory, "/Data/daily/"),
  period = "eight",
  variable = "precipitation"
)


#------------------------------------------------------------------------------
# Step 8: Download and process most recently available MODIS Land Cover dataset
#------------------------------------------------------------------------------

# Generated moving window matrix for resource variable sum across 10km dispersal
# radius of quelea, accounting for area of cell. 
matrix <- dynamicSDM::get_moving_window(radial.distance = 10000,
                                        spatial.res.degrees = spatial_resolution,
                                        spatial.ext = southernafrica)

# MODIS Land Cover typically released at two-year lag but check catalog:
#https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1
# Change number in brackets to 1 if last year's has been released etc.
resource_date <- forecast_intitiation %m-% years(2) 

# Extracting from Google Earth Engine again - check still authenticated and 
# initialise session. 
library(rgee)
library(stars)
library(googledrive)

rgee::ee_Authenticate(user = GEE_email, drive=T)
rgee::ee_check()
rgee::ee_check_credentials()

googledrive::drive_auth()
googledrive::drive_user()

# Total available shrubland, grassland and cereal cropland cells in area
# Extracted at 500m, then aggregated by 12 to 0.05
# Then using moving window, sum of all cells of those categories in surrounding
# area. See manuscript for more details. 

extract_buffered_raster(dates = resource_date,
                        spatial.ext = spatial_extent,
                        datasetname = "MODIS/061/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = GEE_email,
                        categories = c(5,6, 7),
                        agg.factor = 12,
                        temporal.level = "year",
                        varname = "habitat_availability",
                        save.directory =  paste0(directory, "/Data/processed/"))


# Total available surface water
extract_buffered_raster(dates = resource_date,
                        spatial.ext = spatial_extent,
                        datasetname = "MODIS/061/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = GEE_email,
                        categories = 0,
                        temporal.level = "year",
                        agg.factor = 12,
                        varname = "water_availability",
                        save.directory =  paste0(directory, "/Data/processed/"))


# Total available trees
extract_buffered_raster(dates = resource_date,
                        spatial.ext = spatial_extent,
                        datasetname = "MODIS/061/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = GEE_email,
                        temporal.level = "year",
                        categories = 4,
                        agg.factor = 12,
                        varname = "tree_availability",
                        save.directory =  paste0(directory, "/Data/processed/"))

# Land cover type cells for extraction of seed abundance (next section)
extract_dynamic_raster(dates = resource_date,
                       spatial.ext = spatial_extent,
                       datasetname = "MODIS/061/MCD12Q1",
                       bandname = "LC_Type5",
                       spatial.res.metres = 500,
                       GEE.math.fun = "last", 
                       temporal.direction = "prior",
                       temporal.res = 365,
                       user.email = GEE_email,
                       varname = "landcover",
                       save.directory = paste0(directory, "/Data/processed/"))

#------------------------------------------------------------------------------
# Step 9: Extract EVI data from Google Earth Engine 
#------------------------------------------------------------------------------

# EVI available two-week lag from real-time. 
# For seed abundance modelling need data from 52-weeks to 2-weeks prior to
# forecast initation. 

end_date <- as.character(forecast_intitiation - 14)
end_date <- as.Date(paste0(strsplit(end_date, "-")[[1]][1], "-",
                           strsplit(end_date, "-")[[1]][2], "-14"))
start_date <- end_date - 365

# Gets the dates that the MODIS 16-day EVI is available for
numbers <- 1

for (x in 1:floor(365 / 16)) {
  numbers <- c(numbers, 1 + x * 16)
}

modis_16_dates <- as.Date(paste0(lubridate::year(start_date), "-01-01")) + numbers - 1

numbers <- 1

for (x in 1:floor(365 / 16)) {
  numbers <- c(numbers, 1 + x * 16)
}

modis_16_dates <- c(modis_16_dates,
                           as.Date(paste0(lubridate::year(end_date), "-01-01")) + numbers - 1)

selected_dates <- modis_16_dates[modis_16_dates < end_date]

selected_dates <- selected_dates[(length(selected_dates) - 26):length(selected_dates)]


# Create directory for evi extraction
dir.create(paste0(directory, "/Data/evi/"))

# Extract rasters at 16-day intervals for deriving EVI characterstics from
# See catalogue for dataset details:
# https://developers.google.com/earth-engine/datasets/catalog/MODIS_MCD43A4_006_EVI

dynamicSDM::extract_dynamic_raster(dates = selected_dates,
                                   spatial.ext = spatial_extent,
                                   datasetname = "MODIS/MCD43A4_006_EVI", 
                                   bandname = "EVI",
                                   spatial.res.metres = 500,
                                   GEE.math.fun = "mean",
                                   resume = T,
                                   user.email = GEE_email,
                                   varname = "MODIS_historical_EVI",
                                   temporal.res = 16 ,
                                   temporal.direction = "prior",
                                   save.directory = paste0(directory,"/Data/evi/"))

#------------------------------------------------------------------------------
# Step 10: Process, classify and project stages forwards to each interval
#------------------------------------------------------------------------------
library(randomForest)

# Read in MODIS Land Cover Yearly for most recent year extracted in previous step
landcover <- terra::rast(list.files(paste0(directory, "/Data/processed/"),
                                    full.names = T,
                                    pattern = "landcover")[1])

# Run custom function for generating EVI characterstics from 16-day data extracted
# in previous step.
seed_forecast_dataframe <- get_evi_characterstics(directory = paste0(directory, "/Data/evi/"),
                                                  land_cover = landcover,
                                                  spatial_ext = spatial_extent)


# Read in fitted Random Forest classification models for classifying vegetation
# phenology stages based upon EVI characterstics. 

grass_classification_model <- readRDS(paste0(directory,"/Data/models/classification_model_grass.rds"))

cereal_classification_model <- readRDS(paste0(directory,"/Data/models/classification_model_cereal.rds"))


# Read in extracted mean vegetation growth stage lengths, derived from MODIS Land Cover Dynamics Yearly
# See: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1

mean_lengths<-list.files(paste0(directory, "/Data/average_length_phenology/"),
                         pattern = ".tif",
                         full.names = T)

mean_lengths<-terra::rast(list(terra::rast(mean_lengths[1]), # green-up
                              (terra::rast(mean_lengths[4]) + terra::rast(mean_lengths[2])), # mid-greenup + maturity as classes collapsed
                              (terra::rast(mean_lengths[5]) + terra::rast(mean_lengths[6])), # peak + senesence as classes collapsed
                               terra::rast(mean_lengths[3]))) # mid-green-down

# Aggregate by mean to get general length at coarser resolution
mean_lengths <- terra::aggregate(mean_lengths, 12, fun = "mean", na.rm = T)


# Use custom function to project seed availability at every forecast intervals
# based upon classifying EVI characterstics

# Cereal seed
project_seed_availability(EVI_data_frame = seed_forecast_dataframe,
                          mean_lengths = mean_lengths ,
                          forecast_intitiation = forecast_intitiation,
                          type = "cereal",
                          moving.window.matrix = matrix,
                          agg.factor = 12,
                          model = cereal_classification_model,
                          save_dir = paste0(directory, "/Data/processed/"), 
                          forecast_intervals = forecast_intervals)

# Grass seed
project_seed_availability(EVI_data_frame = seed_forecast_dataframe,
                          mean_lengths = mean_lengths ,
                          forecast_intitiation = forecast_intitiation,
                          type = "grass",
                          moving.window.matrix = matrix,
                          agg.factor = 12,
                          model = grass_classification_model, 
                          save_dir = paste0(directory, "/Data/processed/"), 
                          forecast_intervals = forecast_intervals)


# We want to add cereal and grass seed availability for each data to get total

# Get list of seed rasters
all_proc_rasts <- list.files(paste0(directory, "/Data/processed/"), full.names = T)
seed_rasters <- all_proc_rasts[grepl("seed", all_proc_rasts)]

# Iterate through each forecast interval date and combine grass + cereal

for (interval in 1:length(forecast_intervals)) {
  
  seed_raster <- seed_rasters[grepl(forecast_intervals[interval], seed_rasters)]
  seed_raster <- terra::rast(seed_raster[1]) + terra::rast(seed_raster[2])
  writeRaster(seed_raster, file = paste0(directory, "/Data/processed/",
                                     forecast_intervals[interval],
                                     "_seed_abundance.tif"))
  
}

# Remove tifs from temporary directory as may include large file sizes.
file.remove(list.files(tempdir(),pattern=".tif",full.names = T))

#------------------------------------------------------------------------------
# Step 11: Combined data into raster stacks for each forecast interval 
#------------------------------------------------------------------------------

# Create directory for saving covariate stacks to

dir.create(paste0(directory, "/Data/covariates/"))

# Read in the annual resource rasters that remain static for each interval
static <- all_proc_rasts[grepl("tree|habitat|water",all_proc_rasts)]
static <-terra::rast(static)

# For each interval, this function reads in the dynamic weather/seed rasters 
# generated so far; resamples them to same extent/resolution; and then
# combines them with the annual resource variables to save as a raster stack.
dynamicSDM::dynamic_proj_covariates(dates = forecast_intervals,
                                    varnames = c("mean_annual_temperature",
                                                 "mean_eight_temperature",
                                                 "sd_annual_temperature",
                                                 "sd_eight_temperature",
                                                 "sum_annual_precipitation",
                                                 "sum_eight_precipitation",
                                                 "seed_abundance"),
                                    local.directory = paste0(directory, "/Data/processed/") , 
                                    spatial.ext = spatial_extent ,
                                    spatial.mask = spatial_extent,
                                    spatial.res.degrees = 0.05 ,
                                    resample.method = c("bilinear" ,
                                                        "bilinear" ,
                                                        "bilinear" ,
                                                        "bilinear" ,
                                                        "bilinear" ,
                                                        "bilinear" ,
                                                        "near"),
                                    cov.file.type = "tif",
                                    cov.prj =  "+proj=longlat +datum=WGS84",
                                    save.directory = paste0(directory, "/Data/covariates/"), 
                                    static.rasters = static , 
                                    static.varnames = c("tree_availability",
                                                        "habitat_availability",
                                                        "water_availability"), 
                                    static.resample.method = "near")


#------------------------------------------------------------------------------
# Step 12: Species distribution modelling with all available data
#------------------------------------------------------------------------------
library(dplyr)
library(dynamicSDM)
library(readr)
library(gbm)
library(pROC)
# Read in breeding and non-breeding data with associated explanatory variables. 

# See Dobson et al. (2023) for details and code for processing occurrence data
# and extracting dynamic explanatory variables
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14101

quelea_breeding_data <- read.csv("Data/breeding_distribution_quelea.csv")

quelea_nonbreeding_data <- read.csv("Data/nonbreeding_distribution_quelea.csv")

# Extract biome layer for spatial blocking of occurrence data
dynamicSDM::extract_dynamic_raster(dates = "2001-01-01",
                                   spatial.ext = spatial_extent,
                                   datasetname ="OpenLandMap/PNV/PNV_BIOME-TYPE_BIOME00K_C/v01", 
                                   bandname = "biome_type",
                                   spatial.res.metres = 1000 ,
                                   GEE.math.fun = "first",
                                   resume = T,
                                   user.email = GEE_email,
                                   varname = "biome_layer",
                                   temporal.res = 1 ,
                                   temporal.direction = "post",
                                   save.directory = paste0(directory,"/Data/"))

# Read in biome layer
biome_layer <- terra::rast(paste0(directory,"/Data/2001-01-01_biome_layer.tif"))

# Spatial and temporal blocking of occurrence data for model fitting
breed_blocked<-spatiotemp_block(occ.data = quelea_breeding_data,
                                vars.to.block.by = colnames(quelea_breeding_data[11:16]), 
                                spatial.layer = biome_layer,
                                spatial.split.degrees = 3,
                                temporal.block = "year",
                                n.blocks = 6,
                                iterations = 5000)

nonbreed_blocked <- spatiotemp_block(occ.data = quelea_nonbreeding_data,
                                  vars.to.block.by =  colnames(quelea_nonbreeding_data[11:16]),
                                  spatial.layer = biome_layer,
                                  spatial.split.degrees = 3,
                                  temporal.block = "year",
                                  n.blocks = 6,
                                  iterations = 5000)

model_variables <- c("mean_annual_temperature",
                     "sd_annual_temperature",
                     "mean_eight_temperature",
                     "sd_eight_temperature",
                     "sum_annual_precipitation",
                     "sum_eight_precipitation",
                     "seed_abundance",
                     "water_availability",
                     "habitat_availability",
                     "tree_availability")

# Fit Boosted Regression Tree models to each data block, with weights 
# representative of avian sampling effort
breed_SDM_models <- brt_fit(occ.data = breed_blocked,
                            response.col = "presence_absence",
                            varnames = model_variables,
                            distribution = "bernoulli",
                            block.col = "BLOCK.CATS",
                            weights.col = "weights")

nonbreed_SDM_models <- brt_fit(occ.data = nonbreed_blocked,
                            response.col = "presence_absence",
                            varnames = model_variables,
                            distribution = "bernoulli",
                            block.col = "BLOCK.CATS",
                            weights.col = "weights")

# Measure AUC of breeding and non-breeding BRT ensemble.
# These will be used for performance-weighted projections of quelea suitability

breed_sdm_weights = NULL

for (x in 1:length(unique(breed_blocked$BLOCK.CATS))) {
  block <- x
  mod <- breed_SDM_models[[x]]
  test <- breed_blocked[breed_blocked$BLOCK.CATS == block,]
  
  test$pred <- predict(mod, newdata = test, type = 'response')
  auc <- pROC::auc(pROC::roc(presence_absence ~ pred, data = test, quiet = T))
  breed_sdm_weights <- rbind(breed_sdm_weights, auc)
  
}


nonbreed_sdm_weights = NULL

for (x in 1:length(unique(nonbreed_blocked$BLOCK.CATS))) {
  block <- unique(nonbreed_blocked$BLOCK.CATS)[x]
  mod <- nonbreed_SDM_models[[x]]
  test <- nonbreed_blocked[nonbreed_blocked$BLOCK.CATS == block, ]
  
  test$pred <- predict(mod, newdata = test, type = 'response')
  auc <- pROC::auc(pROC::roc(presence_absence ~ pred, data = test, quiet = T))
  nonbreed_sdm_weights <- rbind(nonbreed_sdm_weights, auc)
  
}


#------------------------------------------------------------------------------
# Step 13: Project distribution models onto forecast variables for each interval
#------------------------------------------------------------------------------

# Get forecast intervals for projecting breeding D-SDMs onto. 
breeding_intervals <-  forecast_intervals[lubridate::month(forecast_intervals) %in%  c(12,1,2,3,4,5)]

# Get forecast intervals for projecting non-breeding D-SDMs onto. 
nonbreeding_intervals <- forecast_intervals[!lubridate::month(forecast_intervals) %in%  c(12,1,2,3,4,5)]

# Create folders to save forecast output to
dir.create(paste0(directory,"/Output/"))
dir.create(paste0(directory,"/Output/projections"))

# Functions generate performance weighted forecasts of distribution suitability
# for quelea across southern Africa. 
dynamic_proj(dates = breeding_intervals,
             projection.method = c("proportional"),
             local.directory = paste0(directory,"/Data/covariates/"),
             sdm.mod = breed_SDM_models,
             sdm.weight = as.numeric(breed_sdm_weights),
             cov.file.type = "tif",
             spatial.mask = spatial_extent,
             proj.prj = "+proj=longlat +datum=WGS84",
             save.directory = paste0(directory,"/Output/projections"))

dynamic_proj(dates = nonbreeding_intervals,
             projection.method = c("proportional"),
             local.directory = paste0(directory,"/Data/covariates/"),
             sdm.mod = nonbreed_SDM_models,
             sdm.weight = as.numeric(nonbreed_sdm_weights),
             cov.file.type = "tif",
             spatial.mask = spatial_extent,
             proj.prj = "+proj=longlat +datum=WGS84",
             save.directory = paste0(directory,"/Output/projections"))


#------------------------------------------------------------------------------
# Step 14: Plot seasonal forecasts for sharable images
#------------------------------------------------------------------------------

library(png)
library(ggplot2)

# Create directory to save Seasonal Forecast images to. 
dir.create(paste0(directory,"/Output/images"))

# Get list of projection files, iterate through them, plotting and saving png image
list_of_projections <- list.files(paste0(directory,"/Output/projections"), full.names = T)

for(proj in 1:length(list_of_projections)){

    file_name<-list_of_projections[proj]
    
    # Get date of projection
    date<-strsplit(file_name,"/")[[1]]
    date<-strsplit(date[length(date)],"_")[[1]][1]
    
    # Get lead-time of projection
    leadtime<- as.numeric(forecast_intervals[proj]-forecast_intervals[1])
  
    # Read in tif file and create data frame for plotting
    proj_rast <- terra::rast(paste0(file_name))
    proj_rast <- terra::crop(proj_rast, spatial_extent)
    proj_rast <- terra::mask(proj_rast, spatial_extent)
    proj_rast <- as.data.frame(proj_rast, xy = T)
    colnames(proj_rast) <- c("x", "y", "var")
    
    # Define breaks for suitability categories
    breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
   
    cols <- c("#F0F0F0",
              "#40863A",
              "#9EBD49",
              "#FBF357",
              "#F4C12F",
              "#ED8E07",
              "#DD4704",
              "#cc0000",
              "#990000",
              "#660000")
    
    # Create suitability categories using cut()
    proj_rast$categories <- cut(proj_rast$var, breaks = breaks, labels = FALSE)
    proj_rast$categories <- factor(proj_rast$categories, levels = 1:10)
    
    p1<-ggplot(data = ra) +
      geom_raster(aes(x = x, y = y, fill = categories))+ 
      scale_fill_manual(values = cols,
                        drop = F,
                        name = "Forecast\nsuitability",
                        labels =c("< 0.1",
                                  "0.1 - 0.2",
                                  "0.2 - 0.3",
                                  "0.3 - 0.4",
                                  "0.4 - 0.5",
                                  "0.5 - 0.6",
                                  "0.6 - 0.7",
                                  "0.7 - 0.8",
                                  "0.8 - 0.9",
                                  "0.9 <") )+
      theme_bw()+
      coord_sf(expand = FALSE) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.key.size = unit(1, "cm"),
            legend.text = element_text(size=25),
            legend.title = element_text(size=25,face="bold"),
            plot.title = element_text(size=25,hjust = 0.5,face="bold"),
            axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_text(size=15))+
      geom_sf(data = spatial_extent, fill = "NA",color = "black", size = 15, lwd = 1)
    
    
   # Save plot image 
    ggsave(p1,file =paste0(directory,"/Output/images/",date,"_",leadtime,"_proportional.png"),width=8, height=6)
    
  
}
###############################################################################
# If you have any questions or need help running any of the above code to
# generate seasonal forecasts of red-billed quelea distribution, please contact
# Rachel Dobson at eerdo@leeds.ac.uk
###############################################################################