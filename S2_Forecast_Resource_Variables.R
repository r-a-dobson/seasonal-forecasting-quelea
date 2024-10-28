#-------------------------------------------------------------------------------
# (2) Code for processing seasonal forecasts of resource variable
#------------------------------------------------------------------------------

# Please ensure that you have downloaded the GitHub repository that contains the
# associated data sets and custom functions.
# Found at: https://github.com/r-a-dobson/seasonal-forecasting-quelea

#------------------------------------------------------------------------------
# Step 1: Set up and initiate Google Earth Engine and Google Drive
#------------------------------------------------------------------------------

GEE_email <- "" # Set the email address registered with your GEE account

# Provide path to "seasonal_forecasting_quelea" directory downloaded from GitHub
directory <- "C:/Users/XXXXX/Downloads/seasonal_forecasting_quelea/"

# Set this as your working directory
setwd(directory)

# Read in the custom functions for near-term hindcasting
source(paste0(directory, "/", "Functions_For_Forecast.R"))

# Load required packages
library(dynamicSDM)
library(googledrive)
library(rgee)
library(stars)
library(lubridate)
library(googledrive)

# Authenticate and initialise Google Earth Engine
rgee::ee_Authenticate(user = GEE_email, drive=T)
rgee::ee_check()
rgee::ee_check_credentials()
rgee::ee_Initialize()

googledrive::drive_auth()
googledrive::drive_user()

# Create data frame to iterate through each monthly near-term hindcast date

years <- c(2004:2016)

months <- c(1:12)

days <- c(1)

dataframe <- as.data.frame(expand.grid(years, months, days))

for (i in 1:nrow(dataframe)){

  yr <- dataframe$Var1[i]
  mn <- dataframe$Var2[i]
  day <- dataframe$Var3[i]
  
  forecast_intitiation <- as.Date(paste0(yr, "-", mn, "-", day))
  
  historical_start <- forecast_intitiation %m-% months(12)
  historical_end <- forecast_intitiation %m-% months(3)

  forecast_intervals <- dynamic_proj_dates(startdate = as.character(forecast_intitiation),
                                           enddate =    as.character(forecast_intitiation %m+% months(7)), 
                                           interval = 2, # Represents two of each interval i.e. fortnightly
                                           interval.level = "week") # Change to day or month if desired

# Get daily dates for historical data extraction
dates <- dynamicSDM::dynamic_proj_dates(as.character(historical_start),
                                        as.character(historical_end),
                                        interval = 1,
                                        interval.level = "day")

# Create directory to store near-term hindcasts of resource variables
save_directory <- paste0(directory, "/Data/resource_hindcasts/",forecast_intitiation)

dir.create(save_directory, recursive = T)

#------------------------------------------------------------------------------
# Step 2: Download most recently available MODIS Land Cover dataset
#------------------------------------------------------------------------------

# Generate moving window matrix for summing resource variables across the 10km
# dispersal radius of quelea.
matrix <- dynamicSDM::get_moving_window(radial.distance = 10000,
                                        spatial.res.degrees = spatial_resolution,
                                        spatial.ext = southernafrica)


# MODIS Land Cover Type released at two-year lag from real-time, replicate this
# for near-term hindcasting
resource_date <- forecast_intitiation %m-% years(2) 

# Total available "shrubland", "grassland" and "cereal cropland" cells in area
# Extracted at 0.04 degree spatial resolution, then aggregated by 12 to 0.05
# degree spatial resolution. Then, using the moving window, sum values of each
# cell with surrounding eight cells to get final resource availability value.

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
                        save.directory =  save_directory)

# Total available surface water in accessible area
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
                        save.directory = save_directory)

# Total available trees in accessible area
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
                        save.directory =  save_directory)

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
                       save.directory = save_directory)

#------------------------------------------------------------------------------
# Step 3: Extract EVI data from Google Earth Engine 
#------------------------------------------------------------------------------
# Enhanced vegetation index (EVI) data are available at two-week lag from
# real-time. For seed abundance modelling need data from 52-weeks to 2-weeks
# prior to forecast initiation.

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


# Create directory for EVI extraction
dir.create(paste0(save_directory, "/evi/"))

# Extract EVI rasters at 16-day intervals for deriving EVI characteristics from
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
                                   save.directory = paste0(save_directory, "/evi/"))

#------------------------------------------------------------------------------
# Step 4: Process, classify and project stages forwards to each interval
#------------------------------------------------------------------------------
library(randomForest)

# Read in MODIS Land Cover Type Yearly for the most recent year 
landcover <- terra::rast(list.files(save_directory,
                                    full.names = T,
                                    pattern = "landcover")[1])

# Run custom function for extracting EVI characteristics from 16-day EVI data
# extracted in previous step.
seed_forecast_dataframe <- get_evi_characteristics(directory = paste0(save_directory, "/evi/"),
                                                  land_cover = landcover,
                                                  spatial_ext = spatial_extent)

# Read in fitted Random Forest classification models for classifying vegetation
# phenology stages based upon 16-day EVI characteristics. 
grass_classification_model <- readRDS(paste0(directory,"/Data/models/classification_model_grass.rds"))

cereal_classification_model <- readRDS(paste0(directory,"/Data/models/classification_model_cereal.rds"))


# Read in extracted mean vegetation growth stage lengths, derived from MODIS
# Land Cover Dynamics Yearly See:
# https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1
mean_lengths<-list.files(paste0(directory, "/Data/average_length_phenology/"),
                         pattern = ".tif",
                         full.names = T)

mean_lengths<-terra::rast(list(terra::rast(mean_lengths[1]), # green-up
                               (terra::rast(mean_lengths[4]) + terra::rast(mean_lengths[2])), # mid-greenup + maturity as classes collapsed
                               (terra::rast(mean_lengths[5]) + terra::rast(mean_lengths[6])), # peak + senesence as classes collapsed
                               terra::rast(mean_lengths[3]))) # mid-green-down

# Aggregate by mean to get average durations at coarser resolutions
mean_lengths <- terra::aggregate(mean_lengths, 12, fun = "mean", na.rm = T)


# Use custom function to project seed availability at every forecast interval
# using 16-day EVI characteristics

# Cereal seed
project_seed_availability(EVI_data_frame = seed_forecast_dataframe,
                          mean_lengths = mean_lengths ,
                          forecast_intitiation = forecast_intitiation,
                          type = "cereal",
                          moving.window.matrix = matrix,
                          agg.factor = 12,
                          model = cereal_classification_model,
                          save_dir = save_directory, 
                          forecast_intervals = forecast_intervals)

# Grass seed
project_seed_availability(EVI_data_frame = seed_forecast_dataframe,
                          mean_lengths = mean_lengths ,
                          forecast_intitiation = forecast_intitiation,
                          type = "grass",
                          moving.window.matrix = matrix,
                          agg.factor = 12,
                          model = grass_classification_model, 
                          save_dir = save_directory, 
                          forecast_intervals = forecast_intervals)

# We want to add cereal and grass seed availability for each date to get total
# seed available to quelea 

# Get list of generated seed rasters
all_proc_rasts <- list.files(save_directory, full.names = T)
seed_rasters <- all_proc_rasts[grepl("seed", all_proc_rasts)]

# Iterate through each forecast interval date and combine grass + cereal
for (interval in 1:length(forecast_intervals)) {
  
  seed_raster <- seed_rasters[grepl(forecast_intervals[interval], seed_rasters)]
  
  seed_raster <- terra::rast(seed_raster[1]) + terra::rast(seed_raster[2])
  
  writeRaster(seed_raster, file = paste0(save_directory,
                                         forecast_intervals[interval],
                                         "_seed_abundance.tif"))
}

# Remove files from temporary directory as may contain very large files
file.remove(list.files(tempdir(), pattern = ".tif", full.names = T))

}
