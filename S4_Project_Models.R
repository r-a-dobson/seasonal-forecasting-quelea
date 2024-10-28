#-------------------------------------------------------------------------------
# (4) Code for projecting dynamic species distribution models to hindcast conditions
#------------------------------------------------------------------------------

# Please ensure that you have downloaded the GitHub repository that contains the
# associated data sets and custom functions written for the following R scripts.

library(dynamicSDM)

# Provide path to "seasonal_forecasting_quelea" directory downloaded from GitHub
directory <- "C:/Users/XXXXX/Downloads/seasonal_forecasting_quelea/"

# Set this directory as your working directory for analyses
setwd(directory)

years<-c(2004:2016)

months<-c(1:12)

days<-c(1)

dataframe<-as.data.frame(expand.grid(years,months,days))

model_directory <- paste0(directory,"SDMs/")

results <- NULL

for (x in 1:nrow(dataframe)){
  
  yr <- dataframe$Var1[i]
  mn <- dataframe$Var2[i]
  day <- dataframe$Var3[i]
  
  forecast_intitiation <- as.Date(paste0(yr,"-",mn,"-",day)) 
  
  end_date <- forecast_intitiation %m+% months(7)
  
  forecast_intervals <- dynamic_proj_dates(startdate = as.character(forecast_intitiation),
                                           enddate =    as.character(forecast_intitiation %m+% months(7)),  
                                           interval = 2, # Represents two of each interval i.e. fortnightly
                                           interval.level = "week") # Change to day or month if desired
  
  # Create directory for saving covariate stacks to
  covariate_directory <- paste0(directory, "/Data/",forecast_intitiation,"/covariates/")
  dir.create(covariate_directory, recursive = T)
  
  #-------------------------------------
  # Static resource for each interval
  #------------------------------------
  resource_dir <- paste0(directory, "/Data/resource_hindcasts/",forecast_intitiation)
  resource_files <- list.files(resource_dir, full.names = T)
  static <- resource_files[grepl("tree|habitat|water", resource_files)]
  static <- terra::rast(static)
  
  #--------------------------------------------
  # Combine near-term hindcasts for each interval
  #------------------------------------------
  
  weather_variable_dir <- paste0(working_dir, forecast_intitiation, "/ensemble_medians/")

  dynamicSDM::dynamic_proj_covariates(dates = forecast_intervals,
                                      varnames = c("annual_mean_temperature",
                                                   "eight_mean_temperature",
                                                   "annual_sd_temperature",
                                                   "eight_sd_temperature",
                                                   "annual_sum_precipitation",
                                                   "eight_sum_precipitation",
                                                   "seed_abundance"),
                                      local.directory = c(resource_dir, weather_variable_dir), 
                                      spatial.ext = spatial_extent ,
                                      spatial.mask = spatial_extent,
                                      spatial.res.degrees = spatial_resolution ,
                                      resample.method = c("bilinear" ,
                                                          "bilinear" ,
                                                          "bilinear" ,
                                                          "bilinear" ,
                                                          "bilinear" ,
                                                          "bilinear" ,
                                                          "near"),
                                      cov.file.type = "tif",
                                      cov.prj =  "+proj=longlat +datum=WGS84",
                                      save.directory = covariate_directory, 
                                      static.rasters = static , 
                                      static.varnames = c("tree_availability",
                                                          "habitat_availability",
                                                          "water_availability"), 
                                      static.resample.method = "near")
  
# Get forecast intervals for projecting breeding D-SDMs onto. 
breeding_intervals <-  forecast_intervals[lubridate::month(forecast_intervals) %in%  c(12,1,2,3,4,5)]

# Get forecast intervals for projecting non-breeding D-SDMs onto. 
nonbreeding_intervals <- forecast_intervals[!lubridate::month(forecast_intervals) %in%  c(12,1,2,3,4,5)]


# Create directories for saving projections to
proj_dir <- paste0(directory, "/Output/projections/", forecast_intitiation)
proj_dir_MAX <- paste0(proj_dir, "/maximal")
proj_dir_REAL <- paste0(proj_dir, "/realtime")
dir.create(proj_dir_REAL, recursive = T)
dir.create(proj_dir_MAX, recursive = T)


# Read in real-time/maximal D-SDMs for current hindcast interval 

model_dir <- list.files(paste0(directory,"SDMs/"), full.names = T, pattern = forecast_intitiation)

NBR_MAX <- readRDS(model_dir[grepl("NBR_MAX_.rds", model_dir)])

BR_MAX <- readRDS(model_dir[grepl("_BR_MAX_.rds", model_dir)])

NBR_REAL <- readRDS(model_dir[grepl("NBR_REAL_.rds", model_dir)])

BR_REAL <- readRDS(model_dir[grepl("_NBR_MAX_.rds", model_dir)])

# Project breeding D-SDMs for breeding dates and vice versa

dynamic_proj(dates = breeding_intervals,
             projection.method = c("proportional"),
             local.directory = covariate_directory,
             sdm.mod = BR_REAL,
             cov.file.type = "tif",
             spatial.mask = spatial_extent,
             proj.prj = "+proj=longlat +datum=WGS84",
             save.directory = proj_dir_REAL)

dynamic_proj(dates = breeding_intervals,
             projection.method = c("proportional"),
             local.directory = covariate_directory,
             sdm.mod = BR_MAX,
             cov.file.type = "tif",
             spatial.mask = spatial_extent,
             proj.prj = "+proj=longlat +datum=WGS84",
             save.directory = proj_dir_MAX)

dynamic_proj(dates = nonbreeding_intervals,
             projection.method = c("proportional"),
             local.directory = covariate_directory,
             sdm.mod = NBR_MAX,
             cov.file.type = "tif",
             spatial.mask = spatial_extent,
             proj.prj = "+proj=longlat +datum=WGS84",
             save.directory = proj_dir_MAX)

dynamic_proj(dates = nonbreeding_intervals,
             projection.method = c("proportional"),
             local.directory = covariate_directory,
             sdm.mod = NBR_REAL,
             cov.file.type = "tif",
             spatial.mask = spatial_extent,
             proj.prj = "+proj=longlat +datum=WGS84",
             save.directory = proj_dir_REAL)


}
