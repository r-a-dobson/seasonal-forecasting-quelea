#-------------------------------------------------------------------------------
# (3) Code for generating dynamic species distribution models
#------------------------------------------------------------------------------

# Please ensure that you have downloaded the GitHub repository that contains the
# associated data sets and custom functions.

# Provide path to "seasonal_forecasting_quelea" directory downloaded from GitHub
directory <- "C:/Users/XXXXX/Downloads/seasonal_forecasting_quelea/"

# Set this directory as your working directory for analyses
setwd(directory)

#------------------------------------------------------------------------------
# Step 1: Extract dynamic weather variables for occurrence records
#------------------------------------------------------------------------------

library(dplyr)
library(dynamicSDM)
library(readr)
library(gbm)
library(pROC)

# Read in breeding and non-breeding occurrence data
quelea_br_data <- read.csv("Data/breeding_distribution_quelea.csv")

quelea_nbr_data <- read.csv("Data/nonbreeding_distribution_quelea.csv")

#-------------------------------------------
# Extract 8- and 52-week weather variables
#-------------------------------------------

# First, you need to download CHELSA-W5E5 daily temperature and preciptiation
# data for 1982-2016 from:
# https://chelsa-climate.org/chelsa-w5e5-v1-0-daily-climate-data-at-1km-resolution/

precipitation_dir <- "" # Insert path to CHELSA-W5E5 precipitation files 

temperature_dir <- "" # Insert path to CHELSA-W5E5 temperature files

quelea_br_data <- extract_chelsa(quelea_br_data, precipitation_dir, temperature_dir)

quelea_nbr_data <- extract_chelsa(quelea_nbr_data, precipitation_dir, temperature_dir)


#------------------------------------------------------------------------------
# Step 2: Extract dynamic resource variables for occurrence records
#------------------------------------------------------------------------------

# Create directories to store extracted resource data 
dir.create("dynamic_resource_breed")
dir.create("dynamic_resource_roost")

save_directory_breed_resource <- paste0(directory, "/dynamic_resource_breed")
save_directory_roost_resource <- paste0(directory, "/dynamic_resource_roost")

# Get moving window matrix size
matrix <- dynamicSDM::get_moving_window(radial.distance = 10000,
                                        spatial.res.degrees = 0.05,
                                        spatial.ext = c(-35, -6, 10, 40)  )
#----------------------
# a) Breeding records
#----------------------

# Total available grassland and cereal cropland
extract_buffered_coords(occ.data = quelea_br_data,
                        datasetname = "MODIS/006/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        categories = c(6, 7),
                        agg.factor = 12,
                        varname = "habitat_availability",
                        save.directory = save_directory_breed_resource)

# Median senesence vegetation phenology stage timing
extract_buffered_coords(occ.data = quelea_br_data,
                        datasetname = "MODIS/006/MCD12Q2",
                        bandname = "Senescence_1",
                        spatial.res.metres = 500,
                        GEE.math.fun = "median",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        varname = "senesence",
                        save.directory = save_directory_breed_resource,
                        agg.factor = 12)

# Median dormancy vegetation phenology stage timing
extract_buffered_coords(occ.data = quelea_br_data,
                        datasetname = "MODIS/006/MCD12Q2",
                        bandname = "Dormancy_1",
                        spatial.res.metres = 500,
                        GEE.math.fun = "median",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        varname = "dormancy",
                        save.directory = save_directory_breed_resource,
                        agg.factor = 12)


# Total available surface water
extract_buffered_coords(occ.data = quelea_br_data,
                        datasetname = "MODIS/006/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        categories = 0,
                        agg.factor = 12,
                        varname = "water_availability",
                        save.directory = save_directory_breed_resource)

# Total available trees
extract_buffered_coords(occ.data = quelea_br_data,
                        datasetname = "MODIS/006/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        categories = 4,
                        agg.factor = 12,
                        varname = "tree_availability",
                        save.directory = save_directory_breed_resource)

#----------------------
# b) Roosting records
#----------------------

# Total available grassland and cereal cropland

extract_buffered_coords(occ.data = quelea_nbr_data,
                        datasetname = "MODIS/006/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        categories = c(6, 7),
                        agg.factor = 12,
                        varname = "habitat_availability",
                        save.directory = save_directory_roost_resource)

# Median senesence vegetation phenology stage timing
extract_buffered_coords(occ.data = quelea_nbr_data,
                        datasetname = "MODIS/006/MCD12Q2",
                        bandname = "Senescence_1",
                        spatial.res.metres = 500,
                        GEE.math.fun = "median",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        varname = "senesence",
                        save.directory = save_directory_roost_resource,
                        agg.factor = 12)

# Median dormancy vegetation phenology stage timing
extract_buffered_coords(occ.data = quelea_nbr_data,
                        datasetname = "MODIS/006/MCD12Q2",
                        bandname = "Dormancy_1",
                        spatial.res.metres = 500,
                        GEE.math.fun = "median",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        varname = "dormancy",
                        save.directory = save_directory_roost_resource,
                        agg.factor = 12)


# Total available surface water
extract_buffered_coords(occ.data = quelea_nbr_data,
                        datasetname = "MODIS/006/MCD12Q1",
                        bandname = "LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix = matrix,
                        user.email = user.email,
                        resume = T,
                        save.method = "split",
                        temporal.level = "year",
                        categories = 0,
                        agg.factor = 12,
                        varname = "water_availability",
                        save.directory = save_directory_roost_resource)

# Total available trees
extract_buffered_coords(occ.data=quelea_nbr_data,
                        datasetname = "MODIS/006/MCD12Q1",
                        bandname="LC_Type5",
                        spatial.res.metres = 500,
                        GEE.math.fun = "sum",
                        moving.window.matrix=matrix,
                        user.email=user.email,
                        resume=T,
                        save.method="split",
                        temporal.level="year",
                        categories=4,
                        agg.factor = 12,
                        varname = "tree_availability",
                        save.directory=save_directory_roost_resource)

#------------------------------------------------------------------------------
# Step 3: Combine extracted explanatory variables 
#------------------------------------------------------------------------------
quelea_br_data <- extract_coords_combine(varnames = variablenames,
                                               local.directory = c(save_directory_breed_resource,
                                                                   save_directory_breed),
                                               set_class = TRUE,
                                               col_classes = c(sapply(quelea_br_data,class),
                                                               mean_eight_temperature = "numeric",
                                                               mean_annual_temperature = "numeric",
                                                               sd_eight_temperature = "numeric",
                                                               sd_annual_temperature = "numeric",
                                                               sum_eight_precipitation = "numeric",
                                                               sum_annual_precipitation = "numeric",
                                                               grass_crop_percentage = "integer",
                                                               senescence = "integer",
                                                               dormancy = "integer",
                                                               water_availability = "integer",
                                                               tree_availability = "integer"))

quelea_nbr_data <- extract_coords_combine(varnames = variablenames,
                                                  local.directory = c(save_directory_roost_resource,
                                                                      save_directory_roost),
                                                  set_class = TRUE,
                                                  col_classes = c(sapply(quelea_nbr_data,class),
                                                                  mean_eight_temperature = "numeric",
                                                                  mean_annual_temperature = "numeric",
                                                                  sd_eight_temperature = "numeric",
                                                                  sd_annual_temperature = "numeric",
                                                                  sum_eight_precipitation = "numeric",
                                                                  sum_annual_precipitation = "numeric",
                                                                  grass_crop_percentage = "integer",
                                                                  senescence = "integer",
                                                                  dormancy = "integer",
                                                                  water_availability = "integer",
                                                                  tree_availability = "integer"))

#------------------------------------------------------------------------------
# Step 4: Extract biome layer and split occurrence data into blocks 
#------------------------------------------------------------------------------

variables_to_block_by <- c("mean_annual_temperature",
                           "sd_annual_temperature",
                           "sum_annual_precipitation")

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
breed_blocked<-spatiotemp_block(occ.data = quelea_br_data,
                                vars.to.block.by = variables_to_block_by, 
                                spatial.layer = biome_layer,
                                spatial.split.degrees = 3,
                                temporal.block = "year",
                                n.blocks = 6,
                                iterations = 5000)

nonbreed_blocked <- spatiotemp_block(occ.data = quelea_nbr_data,
                                     vars.to.block.by =  variables_to_block_by,
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

#------------------------------------------------------------------------------
# Step 5: Generate baseline D-SDMs using all available data 
#------------------------------------------------------------------------------
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

#-------------------------------------------------
# Measure AUC to identify top-performing D-SDM
#-------------------------------------------------
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

#-------------------------------------------------------------------------------
# Get parameters of best performing baseline D-SDM to use for hindcasting D-SDMs
#-------------------------------------------------------------------------------

baseline_breeding <- breed_blocked[[which.max(breed_sdm_weights)]]
br_n_trees <- baseline_breeding$n.trees
br_shrinkage <- baseline_breeding$shrinkage
br_interaction <- baseline_breeding$interaction.depth


baseline_nonbreeding <- nonbreed_blocked[[which.max(nonbreed_sdm_weights)]]
nbr_n_trees <- baseline_nonbreeding$n.trees
nbr_shrinkage <- baseline_nonbreeding$shrinkage
nbr_interaction <- baseline_nonbreeding$interaction.depth

#------------------------------------------------------------------------------
# Step 6: Generate real-time and maximal D-SDMs for each hindcast date
#------------------------------------------------------------------------------

years<-c(2004:2016)

months<-c(1:12)

days<-c(1)

dataframe <- as.data.frame(expand.grid(years, months, days))

quelea_br_data$date <- as.Date(with(quelea_br_data, paste(year, month, day, sep = "-")), "%Y-%m-%d")
quelea_nbr_data$date <- as.Date(with(quelea_nbr_data, paste(year, month, day, sep = "-")), "%Y-%m-%d")

save_dir <- paste0(directory,"SDMs/")

results <- NULL


for (x in 1:nrow(dataframe)) {
  
  yr <- dataframe$Var1[i]
  mn <- dataframe$Var2[i]
  day <- dataframe$Var3[i]
  
  forecast_intitiation <- as.Date(paste0(yr, "-", mn, "-", day)) 
  
  end_date <- forecast_intitiation %m+% months(8)
  
  # Generate dataset for maximal/real-time and breeding/non-breeding
  quelea_br_data_real <- quelea_br_data[quelea_br_data$date < forecast_intitiation, ]
  
  quelea_nbr_data_real <- quelea_nbr_data[quelea_nbr_data$date < forecast_intitiation, ]
  
  quelea_br_data_maximal <- quelea_br_data[!(quelea_br_data$date >= start_date & quelea_br_data$date <= end_date), ]
  
  quelea_nbr_data_maximal <- quelea_nbr_data[!(quelea_nbr_data$date >= start_date & quelea_nbr_data$date <= end_date), ]
  
  br_test  <- quelea_br_data[(quelea_br_data$date >= start_date & quelea_br_data$date <= end_date), ]
  
  nbr_test <- quelea_nbr_data[(quelea_nbr_data$date >= start_date & quelea_nbr_data$date <= end_date), ]
  
  
  ###############################################################
  # Fit real-time D-SDMs for breeding and non-breeding
  ###############################################################
 
  quelea_br_data_real_model <- brt_fit(occ.data = quelea_br_data_real,
                                       response.col = "presence_absence",
                                       varnames = model_variables,
                                       distribution = "bernoulli",
                                       interaction.depth = br_interaction,
                                       n.trees = br_n_trees,	
                                       shrinkage = br_shrinkage,
                                       weights.col = "weights")
  
  quelea_nbr_data_real_model <- brt_fit(occ.data = quelea_nbr_data_real,
                                        response.col = "presence_absence",
                                        varnames = model_variables,
                                        distribution = "bernoulli",
                                        interaction.depth = nbr_interaction,
                                        n.trees = nbr_n_trees,	
                                        shrinkage = nbr_shrinkage,
                                        weights.col = "weights")
  
  ###############################################################
  # Fit maximal D-SDMs for breeding and non-breeding
  ###############################################################

  quelea_br_data_maximal_model <- brt_fit(occ.data = quelea_br_data_maximal,
                                             response.col = "presence_absence",
                                             varnames = model_variables,
                                             distribution = "bernoulli",
                                             interaction.depth = br_interaction,
                                             n.trees = br_n_trees,	
                                             shrinkage = br_shrinkage,
                                             weights.col = "weights")
  
  quelea_nbr_data_maximal_model <- brt_fit(occ.data = quelea_nbr_data_maximal,
                                                response.col = "presence_absence",
                                                varnames = model_variables,
                                                distribution = "bernoulli",
                                                interaction.depth = nbr_interaction,
                                                n.trees = nbr_n_trees,	
                                                shrinkage = nbr_shrinkage,
                                                weights.col = "weights")
  
  
  br_test$pred <- predict(quelea_br_data_maximal_model, newdata = br_test, type = 'response')
  BR_MAX_AUC <- pROC::auc(pROC::roc(presence_absence ~ pred, data = br_test, quiet = T))
  
  br_test$pred <- predict(quelea_br_data_real_model, newdata = br_test, type = 'response')
  BR_REAL_AUC <- pROC::auc(pROC::roc(presence_absence ~ pred, data = br_test, quiet = T))
  
  nbr_test$pred <- predict(quelea_nbr_data_maximal_model, newdata = nbr_test, type = 'response')
  NBR_MAX_AUC <- pROC::auc(pROC::roc(presence_absence ~ pred, data = nbr_test, quiet = T))
  
  nbr_test$pred <- predict(quelea_nbr_data_real_model, newdata = nbr_test, type = 'response')
  NBR_REAL_AUC <- pROC::auc(pROC::roc(presence_absence ~ pred, data = nbr_test, quiet = T))
  
  
  # Record AUC for each D-SDM
  result <- data.frame(date = forecast_intitiation, 
                       BR_MAX_AUC = BR_MAX_AUC,
                       BR_REAL_AUC = BR_REAL_AUC,
                       NBR_MAX_AUC = NBR_MAX_AUC,
                       NBR_REAL_AUC = NBR_REAL_AUC
                       )
  
  results <- rbind(results,result)

  # Save D-SDMs for each hindcast date
  saveRDS(quelea_nonbreeding_data_maximal_model,file = paste0(save_dir, forecast_intitiation, "_NBR_MAX_.rds"))
  saveRDS(quelea_breeding_data_maximal_model,file = paste0(save_dir, forecast_intitiation, "_BR_MAX_.rds"))
  saveRDS(quelea_nonbreeding_data_real_model,file = paste0(save_dir, forecast_intitiation, "_NBR_REAL_.rds"))
  saveRDS(quelea_breeding_data_real_model,file = paste0(save_dir, forecast_intitiation, "_BR_REAL_.rds"))
  
}

# Save D-SDM performances
write.csv(results, file = paste0(save_dir, "DSDM_PERFORMANCES.csv"))
