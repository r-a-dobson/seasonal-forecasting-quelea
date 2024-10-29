# Near-term forecasting of red-billed quelea distribution

# Summary

In this repository, we provide code and data for generating near-term hindcasts of red-billed quelea (*Quelea* *quelea*) distribution suitability up to seven months ahead. 

Our code is flexible, you can specify the: 
 - `Spatial extent:` country or countries to generate seasonal hindcasts across (or provide your own custom
`sf` polygon). Please note all species distribution models (SDMs) are trained on *Quelea quelea lathamii*
 records from southern Africa.

- `Spatial resolution:` given in degrees, this value specifies the spatial
 resolution for the output seasonal forecasts.

Before beginning the script, please ensure that all necessary packages have been installed and that you have registered for free [Climate Data Store](https://cds.climate.copernicus.eu/#!/home),
[Google Earth Engine](https://developers.google.com/earth-engine/) and [Google Drive](https://www.google.co.uk/intl/en-GB/drive/) accounts. These are required for download 
of seasonal forecast and historical datasets that are used to generate the seasonal forecasts.

The final outputs of the script are seasonal hindcasts of quelea distribution suitability in `tif` format.

# Repository structure

## `seasonal_forecasting_quelea`

- `packages_to_install.txt` - list of packages required for generating near-term hindcasts- you do not need to open or edit this file;

- `Functions_For_Forecast.R` - custom functions for generating near-term hindcasts - you do not need to open or edit this code;

- `S1_Forecast_Weather_Variables.R` - the code for generating near-term hindcasts of 8- and 52-week weather variables. You will need to open and run this script. 

- `S2_Forecast_Resource_Variables.R` - the code for generating near-term hindcasts of resource variables. You will need to open and run this script. 

- `S3_Fit_Models.R` - the code for generating dynamic species distribution models for projecting red-billed quelea distribution suitability. You will need to open and run this script. 

- `S4_Project_Models.R` - the code projecting near-term distribution suitability for quelea using hindcast weather and resource variables. You will need to open and run this script. 


## `/Data`

- `breeding_distribution_quelea.csv` - filtered breeding season occurrence records for *Quelea* *quelea* *lathamii* in southern Africa; 

- `nonbreeding_distribution_quelea.csv` - filtered non-breeding season occurrence records for *Quelea* *quelea* *lathamii* in southern Africa.

- `README.txt` - description of occurrence data frames.

## `/Data/average_length_phenology`

- `Greenup_mean.tif` - average duration (number of days) of the vegetation stage "green-up" extracted from MODIS Land Cover Dynamics Yearly dataset (Friedl et al., 2019);

- `Maturity_mean.tif` - average duration (number of days) of the vegetation stage "maturity" extracted from MODIS Land Cover Dynamics Yearly dataset; 

- `MidGreendown_mean.tif` - average duration (number of days) of the vegetation stage "mid-green-down" extracted from MODIS Land Cover Dynamics Yearly dataset; 

- `Peak_mean.tif` - average duration (number of days) of the vegetation stage "peak" extracted from MODIS Land Cover Dynamics Yearly dataset; 

- `Senescence_mean.tif` - average duration (number of days) of the vegetation stage "senescence" extracted from MODIS Land Cover Dynamics Yearly dataset; 

FRIEDL, M., GRAY, J. & SULLA-MENASHE, D. 2019. MCD12Q2 MODIS/Terra+Aqua Land Cover Dynamics Yearly L3 Global 500m SIN Grid V006 [MCD12Q2].


## `/Data/models`

- `classification_model_cereal.rds` - Random Forest classification model for classifying vegetation growth stages in cereal croplands land cover cells (Friedl and Sulla-Menashe, 2019) based upon remote-sensed enhanced vegetation index (EVI) characteristics extracted from Didan (2021);

- `classification_model_grass.rds` - Random Forest classification model for classifying vegetation growth stages in grasslands based upon remote-sensed EVI characteristics.

DIDAN, K. 2021. MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061 [Data set]. NASA EOSDIS Land Processes DAAC. 

FRIEDL, M. & SULLA-MENASHE, D. 2019. MCD12Q1 MODIS/Terra+ Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V006 [MCD12Q1].


# Help with Seasonal Forecasting 

If you encounter an error or bug when generating seasonal forecasts using our code and data, please post a comment
[here](https://github.com/r-a-dobson/seasonal-forecasting-quelea/issues) for guidance and
support from us.
