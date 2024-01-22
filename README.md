# Seasonal forecasting of red-billed quelea distribution

# Summary

In this repository, we provide code and data for generating seasonal forecasts of red-billed quelea (*Quelea* *quelea*) distribution suitability up to seven months ahead. 

Our code is highly flexible. In this first step, you can specify the: 
 - `Spatial extent:` country or countries to forecast across (or your own custom
`sf` polygon). Please note all SDMs are trained on Quelea quelea lathamii
 records from southern Africa countries.

- `Spatial resolution:` given in degrees, this value specifies the spatial
 resolution for the output seasonal forecasts.

- `Temporal extent:` the month and year to initiate the seasonal forecast. ECMWF
 SEAS5 forecasts have an initial date of the 1st of each month, and run for 7
 months. These data are released on the 5th of each month. We recommend
 choosing the closest available month.

 - `Temporal resolution:` the intervals to generate seasonal forecasts at between
 zero and seven months ahead of the initiation date. This can range from daily
 to weekly to monthly intervals depending on your needs.

Once you have specified these details in the first step, the rest of the code
 does not need to be edited aside from giving your Climate Data Store (CDS)
 username and key (Step 2), and your Google Earth Engine/Google Drive email
 (Step 4).

Before beginning the script, please ensure that all necessary packages have been installed and that you have registered for free Climate Data Store,
Google Earth Engine and Google Drive accounts. These are required for download 
of seasonal forecast and historical data used to generate forecasts.

The final output of this script are seasonal forecasts of quelea distribution suitability, exported in both `tif` and `png` format.


# Repository structure

## `seasonal_forecasting_quelea`

- `packages_to_install.txt` - list of packages required for generating seasonal forecasts - you do not need to open or edit this file;

- `Functions_For_Forecast.R` - custom functions for generating forecasts - you do not need to open or edit this code;

- `Seasonal_Forecast_Quelea.R` - This script includes the code for every step in generating seasonal forecasts of red-billed quelea distribution suitability. You will need to open and run this script. 


## `/Data`

- `breeding_distribution_quelea.csv` - filtered breeding season occurrence records for Quelea quelea lathamii in southern Africa with associated dynamic explanatory variables; 

- `nonbreeding_distribution_quelea.csv` - filtered non-breeding season occurrence records for Quelea quelea lathamii in southern Africa with associated dynamic explanatory variables. 


## `/Data/average_length_phenology`

- `Greenup_mean.tif` - average duration (number of days) of the vegetation stage "green-up" extracted from MODIS Land Cover Dynamics Yearly dataset;

- `Maturity_mean.tif` - average duration (number of days) of the vegetation stage "maturity" extracted from MODIS Land Cover Dynamics Yearly dataset; 

- `MidGreendown_mean.tif` - average duration (number of days) of the vegetation stage "mid-green-down" extracted from MODIS Land Cover Dynamics Yearly dataset; 

- `Peak_mean.tif` - average duration (number of days) of the vegetation stage "peak" extracted from MODIS Land Cover Dynamics Yearly dataset; 

- `Senescence_mean.tif` - average duration (number of days) of the vegetation stage "senescence" extracted from MODIS Land Cover Dynamics Yearly dataset; 

FRIEDL, M., GRAY, J. & SULLA-MENASHE, D. 2019. MCD12Q2 MODIS/Terra+Aqua Land Cover Dynamics Yearly L3 Global 500m SIN Grid V006 [MCD12Q2].


## `/Data/models`

- `classification_model_cereal.rds` - Random Forest classification model for projecting vegetation growth stages in cereal croplands land cover cells (Friedl and Sulla-Menashe, 2019) based upon remote-sensed enhanced vegetation index (EVI) characteristics extracted from Didan (2021);

- `classification_model_grass.rds` - Random Forest classification model for projecting vegetation growth stages in grasslands based upon remote-sensed EVI characteristics.

DIDAN, K. 2021. MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061 [Data set]. NASA EOSDIS Land Processes DAAC. 

FRIEDL, M. & SULLA-MENASHE, D. 2019. MCD12Q1 MODIS/Terra+ Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V006 [MCD12Q1].


# Help with Seasonal Forecasting 

If you encounter an error or bug when generating seasonal forecasts using our code and data, please post a comment
[here](https://github.com/r-a-dobson/seasonal-forecasting-quelea/issues) for guidance and
support from us.
