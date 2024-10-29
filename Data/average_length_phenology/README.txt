
Documentation for vegetation phenology stage mean duration TIF files

File Names:
- Greenup_mean.tif: Mean duration of the "Greenup" stage (in days, integer).
- Maturity_mean.tif: Mean duration of the "Maturity" stage (in days, integer).
- MidGreendown_mean.tif: Mean duration of the "MidGreendown" stage (in days, integer).
- MidGreenup_mean.tif: Mean duration of the "MidGreenup" stage (in days, integer).
- Peak_mean.tif: Mean duration of the "Peak" stage (in days, integer).
- Senescence_mean.tif: Mean duration of the "Senescence" stage (in days, integer).

Description:
These .tif files provide the average duration for specific vegetation phenology stages across the extent of southern Africa. 
The mean durations were derived from the MODIS Land Cover Dynamics dataset (Friedl et al., 2019), which classifies vegetation growth stages based on 
the Enhanced Vegetation Index. Each file represents the average number of days for a specific phenology stage based on 
historical data (2002-2019). 

Spatial projection and coverage:
- Projection: Geographic Coordinate System (longitude and latitude) 
- Coverage: Southern Africa

References:
- Didan, K. (2021). MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061 [Data set]. NASA EOSDIS Land Processes DAAC. 
  https://doi.org/10.5067/MODIS/MOD13Q1.061
- Friedl, M., Gray, J., & Sulla-Menashe, D. (2019). MCD12Q2 MODIS/Terra+Aqua Land Cover Dynamics Yearly L3 Global 500m SIN Grid V006 [Data set]. 
  In N. E. L. P. DAAC (Ed.). https://doi.org/10.5067/MODIS/MCD12Q2.061

Usage:
These files can be loaded into R through the `terra` package using `terra::rast("path/to/file.tif")`.
