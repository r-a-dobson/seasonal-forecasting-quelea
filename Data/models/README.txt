
Documentation for Random Forest models

File Names:
- classification_model_cereal.rds - Model trained to classify vegetation phenology stages for "cereal cropland" cells.
- classification_model_grass.rds - Model trained to classify vegetation phenology stages for "grassland" cells.

Description:
These Random Forest (RF) models were trained to predict vegetation phenology stages based on dynamic EVI characteristics extracted from 
the MODIS 16-day EVI dataset (Didan, 2021).

Model Training: 
Historical seed abundance values were extracted from the MODIS Land Cover Dynamics dataset (Friedl et al., 2019), which provides global data on 
seven vegetation growth stages derived from EVI data for each 500m cell. Using Random Forests, the phenology stages were classified based on 
EVI characteristics. Training data consisted of 5000 randomly sampled cells from each year (2002–2019), with vegetation growth stages and 
associated EVI data extracted for each cell.

Model Explanatory Variables:
- Peak in last year: Maximum EVI in the past fifty-two weeks.
- Trough in last year: Minimum EVI in the past fifty-two weeks.
- Time since peak: Weeks since the last peak EVI.
- Time since trough: Weeks since the last trough EVI.
- Current EVI: Most recent 16-day EVI value before forecast date.
- Change in past four weeks: Difference between current EVI and the EVI four weeks prior.
- Change in past eight weeks: Difference between current EVI and the EVI eight weeks prior.
- Change in past sixteen weeks: Difference between current EVI and the EVI sixteen weeks prior.
- Percentage of amplitude: Ratio of current EVI to EVI amplitude (peak minus trough from the last year) multiplied by 100.

References:
- Didan, K. (2021). MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061 [Data set]. NASA EOSDIS Land Processes DAAC. 
  https://doi.org/10.5067/MODIS/MOD13Q1.061
- Friedl, M., Gray, J., & Sulla-Menashe, D. (2019). MCD12Q2 MODIS/Terra+Aqua Land Cover Dynamics Yearly L3 Global 500m SIN Grid V006 [Data set]. 
  In N. E. L. P. DAAC (Ed.). https://doi.org/10.5067/MODIS/MCD12Q2.061
- Friedl, M., & Sulla-Menashe, D. (2019). MCD12Q1 MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V006 [MCD12Q1]. 
  In N. E. L. P. DAAC (Ed.). https://doi.org/10.5067/MODIS/MCD12Q1.061

Usage:
Load these models into R using readRDS("path/to/model.rds") and apply them to new data using the predict() function.

