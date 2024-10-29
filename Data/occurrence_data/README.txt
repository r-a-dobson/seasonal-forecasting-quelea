
Documentation for occurrence datasets

File names: 
- breeding_distribution_quelea.csv - Occurrence records from the breeding season
- nonbreeding_distribution_quelea.csv - Occurrence records from the non-breeding season

Description:
Occurrence records for the red-billed quelea (Quelea quelea) in the breeding and non-breeding seasons between 2003 and 2019. Records were obtained from GBIF (2021) or control operation data sources in southern Africa. Each record represents either an observed species presence or a generated pseudo-absence with associated coordinates and timestamps, as well as calculated weights representing sampling effort in the e-Bird dataset (Auer et al., 2022). 

Columns:

1. day
   - Description: Day of the species occurrence record.
   - Type: Integer
   - Range: 1 to 31

2. month
   - Description: Month of the species occurrence record.
   - Type: Integer
   - Range: 1 to 12

3. year
   - Description: Year of the species occurrence record.
   - Type: Integer
   - Range: 2003 to 2019

4. x
   - Description: Longitude of the species occurrence record.
   - Type: Float
   - Units: Decimal degrees

5. y
   - Description: Latitude of the species occurrence record.
   - Type: Float
   - Units: Decimal degrees

6. data_source
   - Description: Primary source of the occurrence data.
   - Type: String
   - Examples: 'SABAP', 'eBird'
   - Details: The original dataset from which the occurrence data was obtained.

7. data_source2
   - Description: Classification of the data source as control or GBIF data.
   - Type: String
   - Possible Values: 'controldata', 'gbifdata'
   - Details: Further categorises the source of the data to differentiate between control data and GBIF data sources.

8. presence_absence
   - Description: Indicator of species presence (1) or absence (0).
   - Type: Integer (Binary)
   - Possible Values: 1 (Presence), 0 (Absence)
   - Details: Presence indicates that the red-billed quelea was observed at the location and date. Absence indicates a pseudo-absence record.

9. weights
   - Description: Assigned weight for each record.
   - Type: Float
   - Details: The weight is calculated based on the total number of avian e-Bird sampling events within 
              a spatiotemporal buffer around the occurrence record location and date.
   - Units: Relative (no specific units)

Purpose:
These datasets are used in the associated code to extract dynamic explanatory variables and generate dynamic species distribution models (D-SDMs) for the red-billed quelea. 

Notes:
- In both datasets, each row represents a unique species occurrence (2003-2019) or pseudo-absence record with associated metadata.
- Weights can be applied in the D-SDMs to correct for spatiotemporal sampling bias in occurrence records. 

References:
- Auer, T., Barker, S., Borgmann, K., Charnoky, M., Childs, D., Curtis, J., Davies, I., Downie, I., Fink, D., Fredericks, T., Ganger, J., Gerbracht, J., Hanks, C., Hochachka, W., Iliff, M., Imani, J., Johnston, A., Lenz, T., Levatich, T., Ligocki, S., Long, M. T., Morris, W., Morrow, S., Oldham, L., Padilla Obregon, F., Robinson, O., Rodewald, A., Ruiz-Gutierrez, V., Strimas-Mackey, M., Wolf, H., & Wood, C. (2022). EOD – eBird Observation Dataset. Cornell Lab of Ornithology. Occurrence dataset. https://doi.org/10.15468/aomfnb (accessed via GBIF.org on 2021-07-21).

- GBIF Occurrence Download. (2021, July 6). https://doi.org/10.15468/dl.qza9ty
