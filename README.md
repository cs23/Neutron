# Neutron
Neutron Flux Modeling

Main Purpose: Development of spatial model for neutron count/flux monitoring in STAT 575 Final Project.

There will be two main pipelines, data harvest and statistical analysis.

###Data Harvest Pipeline: *COMPLETE* 
`source("data_harvest_pipeline.R")` will gather neutron counts aggregated by day for the 8 Bartol stations in 2014.

####1. `get_neutron_count_data()` 

#####Purpose
Grabs ASCII format data from individual links at Bartol ftp site. Also cleans the data up according to Hadley Wickam standards.

#####Required Input
Location of Bartol Data, [Bartol Research Institute Neutron Monitor Data](http://neutronm.bartol.udel.edu/~pyle/bri_table.html)  

#####Optional Input
Year, (default is 2014), accepts a single year `get_neutron_count_data(years = 2000)` or continous range of years `get_neutron_count_data(years = 2010:2014)`
  
####2. `agg_neutron_count_data()` 

#####Purpose
Takes hourly count data and aggregates by specified time unit.

#####Required Input
Output from `get_neutron_count_data()`

#####Optional Input
Aggregation unit, (default is Day), accepts month `agg_neutron_count_data(agg_unit = "Month")` or year `agg_neutron_count_data(agg_unit = "Year")`

_Future improvements will allow aggregation by partial day units (3, 6, or 12 hours). Not part of this class!_
  
####3. `join_agg_count_to_spatial()`

#####Purpose
Takes aggregated count data and joins it with "clean" geospatial station information.

#####Required Input
Output from `agg_neutron_count_data()`
Location of Station Data, `Delaware_Stations.csv`

####4. `join_agg_countspatial_to_weather()`

#####Purpose
Takes aggregated count and spatial data and finds the weather for the specified times. Next, it merges weather data with the count data for a final data frame for analysis. *Requires [Forecast.io API Key](https://developer.forecast.io/)*. Also checks to see that package `Rforecastio` is installed, see [Rforecastio on GitHub](https://github.com/hrbrmstr/Rforecastio) for further information.

#####Required Input
Output from `join_agg_count_to_spatial()`

#####Optional Input
If using data from this website for 2014 aggregated by day, the previously saved `2014_Daily_Weather.csv` will automatically be loaded. For different years, new weather data will have to be downloaded and saved locally. Please set arguments `use_saved_weather = FALSE` and adjust `get_saveName` to reflect the year choice. Afterwards, set `used_saved_weather = TRUE`.

###Statistical Analysis Pipeline: *COMING SOON*

####General Ideas
#####1. Geostatistical analysis with kriging
Test prediction capabilities with cross-validation (leave 1 station out) for all 8 stations.

#####2. Visualizations
For usage in final paper and presentation.
