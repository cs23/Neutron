# Neutron
Neutron Flux Modeling

Main Purpose: Development of spatial model for neutron count monitoring in STAT 575 Final Project

##Data Harvest Pipeline: *COMPLETE*. 
Running this (`source("data_harvest_pipeline.R")`) will gather neutron counts aggregated by day for the 8 Bartol stations in 2014.
This contains:
  * `get_neutron_count_data()`
  * Weather data `get_neutron_weather_data()`
  * Spatial data

Statistical Analysis Pipeline
  * Geostatistical analysis with kriging
    - Test prediction capabilities with cross-validation (leave 1 station out) for all 8 stations then average
  * Visualizations
