# Neutron
Neutron Flux Modeling

Main Purpose: Development of spatial model for neutron count monitoring.

Data Harvest Pipeline
  * Count data `get_neutron_count_data()`
  * Weather data `get_neutron_weather_data()`
  * Spatial data

Aggregation Pipeline
  * Choose time interval to aggregate by (default is day)

Statistical Analysis Pipeline
  * Geostatistical analysis with kriging
    - Test prediction capabilities with cross-validation (leave 1 station out)
  * Visualizations
