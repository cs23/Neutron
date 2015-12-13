# Data Harvest Pipeline
rm(list = ls())
source("get_neutron_count_data.R")
source("agg_neutron_count_data.R")
source("join_agg_count_to_spatial.R")
source("join_agg_countspatial_to_weather.R")
bartol_url = "http://neutronm.bartol.udel.edu/~pyle/bri_table.html"
final_data = get_neutron_count_data(bartol_url) %>% 
  agg_neutron_count_data %>% 
  join_agg_count_to_spatial("Delaware_Stations.csv") %>%
  join_agg_countspatial_to_weather
