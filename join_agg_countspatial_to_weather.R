source("get_weather_data.R")
# This function takes the aggregated spatial neutron counts and gets the weather in two ways:
# If saved previously, load_previous_weather_data() is used. Otherwise, get_weather_data() calls
# the forecast.io API based on the level of aggretation of the spatial neutron counts. Finally,
# the weather data is binded to aggregated spatial neutron counts for statistical analysis.
join_agg_countspatial_to_weather = function(neutron_counts_agg_spatial, use_saved_weather = TRUE, 
                                            previous_weather_data_location = "2014_Daily_Weather.csv",
                                            get_saveName = "2014_Daily_Weather", get_saveOption = "csv") {
  if (use_saved_weather == TRUE) {
    weather_data = load_previous_weather_data(previous_weather_data_location)
  } else {
    weather_data = get_weather_data(joined_agg_data, agg_unit = "Day",
                                    saveName = get_saveName, saveOption = get_saveOption) 
  }
  return(bind_cols(neutron_counts_agg_spatial, weather_data))
}

# Testing
# final_data = join_agg_countspatial_to_weather(neutron_counts_agg_spatial)