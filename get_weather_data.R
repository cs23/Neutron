# Check to see if Rforecastio package is installed. This package allows API access from R
# to the forecast.io database.
# https://github.com/hrbrmstr/Rforecastio
if ("Rforecastio" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("hrbrmstr/Rforecastio")
}
library(Rforecastio)
# User will need to enter a forecast.io API key which can be generated 
# from the following site: https://developer.forecast.io/
# 1000 API calls per day are free, but subsequent calls cost $0.00001 
# (10k per dollar). Setting up an account requires an email. If user 
# does not provide credit card information, then API calls are blocked.

# This function gets weather data for aggregated counts. Daily weather information is
# given by "Day". Future updates will include getting hourly weather data. Saving the
# weather data from forecast.io API calls has two options: csv (default) or RData 
# (binary) along with the choice of a name.
get_weather_data = function(joined_agg_data, agg_unit = "Day",
                            saveName = "2014_Daily_Weather", saveOption = "csv") {
  agg_unit = "currently,minutely,hourly,alerts,flags"
  joined_agg_data = joined_agg_data %>% 
    mutate(Forecastio_Time = paste0(Year, "-", Month, "-", Day, "T00:00:00-0000"))
  forecast_for_agg_data = tbl_df(data.frame())
  for (i in 1:nrow(joined_agg_data)) {
    cat("Grabbing daily forecast", i, "of", nrow(joined_agg_data), "\n")
    this_forecast = get_forecast_for(joined_agg_data$Lat[i], 
                                     joined_agg_data$Lon[i], 
                                     joined_agg_data$Forecastio_Time[i], 
                                     units = "si", exclude = agg_unit)
    forecast_for_agg_data = bind_rows(forecast_for_agg_data, this_forecast$daily)
  }
  if (saveOption == "csv") {
    write_csv(forecast_for_agg_data, paste(saveName, saveOption, sep = "."))
  } else {
    save(forecast_for_agg_data, file = paste0(saveName, ".RData"))
  }
  return(forecast_for_agg_data)
}

# Loads previously saved weather data saved by get_weather_data()
load_previous_weather_data = function(previous_weather_data_location) {
  if (grepl("csv$", previous_weather_data_location)) {
    forecast_for_agg_data = read_csv(previous_weather_data_location)
  } else {
    forecast_for_agg_data = loadRData(previous_weather_data_location)
  }
  return(forecast_for_agg_data)
}
# Testing
# previous_weather_data_location = "2014_Daily_Weather.csv"
# test1 = load_previous_weather_data(previous_weather_data_location) # Works correctly!
# previous_weather_data_location = "2014_Daily_Weather.RData"
# test2 = load_previous_weather_data(previous_weather_data_location) # Works correctly!
# identical(test1, test2) # Both are identical!

# Helper function that loads an RData file, and returns it. Used in 
# load_previous_weather_data()
loadRData = function(filename) {
  load(filename)
  get(ls()[ls() != "filename"])
}


# test_joined_agg_data = joined_agg_data


# agg_unit = "currently,minutely,hourly,alerts,flags"
# forecast_test2 = get_forecast_for(-77.90, 166.60, "2014-01-01T00:00:00-0100", units = "si", exclude = agg_unit)


# joined_agg_data = joined_agg_data %>% mutate(Forecastio_Time = paste0(Year, "-", Month, "-", Day, "T00:00:00-0000"))

# all_forecast = tbl_df(data.frame())
# for (i in 1:nrow(joined_agg_data)) {
#   cat("Grabbing daily forecast", i, "of", nrow(joined_agg_data), "\n")
#   this_forecast = get_forecast_for(joined_agg_data$Lat[i], joined_agg_data$Lon[i], joined_agg_data$Forecastio_Time[i], units = "si", exclude = agg_unit)
#   all_forecast = bind_rows(all_forecast, this_forecast$daily)
# }
# save(all_forecast, file = "2014_Daily_Weather.RData")
# write_csv(all_forecast, "2014_Daily_Weather.csv")
# load("2014_Daily_Weather.RData")



# now <- get_current_forecast(43.2672, -70.8617)


# specific_time = "2013-05-06T12:00:00-0400"
# 
# tmp = get_forecast_for(-90, 0, specific_time)
# tmp$hourly
# 
# then1 <- get_forecast_for(43.07, -89.401, "2013-05-06T12:00:00-0400")
# print(then$daily)
# 
# 
# then = get_forecast_for(-77.90, 166.60, "2014-01-01T00:00:00-0000")
# print(then$daily)

# test2 = get_forecast_for(43.267)
# 
# more_than_one <- data.frame(loc=c("Maine", "Seattle"),
#                             lon=c(43.2672, 47.6097),
#                             lat=c(70.8617, 122.3331),
#                             when=c("2013-05-06T12:00:00-0400",
#                                    "2013-05-06T12:00:00-0400"),
#                             stringsAsFactors=FALSE)
# 
# bigger_list <- mapply(get_forecast_for, 
#                       more_than_one$lon, more_than_one$lat, more_than_one$when,
#                       SIMPLIFY=FALSE)
# names(bigger_list) <- more_than_one$loc
# 
# bigger_list$Seattle$daily
