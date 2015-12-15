##########UNDER CONSTRUCTION#############

source("data_harvest_pipeline.R")


summary(final_data)

selected_data = final_data %>% select(Year:Pressure, Lat:Alt, precipIntensity, dewPoint,
                                      humidity, visibility, cloudCover, temperatureMin, temperatureMax)

selected_data %>% na.omit %>% group_by(Station) %>% tally()
selected_data %>% group_by(Station) %>% tally()


final_data %>% filter(precipAccumulation == 0)


# precipIntensity, dewPoint humidity, visibility, cloudCover, temperatureMin, temperatureMax

# Description of weather variables from forecast.io API documentation
# https://developer.forecast.io/docs/v2

# precipIntensity: A numerical value representing the average expected intensity (in inches of liquid
# water per hour) of precipitation occurring at the given time conditional on probability (that is, 
# assuming any precipitation occurs at all). A very rough guide is that a value of 0 in./hr. 
# corresponds to no precipitation, 0.002 in./hr. corresponds to very light precipitation, 0.017 in./hr.
# corresponds to light precipitation, 0.1 in./hr. corresponds to moderate precipitation, and 0.4 in./hr.
# corresponds to heavy precipitation.

# dewPoint: A numerical value representing the dew point at the given time in degrees Fahrenheit.

# humidity: A numerical value between 0 and 1 (inclusive) representing the relative humidity.

# visibility: A numerical value representing the average visibility in miles, capped at 10 miles.

# cloudCover: A numerical value between 0 and 1 (inclusive) representing the percentage of sky 
# occluded by clouds. A value of 0 corresponds to clear sky, 0.4 to scattered clouds, 0.75 to 
# broken cloud cover, and 1 to completely overcast skies.

# temperatureMin, temperatureMinTime, temperatureMax, and temperatureMaxTime (only defined on daily 
# data points): numerical values representing the minimum and maximumum temperatures (and the UNIX 
# times at which they occur) on the given day in degrees Fahrenheit.

# General analysis strategy
# 5 stations that have complete values kriging for single day, do leave one out analysis and predict missing one
# 8 stations that have only pressure and counts kriging for single day, do leave one out analysis and predict missing one
# 8 stations do counts kriging for single day, do leave one out analysis and predict missing one

library(gstat)
library(sp)
# library(spacetime)
# library(raster)
# library(rgdal)
# library(rgeos)

first_day = selected_data %>% filter(Year == "2014", Month == "01", Day == "01")
spatial_selected_data = as.data.frame(first_day)


coordinates(spatial_selected_data) = ~ Lon + Lat

bbox(spatial_selected_data)
proj4string(spatial_selected_data)

counts_vgm = variogram(Corr_Count ~ Pressure + Alt, spatial_selected_data)
plot(counts_vgm)

selected_data

REG1 = lm(Corr_Count ~ Alt + Pressure, data = selected_data)
summary(REG1)

plot_data = selected_data %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% select(Date, Station, Corr_Count, Pressure, Alt)

ggplot(plot_data, aes(x = Date, y = Corr_Count)) + geom_line() + scale_x_date() + facet_wrap(~ Station)

ggplot(plot_data, aes(x = Pressure, y = Corr_Count, color = Station, group = Station)) + geom_point() + coord_flip()

library(maps)
world = map_data("world")
p = ggplot(legend=FALSE) +
  geom_polygon(data=world, aes(x=long, y=lat,group=group), color = "gray60", fill = "gray95") + 
  theme_bw() + 
  labs(x = "Longitude",
       y = "Latitude",
       title = "Location of Bartol Neutron Monitors")

locations = selected_data %>% group_by(Station) %>% distinct() %>% select(Station, Lat:Lon)

graphic1 = p + geom_point(data = locations, aes(x = Lon, y = Lat), color = "red", size = 4) +
  geom_text(data = locations, 
            aes(label = Station, x = Lon, y = Lat-5, size = 10, color = "darkred"), show_guide = FALSE)

ggsave("neutron_map.pdf", graphic1)
