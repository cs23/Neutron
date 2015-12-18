##########UNDER CONSTRUCTION#############

source("data_harvest_pipeline.R")


# summary(final_data)

# selected_data = final_data %>% 
#   select(Year:Pressure, Lat:Alt, precipIntensity, dewPoint,
#          humidity, visibility, cloudCover, temperatureMin, temperatureMax) %>%
#   mutate(temperature = (temperatureMax + temperatureMin) / 2) %>% 
#   select(Year:cloudCover, temperature) %>% 
#   mutate(Date = as.POSIXct(paste(Year, Month, Day, sep = "-")))


selected_data = final_data %>% 
  select(Year:Pressure, Lat:Alt, precipIntensity, dewPoint,
         humidity, visibility, cloudCover, temperatureMin, temperatureMax) %>%
  mutate(temperature = (temperatureMax + temperatureMin) / 2) %>% 
  select(Year:cloudCover, temperature) %>% 
  mutate(Date = as.POSIXct(paste(Year, Month, Day, sep = "-"))) %>% mutate(Date = factor(Date))


library(maps)
world = map_data("world")
p = ggplot(legend=FALSE) +
  geom_polygon(data=world, aes(x=long, y=lat,group=group), color = "gray60", fill = "gray95") + 
  theme_bw() + 
  labs(x = "Longitude",
       y = "Latitude")

locations = selected_data %>% group_by(Station) %>% distinct() %>% select(Station, Lat:Lon)

graphic1 = p + geom_point(data = locations, aes(x = Lon, y = Lat), color = "red", size = 4) +
  geom_text(data = locations, 
            aes(label = Station, x = Lon, y = Lat-5, size = 10, color = "darkred"), show_guide = FALSE)
# ggsave("neutron_map.pdf", graphic1)

plot_data = selected_data %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% select(Date, Station, Corr_Count, Pressure, Alt)

graphic2 = ggplot(plot_data, aes(x = Date, y = Corr_Count, color = Station, group = Station)) + geom_line() + scale_x_date() +
  labs(y = "Correlated Counts (per hour)")
# ggsave("counts_timeseries.pdf", graphic2)

library(ape)
library(geosphere)
all_moran = NULL
for (i in 1:nlevels(selected_data$Date)) {
  subset_data = selected_data %>% filter(Date == levels(Date)[i])
  subset_data = as.data.frame(subset_data)
  coordinates(subset_data) = ~ Lon + Lat
  proj4string(subset_data) = CRS("+init=epsg:3348")
  stations_dists.inv = 1/distm(subset_data@coords)
  diag(stations_dists.inv) = 0
  moran_output = Moran.I(subset_data@data$Corr_Count, stations_dists.inv)
  all_moran = rbind(all_moran, unlist(moran_output))
}
all_dates = tbl_df(data.frame(date = as.vector(levels(selected_data$Date)), stringsAsFactors = FALSE))
all_moran = tbl_df(as.data.frame(all_moran))
final_moran = bind_cols(all_moran, all_dates) %>% mutate(date = as.Date(date))

graphic3 = ggplot(final_moran, aes(x = date, y = p.value)) + 
  geom_point() + 
  scale_x_date() + 
  geom_hline(aes(yintercept = .05), color = "red") +
  labs(x = "Date",
       y = "Moran's I P-value")

library(xtable)
selected_data = selected_data %>% rename(pressure = Pressure) %>% rename(altitude = Alt)

summary(test <- lm(Corr_Count ~ pressure + altitude + precipIntensity + dewPoint + humidity + visibility + temperature, data = selected_data))
selected_data2 = selected_data %>% filter(Corr_Count != 0)
summary(test2 <- lm(Corr_Count ~ pressure + altitude + precipIntensity + dewPoint + humidity + visibility + temperature, data = selected_data2))
summary(test2_null <- lm(Corr_Count ~ altitude + pressure + precipIntensity + humidity + visibility, data = selected_data2))
anova(test2_null, test2)
# qqnorm(resid(test2))
# qqline(resid(test2))
selected_data3 = selected_data %>% na.omit %>% filter(Corr_Count != 0)
summary(test3 <- lm(Corr_Count ~ altitude + pressure + precipIntensity + dewPoint + humidity + visibility + temperature, data = selected_data3))
summary(test3_null <- lm(Corr_Count ~ altitude + pressure + precipIntensity + humidity + visibility, data = selected_data3))
anova(test3_null, test3)

# Do predictions on Madison, WI
# as.vector(levels(selected_data$Date))
madlat = 43.071240
madlon = -89.400137
madalt = 265
new_time = selected_data %>% select(Year:Day)

madison_df = new_time %>% group_by(Year, Month, Day) %>% distinct %>% mutate(Lat = madlat, Lon = madlon)

# Run the get_weather_data only once! Then use load_previous_weather_data
# madison_weather = madison_df %>% get_weather_data(saveName = "2014_Madison_Weather")

madison_weather = load_previous_weather_data("2014_Madison_Weather.csv")

madison_prediction = madison_weather %>% select(pressure, precipIntensity, dewPoint, humidity, visibility, temperatureMin, temperatureMax) %>% 
  mutate(temperature = (temperatureMax + temperatureMin) / 2) %>% select(pressure:visibility, temperature) %>% mutate(altitude = madalt) %>%
  select(altitude, pressure, precipIntensity:temperature)

madison_prediction_output = predict(test3, madison_prediction)
madison_prediction_output_df = tbl_df(data.frame(Corr_Count = unname(madison_prediction_output), Date = as.Date(levels(selected_data$Date))))
madison_prediction_output_df = madison_prediction_output_df %>% mutate(Station = "Madison") 

# plot_data = selected_data %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>% select(Date, Station, Corr_Count, Pressure, Alt)

with_madison_predictions = bind_rows(plot_data %>% select(Corr_Count, Date, Station), madison_prediction_output_df)

graphic4 = ggplot(with_madison_predictions, aes(x = Date, y = Corr_Count, color = Station, group = Station)) + geom_line() + scale_x_date() +
  labs(y = "Correlated Counts (per hour)")


# qqnorm(resid(test3))
# qqline(resid(test3))

# library(sp) 
# library(gstat)
# m = vgm(psill=5e6,"Mat",range=10,nugget=0,kappa=2)
# subset_data = selected_data %>% filter(Date == "2014-01-01", Station %in% canada)
# subset_data = as.data.frame(subset_data)
# coordinates(subset_data) = ~ Lon + Lat
# proj4string(subset_data) = CRS("+init=epsg:3348")
# # 
# counts_first_vgm = variogram(Corr_Count ~ 1, subset_data)
# 
# counts_first_vgm = variogram(Corr_Count ~Lon+Lat, subset_data, alpha=c(0,90))
# 
# plot(counts_first_vgm,  vgm(psill=5e6,"Mat",range=500,nugget=0,kappa=2))
# plot(counts_first_vgm)
# 


# 
# vgm(psill=.2,"Mat",range=10,nugget=.10,kappa=1.5)
# 
# 
# 
# variogram(Corr_Counts~1, subset_data)
# 
# 
# x = krige.cv(Corr_Counts~1, subset_data, m, nfold=8)
# bubble(x, "residual", main = "log(zinc): 5-fold CV residuals")

# stations_dists.inv = 1/distm(subset_data@coords)
# diag(stations_dists.inv) = 0





# used in lighting talk
# original data
selected_data %>% 
  group_by(Station) %>% 
  tally
# used in lighting talk
# remove missing weather data
selected_data %>% 
  na.omit %>% 
  group_by(Station) %>%
  tally
# used in lighting talk
# keep only temperature from forecast.io -- complete cases
selected_data2 = selected_data %>% 
  select(-cloudCover, -visibility, 
         -precipIntensity, -dewPoint, 
         -humidity)

# used in lighting talk
selected_data2 %>% 
  na.omit %>% 
  group_by(Station) %>%
  tally


# Complete cases dataset where Station "Nain" has the smallest amount of observations (142)
nain_only_dates = selected_data %>% na.omit %>% filter(Station == "Nain") %>% select(Date) %>% distinct
complete_cases = selected_data %>% filter(Date %in% as.vector(nain_only_dates$Date))


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

# library(gstat)
# library(sp)
# library(spacetime)
# library(raster)
# library(rgdal)
# library(rgeos)

# test_selected_data = selected_data %>% select(Date, Station:Alt)
# first = test_selected_data %>% filter(Date == "2014-01-01")
# spatial_first = as.data.frame(first)
# coordinates(spatial_first) = ~ Lon + Lat


# canada = c("Fort Smith", "Inuvik", "Nain", "Peawanuck", "Thule")
# 
# k_selected_data = selected_data %>% filter(Station %in% canada)
# 
# for (i in 1:nlevels(k_selected_data$Date)) {
#   sub_k_selected_data = selected_data %>% filter(Date == levels(Date)[i])
#   sub_k_selected_data = as.data.frame(sub_k_selected_data)
#   coordinates(sub_k_selected_data) = ~ Lon + Lat
#   proj4string(sub_k_selected_data) = CRS("+init=epsg:3348")
#   # NAD83 / Statistics Canada Lambert Projection
#   # http://spatialreference.org/ref/epsg/3348/
#   sub_k_selected_data.clambert = spTransform(sub_k_selected_data, CRS("+init=epsg:3348"))
#   # five-fold cross validation:
#   x <- krige.cv(Corr_Count~Lon + Lat, sub_k_selected_data, m, nmax = 40, nfold=155)
#   bubble(x, "residual", main = "log(zinc): 5-fold CV residuals")
#   
# }
# 
# counts_vgm = variogram(Corr_Count ~ 1, sub_k_selected_data)
# plot(counts_vgm)
# plot(counts_vgm)




# k_selected_data = as.data.frame(k_selected_data)
# coordinates(k_selected_data) = ~ Lon + Lat
# proj4string(k_selected_data) = CRS("+init=epsg:4326")
# k_selected_data.mercator = spTransform(k_selected_data, CRS("+init=epsg:32663"))
# k_selected_data.miller = spTransform(k_selected_data, CRS("+init=epsg:3348"))
# k_selected_data.miller = spTransform(k_selected_data, CRS("+init=esri:54003"))
# countsSP = SpatialPoints(k_selected_data.miller@coords, CRS("+init=esri:54003"))
# countsSP = SpatialPoints(k_selected_data.miller@coords, CRS("+init=epsg:3348"))
# dupl = zerodist(countsSP)
# countsDF = data.frame(Counts = k_selected_data.miller@data$Corr_Count)
# countsTM = as.POSIXct(k_selected_data.miller@data$Date)
# timeDF = STIDF(countsSP, countsTM, data = countsDF)
# stplot(timeDF)

# var = variogramST(Counts ~ 1, data = timeDF, tunit = "hours", assumeRegular = FALSE, na.omit = TRUE)


# projection(spatial_first)=CRS("+init=epsg:3857") #WGS84 Web Mercator (Auxiliary Sphere) Google Maps projection
# counts_first_vgm = variogram(Corr_Count ~ Pressure + Alt, spatial_first)
# plot(counts_first_vgm,  vgm(psill=5e6,"Mat",range=10,nugget=0,kappa=2))
# 
# vgm(psill=.2,"Mat",range=10,nugget=.10,kappa=1.5)


# selected_data = selected_data %>% select(Date, Station:temperature)


# spatial_selected_data = as.data.frame(selected_data)
# coordinates(spatial_selected_data) = ~ Lon + Lat
# 
# projection(spatial_selected_data)=CRS("+init=epsg:4326") #WGS84 projection
# projection(spatial_selected_data)=CRS("+init=epsg:3857") #WGS84 Web Mercator (Auxiliary Sphere) Google Maps projection
# projection(spatial_selected_data)=CRS("+init=epsg:3395") #WGS84 World Mercator
# 
# sp = cbind(x = c(0,0,1), y = c(0,1,1))
# row.names(sp) = paste("point", 1:nrow(sp), sep="")
# library(sp)
# sp = SpatialPoints(sp)
# time = as.POSIXct("2010-08-05")+3600*(10:13)
# m = c(10,20,30) # means for each of the 3 point locations
# mydata = rnorm(length(sp)*length(time),mean=rep(m, 4))
# IDs = paste("ID",1:length(mydata))
# mydata = data.frame(values = signif(mydata,3), ID=IDs)
# 
# 
# stidf = as(STFDF(sp, time, mydata), "STIDF")
# stidf[1:2,]
# all.equal(stidf, stidf[stidf,])
# 
# stfdf[1:2,]
# stfdf[,1:2]
# stfdf[,,2]
# stfdf[,,"values"]
# stfdf[1,]
# stfdf[,2]

# (spatial_selected_data) = CRS("+init=epsg:4326")
# neutronSP = SpatialPoints(spatial_selected_data@coords)
# dupl = zerodist(neutronSP)
# 
# neutronDF = data.frame(counts = )
# 
# 
# bbox(spatial_selected_data)
# proj4string(spatial_selected_data)
# 
# counts_vgm = variogram(Corr_Count ~ Pressure + Alt, spatial_selected_data)
# plot(counts_vgm)

# selected_data
# 
# REG1 = lm(Corr_Count ~ Alt + Pressure, data = selected_data)
# summary(REG1)



# ggplot(plot_data, aes(x = Pressure, y = Corr_Count, color = Station, group = Station)) + geom_point() + coord_flip()



# library(sp)
# library(spacetime)
# sumMetricVgm <- vgmST("sumMetric",
#                       space=vgm( 4.4, "Lin", 196.6,  3),
#                       time =vgm( 2.2, "Lin",   1.1,  2),
#                       joint=vgm(34.6, "Exp", 136.6, 12),
#                       stAni=51.7)

# data(air)
# if (!exists("rural")) rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
# 
# rr <- rural[,"2005-06-01/2005-06-03"]
# rr <- as(rr,"STSDF")
# 
# x1 <- seq(from=6,to=15,by=1)
# x2 <- seq(from=48,to=55,by=1)
# 
# DE_gridded <- SpatialPoints(cbind(rep(x1,length(x2)), rep(x2,each=length(x1))), 
#                             proj4string=CRS(proj4string(rr@sp)))
# gridded(DE_gridded) <- TRUE
# DE_pred <- STF(sp=as(DE_gridded,"SpatialPoints"), time=rr@time)
# DE_kriged <- krigeST(PM10~1, data=rr, newdata=DE_pred,
#                      modelList=sumMetricVgm)
# gridded(DE_kriged@sp) <- TRUE
# stplot(DE_kriged)


# gcd_hf = function(lon1, lat1, lon2, lat2) {
#   R = 6371 # Earth mean radius [km]
#   delta_lon = (lon2 - lon1)
#   delta_lat = (lat2 - lat1)
#   a = sin(delta_lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta_lon/2)^2
#   c = 2 * asin(min(1, sqrt(a)))
#   d = R * c
#   return(d) # Distance in km
# }
# 
# deg2rad = function(deg) return(deg*pi/180)
# 
# gcd_hf(-111.93, 60.02, -133.72, 68.36)
# 
# gcd = function(lon1, lat1, lon2, lat2) {
#   # Convert degrees to radians
#   lon1 = deg2rad(lon1)
#   lat1 = deg2rad(lat1)
#   lon2 = deg2rad(lon2)
#   lat2 = deg2rad(lat2)
#   return(gcd_hf(lon1, lat1, lon2, lat2))
# }

# gcd(-111.93, 60.02, -133.72, 68.36)

# proj4string(df) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


# library(sp)
# data(meuse)
# coordinates(meuse) <- ~x+y
# m <- vgm(.59, "Sph", 874, .04)
# # five-fold cross validation:
# x <- krige.cv(log(zinc)~1, meuse, m, nmax = 40, nfold=155)
# bubble(x, "residual", main = "log(zinc): 5-fold CV residuals")

# meuse.g <- gstat(id = "zn", formula = log(zinc) ~ 1, data = meuse)
# meuse.g <- gstat(meuse.g, "cu", log(copper) ~ 1, meuse)
# meuse.g <- gstat(meuse.g, model = vgm(1, "Sph", 900, 1), fill.all = TRUE)
# x <- variogram(meuse.g, cutoff = 1000)
# meuse.fit = fit.lmc(x, meuse.g)
# out = gstat.cv(meuse.fit, nmax = 40, nfold = 5) 
# summary(out)
# out = gstat.cv(meuse.fit, nmax = 40, nfold = c(rep(1,100), rep(2,55))) 
# summary(out)
# # mean error, ideally 0:
# mean(out$residual)
# # MSPE, ideally small
# mean(out$residual^2)
# # Mean square normalized error, ideally close to 1
# mean(out$zscore^2)
# # correlation observed and predicted, ideally 1
# cor(out$observed, out$observed - out$residual)
# # correlation predicted and residual, ideally 0
# cor(out$observed - out$residual, out$residual)


# library(lme4)
# lmselected = selected_data %>% mutate(Station = factor(Station))
# 
# summary(lmer(Corr_Count ~ Pressure + Alt + (1|Date), lmselected))
# 
# 
# filtered_data = selected_data %>% na.omit %>% filter(Corr_Count != 0)
# 
# summary(lm(Corr_Count ~ Pressure + Alt + precipIntensity + dewPoint + humidity + visibility + cloudCover + temperature + Date, filtered_data))
# 
