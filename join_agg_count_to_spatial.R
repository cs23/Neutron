# This function joins the aggregated count data to the cleaned spatial data. The output
# data table is in a format which weather data can be gathered.
join_agg_count_to_spatial = function(agg_count_data, raw_spatial_data_location) {
  cleaned_spatial_data = clean_spatial_data(raw_spatial_data_location)
  joined_agg_data = left_join(agg_count_data, cleaned_spatial_data)
  return(joined_agg_data)
}
# Testing
# raw_spatial_data_location = "Delaware_Stations.csv"
# neutron_counts_agg_spatial = join_agg_count_to_spatial(neutron_counts_agg, raw_spatial_data_location)

# Helper function for join_agg_count_to_spatial that cleans spatial data by converting
# latitude and longitude to positive (North and East) and negative (South and West)
# along with separating station name from the state or country. Input is location of
# csv file.
clean_spatial_data = function(spatial_data_location) {
  locations = read_csv(spatial_data_location, col_types = "cccdc")
  locations = suppressWarnings(locations %>% 
                                 mutate(Lat = ifelse(grepl("S$", Lat), -as.numeric(gsub("S$", "", Lat)), as.numeric(gsub("N$", "", Lat))),
                                        Lon = ifelse(grepl("W$", Lon), -as.numeric(gsub("W$", "", Lon)), as.numeric(gsub("E$", "", Lon)))))
  locations %>% separate(Station, c("Station", "State/Country"), sep = ",")
}
# Testing
# cleaned_spatial_data = clean_spatial_data(raw_spatial_data_location)
