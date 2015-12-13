# This function aggregates the neutron count data by choice of time for each station
# in the data set then calculates the means for correlated and uncorrelated counts 
# and pressure. The user can choose the aggregation unit from choices "Day" (default),
# "Month", or "Year".
agg_neutron_count_data = function(neutron_count_data, agg_unit = "Day") {
  if (agg_unit == "Day") {
    grouping_pattern = c("Year", "Month", "Day")
  } else if (agg_unit == "Month") {
    grouping_pattern = c("Year", "Month")
  } else {
    grouping_pattern = "Year"
  }
  grouping_pattern = c(grouping_pattern, "Station")
  dots = lapply(grouping_pattern, as.symbol)
  neutron_count_data %>% group_by_(.dots = dots) %>% summarize_each(funs(mean), Corr_Count, Pressure, Uncorr_Count)
}
# Testing
# neutron_counts_agg = agg_neutron_count_data(neutron_counts_tidy)
# neutron_counts_agg = agg_neutron_count_data(neutron_counts_tidy, by = "Month")
# neutron_counts_agg = agg_neutron_count_data(neutron_counts_tidy, by = "Year")
