source("data_grab_delaware.R")
source("extract_bartol_years_url.R")
source("tidy_neutron_data.R")

# Gets neutron data from Bartol website for given years with nice tidy style formatting.
# This function is a main component of the data harvest pipeline
get_neutron_count_data = function(bartol_table_url, years = 2014){
  year_url = extract_bartol_year_url(bartol_table_url, years)
  all_data = tbl_df(data.frame())
  for (i in 1:length(year_url)) {
    this_data = data_grab_delaware(year_url[i]) %>% tidy_neutron_data
    all_data = bind_rows(all_data, this_data)
  }
  return(all_data)
}
# Testing
# bartol_url = "http://neutronm.bartol.udel.edu/~pyle/bri_table.html"
# data_test = get_neutron_count_data(bartol_url, 2012:2014) # Gets 2012, 2013, and 2014 data
# data_test = get_neutron_count_data(bartol_url) # Gets 2014 data only
# data_test = get_neutron_count_data(bartol_url, 1958) # Gets 1958 data only
