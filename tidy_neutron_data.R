library(tidyr)
# Turns messy webscraped Neutron data into "tidy" format.
tidy_neutron_data = function(this_neutron_data) {
  corr_data = this_neutron_data %>% 
    select(Year:Minute, starts_with("Corr")) %>% 
    gather("Station", "Corr_Count", 6:9) %>% 
    mutate(Station = gsub("Corr_", "", Station))
  press_data = this_neutron_data %>% 
    select(Year:Minute, starts_with("Press")) %>%
    gather("Station", "Pressure", 6:9) %>%
    mutate(Station = gsub("Press_", "", Station))
  uncorr_data = this_neutron_data %>% 
    select(Year:Minute, starts_with("Uncorr")) %>%
    gather("Station", "Uncorr_Count", 6:9) %>%
    mutate(Station = gsub("Uncorr_", "", Station))
  suppressMessages(left_join(left_join(corr_data, press_data), uncorr_data))
}
# Testing
# test1 = tidy_neutron_data(this_data)