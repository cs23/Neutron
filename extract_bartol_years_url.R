# Extracts the individual urls according to years chosen from main Bartol Table URL.
extract_bartol_year_url = function(bartol_url, years = 2014) {
  stopifnot(is.numeric(years))
  years = years[order(years)]
  start_end_limits = c(1957, 2016)
  if (years[1] < start_end_limits[1]) stop("Choose a starting year after 1956!")
  if (years[length(years)] > start_end_limits[2]) stop("Choose an ending year before 2017!")
  lines = readLines(bartol_url)
  these_years = paste0(years, collapse = "|")
  all_urls = c(clean_url_lines(lines), clean_url_lines(lines, new_sites = TRUE))
  these_years = grep(these_years, all_urls, value = TRUE)
  these_years = these_years[order(these_years)]
  return(these_years)
}
# Testing
# bartol_url = "http://neutronm.bartol.udel.edu/~pyle/bri_table.html"
# extract_bartol_year_url(bartol_url, 2012:2014) # Works!
# extract_bartol_year_url(bartol_url, 2012:2018) # Fails correctly!
# extract_bartol_year_url(bartol_url, 2012:2010) # Works!
# extract_bartol_year_url(bartol_url, 1953:2014) # Fails correctly!
# extract_bartol_year_url(bartol_url, 2014:1953) # Fails correctly!
# extract_bartol_year_url(bartol_url, 1957) # Works correctly!
# extract_bartol_year_url(bartol_url, 2016) # Works correctly!

# Helper function for extract_bartol_year_url that cleans URL links.
clean_url_lines = function(dirty_url_lines, new_sites = FALSE) {
  target_pattern = ifelse(new_sites, ".*ftp.*/BRIData/BRI\\d{4}B.txt*", ".*ftp.*/BRIData/BRI\\d{4}.txt*")
  neutron_count_lines = grep(pattern = target_pattern, x = dirty_url_lines, value=TRUE) %>% 
    sub(pattern = ".*href=.", replacement = "") %>% 
    sub(pattern = ".txt.*", replacement = ".txt")
}
# For reference in clean_url_lines
# Data for McMurdo, Swarthmore/Newark, South Pole, and Thule see ftp://ftp.bartol.udel.edu/pyle/BRIData/BRI2014.txt
# Data for Fort Smith, Peawanuck, Nain and Inuvik see ftp://ftp.bartol.udel.edu/pyle/BRIData/BRI2014B.txt