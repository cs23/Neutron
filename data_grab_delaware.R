library(readr)
library(stringr)
library(dplyr)

data_grab_delaware = function(file_location) {
  main_data = read_lines(file_location)
  parse_limits = grep("\\*", main_data)
  header = main_data[(parse_limits[1]-2):(parse_limits[1]-1)]
  main_data = main_data[(parse_limits[1]+1):(parse_limits[2]-1)]
  main_data = paste0(main_data, "\n")
  header[1] = gsub("South Pole", "South_Pole", header[1])
  header[1] = gsub("Fort Smith", "Fort_Smith", header[1])
  header_top_split = unlist(str_split(header[1], " "))
  header_top_split = header_top_split[-which(header_top_split == "")]
  header_top_split = header_top_split[-(1:3)]
  header_bottom_split = unlist(str_split(header[2], " "))
  header_bottom_split = header_bottom_split[-which(header_bottom_split == "")]
  header_bottom_split[1:5] = c("Year", "Month", "Day", "Hour", "Minute")
  index = grep("Corr", header_bottom_split)
  for (i in 1:4) {
    header_bottom_split[index[i]] = paste0(header_bottom_split[index[i]], "_", header_top_split[i])
    header_bottom_split[index[i]+1] = paste0(header_bottom_split[index[i]+1], "_", header_top_split[i])
    header_bottom_split[index[i]+2] = paste0(header_bottom_split[index[i]+2], "_", header_top_split[i])
  }
  main_data = tbl_df(read.table(textConnection(main_data)))
  colnames(main_data) = header_bottom_split
  return(main_data)
}