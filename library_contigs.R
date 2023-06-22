################################################################################
###Read excel files and create a gap file to blast. 
################################################################################


  setwd("/home/pallavi/git_repositories/M_hapla/Marker_oldtonew/Thomas_Markers_wo_AFLP/")
  workwd <- getwd()

library(tidyverse)
library(dplyr)
  
## This function reads the excel file and returns the modified excel files on the list. 

  ## read_excel_files will be the name of the function
  ## list.files will list all the files present in the directory that matches the pattern with xlsx extension. 
  ## full.names=TRUE will include the full path of each file. 
  ## data is an empty list to store the data from each excel file. 
  ## for loop loops over each excel file sheet. 
  ## read_excel function reads the data from the current sheet of the current file and saves it under sheet_data
  ## data is modified to include 250 bp above and below the markers.
  ## modified data is assigned to data using the current file's full path as key. 
  ## last function returns data list containg the modified data from all the excel files in the directory. 
  
read_excel_files <- function(directory) {
  files <- list.files(directory, pattern = "\\.xlsx$", full.names = TRUE)
  data <- list()
  
  for (file in files) {
    sheet_data <- readxl::read_excel(file)
    
    modified_data <- sheet_data %>%
      select(-c(4,5,6,7)) %>%                                           ### bed files only need the coordinates
      dplyr::mutate(start_use = Start - 250, end_use = Stop + 250) %>%
      select(-c(2,3))
    
    
    data[[file]] <- modified_data
  }
  
  return(data)
}  
  ### Calling the function to do its job:
directory <- workwd 

result <- read_excel_files(directory)

result[[1]] ##checking if everything is okay. 

## result is a list of excel files. We need to save it in a text format if possible. 
## in for loop we get the modified data from current file using result[[file_path]]
## create output file by appending _modified.txt to the original file path. 
## write.delim used to write the modified data to the output file and delimiter is set to tab-separated values. 

  for (file_path in names(result)) {
    modified_data <- result[[file_path]]
    output_file <- paste0(file_path, "_modified.bed")
    write_delim(modified_data, file = output_file, delim = "\t", col_names = FALSE)
  }
    
    
    
    
    
    