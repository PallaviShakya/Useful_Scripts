# Read the geneticmap excel file

# Excel file info: 
### Column1: Scaffold number,
### Column2: The first position of the BLAST window
### Column3: First position +250 (because this is where the Marker starts)
### Column4: Position where marker ends
### Column5: Linkage group
### Cloumn6: Marker type

## To do:
# Remove column 2 and 7
# merge column 5 and 6 with : in the middle. 

setwd("/home/pallavi/git_repositories/M_hapla/Marker_oldtonew/ALLMAPs_input_map/")
workwd <- getwd()

library(tidyverse)
library(dplyr)
library(readxl)
library(openxlsx)


genetic_map <- read_excel("geneticmap_input_woAFLP.xlsx", col_names = FALSE)

genetic_map1 <- genetic_map[, -c(2,7)]


genetic_map1$distance <- paste0(genetic_map$...1, ":", genetic_map$...6)
genetic_map1$new_column <- "GeneticMap"
genetic_map1$dummy <- paste0(genetic_map1$new_column, "-", genetic_map1$distance)

genetic_map_final <- genetic_map1[, -c(4,5,6,7)]


write_delim(genetic_map_final, file = "final_genetic_map_thomas_wo_AFLP.bed", delim = "\t", col_names = FALSE)

write_tsv(genetic_map_final, file = "final_genetic_map_thomas_wo_AFLP.tsv", col_names = FALSE)

####################################################################################
separate_excel_by_column <- function(input_file, column_name, output_file) {
  # Read the input Excel file
  excel_data <- read_excel(input_file)
  
  # Get unique values in the specified column
  unique_values <- unique(excel_data[[column_name]])
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Iterate over unique values and create a sheet for each group
  for (value in unique_values) {
    group_data <- excel_data[excel_data[[column_name]] == value, ]
    
    ## Remove columns we don't need
    
    group_data <- group_data[, -c(2, 4, 7)]
    
    addWorksheet(wb, sheetName = as.character(value))
    writeData(wb, sheet = value, group_data)
  }
  
  # Save the workbook to the output file
  saveWorkbook(wb, output_file)
  
  cat("Excel file separated and saved successfully!")
}

input_file <- "geneticmap_input.xlsx"
column_name <- "LG"
output_file <- "VW8XVW9.xlsx"

separate_excel_by_column(input_file, column_name, output_file)

