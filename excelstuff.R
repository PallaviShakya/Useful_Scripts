## Merge to excel files on the basis of a column. 

## To do: 
### In Thomas marker excel files, merge columns 1 and two with a dot in the middle.
### Merge all the exccel files into one big excel file with sheet name as lg1, lg2.. etc. 
### Merge the merged column with another excel file where there are 19 sheets based on the sheet number?

library(readxl)
library(openxlsx)
library(fs)


merge_excel_columns<- function(input_dir, output_file) {
  
  #Getting the list from the directory
  
  excel_files <- list.files("/home/pallavi/git_repositories/M_hapla/Marker_oldtonew/Thomas_Markers/", pattern = "\\.xlsx$", full.names = TRUE)
  
  # Initialize empty directory
  
  merged_data <- data.frame()
  
  #Loop
  
  for (file in excel_files) {
    
    data <- read_excel(file)
    
    # merge contigs with start wit . in the middle
    data <- data %>% 
      unite(new_column, Contig, Start, sep = ".")
    
    merged_data <- bind_rows(merged_data, data)
  }
  
  ## Write to a new excel file 
  
  write.xlsx(merged_data, output_file)
  
    
}

input_dir = "/home/pallavi/git_repositories/M_hapla/Marker_oldtonew/Thomas_Markers/"
output_file = "/home/pallavi/git_repositories/M_hapla/Marker_oldtonew/"

merge_excel_columns(input_dir, output_file)
