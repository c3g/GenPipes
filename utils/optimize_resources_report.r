#!/usr/bin/env Rscript
#
# This script is made to create a summary of a job-output folder by creating graphs and a table.
# Graphs and table are used to compare the amount of resources asked and the amount used.
#
# IN : job_output folder path (str)
# OUT : folder containing graphs (.pdf) and a table containing all values (.csv)
#

# ask for job_output path
#cat("Job_output path: ");
#job_output_path <- readLines("stdin",n=1);

#à modifier, permet de pas avoir à rentrer le nom du fichier à chaque fois
job_output_path = "job_output"


#list all folders path in job_output folder as list
job_output_list_folder <- as.list(list.dirs(path = job_output_path,
									recursive = FALSE,
									full.names = TRUE))

#Same parsing and values compilation for each folder
o_files_list = list() 


#Create list of .o file
	# sub list for each folder

for (folder_path in job_output_list_folder) {
	# open folder and parse .o files 

	o_files_list <- append(o_files_list, list(list.files(folder_path,
								pattern = "\\.o$",
								recursive = TRUE,
								full.names = TRUE)))
}

parsed_folder <- function(file_path_list){
	#IN : list of .o file path
	#OUT : dataset containing wanted information from these files
	
}