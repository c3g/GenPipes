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
job_output_path = "~/Documents/local/projet/optimize_resources_report/job_output"


#list all folders path in job_output folder as list
job_output_list_folder <- (list.dirs(path = job_output_path,
									recursive = FALSE,
									full.names = TRUE))

#list containing sub list for each folder
o_files_list = list() 

for (i in 1:length(job_output_list_folder)) {
	# open each folder from folder list 
  folder_path <- job_output_list_folder[i]
  
  # find all .o file in each sub folder
  new_folder <- list(list.files(folder_path,
                           pattern = "\\.o$",
                           recursive = TRUE,
                           full.names = TRUE))
  
  # some folder doesn't have any .o file, condition to avoid type error
  if (length(new_folder) == 0 ){
    o_files_list[i] <- NULL
  }else{
    o_files_list[i] <- new_folder
  }
}


parsed_folder <- function(file_path_list){
	#IN : list of .o file path
	#OUT : dataset containing wanted information from these files
	for (i in file_path_list){
		for (j in i){
		  #file reading
		  FileInput = readLines(j)
		  
		  # Replace line feeds with spaces
		  FileInput_Space <- gsub(pattern = "\\n", replacement = " ", x = FileInput)

		  #cut before FAKE EPILOGUE
		  #EndFileInput_Space <- strsplit(x = FileInput_Space, split = "SLURM FAKE EPILOGUE (MUGQIC)")
		  
		  #split file to extract value
		  FileInput_List <- strsplit(x = FileInput_Space, split = " ")#[[1]]
		  

		  #research EligibleTime
		  EligibleTime <- research_Element(FileInput_List, "EligibleTime")
		  
		  #research StartTime
		  StartTime <- research_Element(FileInput_List, "StartTime")

		  #research RunTime
		  RunTime <- research_Element(FileInput_List, "RunTime")

		  #research TimeLimit
		  TimeLimit <- research_Element(FileInput_List, "TimeLimit")

		  #research NumCPUs
		  NumCPUs <- research_Element(FileInput_List, "NumCPUs")
		  
		  return (typeof(FileInput_List))
		}
	}
}


research_Element <- function(Input_file, Researched_element){
	#IN : pre-treated file : \n replaced by space then space splited into list, researched element (string)
	#OUT : 
	out <- tryCatch
		{
			#research element StartTime
			Pos_eli_time <- grep(Researched_element,Input_file)
			Pos1 <- as.integer(strsplit(x = as.character(Pos_eli_time), split = " ")[[2]]) #position of list containing Start Time

			#list containing research element
			list <- Input_file[Pos1]

			#position of element inside sub list
			Pos2 <- as.integer(grep(Researched_element,list[[1]]))
			Complet_Element <- Input_file[[Pos1]][Pos2]

			Only_Element <- strsplit(x = Complet_Element, split = "=")[[1]][2]

			return (Only_Element)
		},
		error = function(cond){

		}
	return (out)
}

#execution time

#memory


#function call with complete list of .o file
result <- parsed_folder(o_files_list)
print(result)
