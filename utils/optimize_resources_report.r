#!/usr/bin/env Rscript
#
# This script is made to create a summary of a job-output folder by creating graphs and a table.
# Graphs and table are used to compare the amount of resources asked and the amount used.
#
# IN : job_output folder path (str)
# OUT : folder containing graphs (.pdf) and a table containing all values (.csv)
#

#install.packages("dplyr")    
library(dplyr, warn.conflicts = FALSE) #avoid conflict message for duplicate functions
#library(lubridate)

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

#look to every folder and keep .o file path
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

	#creating df for parsed informations
	Info_df <<- data.frame(#JobId = as.integer(),
							#EligibleTime = as.character(),
		  					#StartTime = as.character(),
							WaitingTime = as.character(),
							RunTime = as.character(),
							TimeLimit = as.character(),
							NumCPUs = as.integer())


	for (i in file_path_list){
		for (j in i){
		  #file reading
		  FileInput = readLines(j)
		  
		  # Replace line feeds with spaces
		  FileInput_Space <- gsub(pattern = "\\n", replacement = " ", x = FileInput)
		  
		  #split file to extract value
		  FileInput_List <- strsplit(x = FileInput_Space, split = " ")
		  
		  #Create (global) temporary dataframe 
		  Info_df_temp <<- Info_df

		  #research JobId
		  #JobId <- research_Element(FileInput_List, "JobId")

		  #research EligibleTime
		  EligibleTime <- String_to_Date(research_Element(FileInput_List, "EligibleTime"))
		  
		  #research StartTime
		  StartTime <- String_to_Date(research_Element(FileInput_List, "StartTime"))

		  #calculation WaitingTime
		  WaitingTime <- (difftime(StartTime, EligibleTime, units = "mins"))

		  #research RunTime
		  RunTime <- research_Element(FileInput_List, "RunTime")

		  #research TimeLimit
		  TimeLimit <- research_Element(FileInput_List, "TimeLimit")

		  #research NumCPUs
		  NumCPUs <- research_Element(FileInput_List, "NumCPUs")

		  #add informations in Info_df_temp
		  Info_df_temp[nrow(Info_df_temp) + 1,] = list(WaitingTime, RunTime, TimeLimit, NumCPUs) 
		  
		}
		#compute maximum values if there is more than one .o file per folder
		if (nrow(Info_df_temp) > 1){
			for (var in colnames(Info_df_temp)){
				
				Info_df_temp[var] <- max_hour(Info_df_temp[var]) 
			}
		}
		Info_df_temp <- Info_df_temp %>% filter(row_number()==1)
		
		Info_df <- rbind(Info_df, Info_df_temp)

	}



	#Change WaitingTime format
	Info_df$WaitingTime <- round(as.double(Info_df$WaitingTime), digits = 2) %>% 
	 						gsub(pattern = "\\.", replacement = ":")
	for (i in 1:length(Info_df$WaitingTime)){
		time <- as.character(Info_df$WaitingTime[i])
		min <- strsplit(x = time, split = ":")[[1]][1]
		sec <- strsplit(x = time, split = ":")[[1]][2]

		# transforme minutes into hours and min
		h <- as.numeric(min) %/% 60
		#return(typeof(h))
		#h <- add_0_time(h)

		min <- as.numeric(min) %% 60
		#return(typeof(min))
		#return(floor(log10(h)) + 1)
		min <- add_0_time(min)

		#re-create the full WaintingTime value
		Info_df$WaitingTime[i] <- paste(c(h, min, sec), collapse=":")

	}

	

	#Convert day into hours TimeLimit
	for (i in 1:length(Info_df$TimeLimit)){

		time <- as.character(Info_df$TimeLimit[i])
		day_in_hours <- as.integer(strsplit(x = time, split = "-")[[1]][1]) * 24
		hminsec <- strsplit(x = time, split = "-")[[1]][2]
		h <- strsplit(x = hminsec, split = ":")[[1]][1]
		min <- strsplit(x = hminsec, split = ":")[[1]][2]
		sec <- strsplit(x = hminsec, split = ":")[[1]][3]
		day_and_hours <- as.numeric(day_in_hours) + as.numeric(h)

		# re-create the full TimeLimit value
		Info_df$TimeLimit[i] <- paste(c(day_and_hours, min, sec), collapse=":")
	}

	#return complete dataframe
	return (Info_df)
}

add_0_time <- function(number){
		if (length(number) == 1){
			number <- paste(c(0,number), collapse="")
		}
		return(number)
}

research_Element <- function(Input_file, Researched_element){
	#IN : pre-treated file : \n replaced by space then space splited into list, researched element (character)
	#OUT : wanted value associated with researched_element given (character)
	out <- tryCatch(
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
			### MESSAGE A AJOUTER POUR CAUSE STOP ###
		}
	)
	return (out)
}

String_to_Date <- function(stringD){
	#IN : date as character type
	#OUT : date as Date type
	only_date_String <- strsplit(x = stringD, split = "T")[[1]][2] #keep hour:min:sec
	only_date_String <- hms::as_hms(only_date_String)			   #change type to hms

	return(only_date_String)
}

max_hour <- function(df){
	#IN : list column containing Date type values
	#OUT : one Date type value
	return (max(unlist(df)))
}



#############
### TO DO ###

#memory

#############

system2("ls", stdout = TRUE, stderr = TRUE)


#function call with complete list of .o file
result <- parsed_folder(o_files_list)
print(result)

