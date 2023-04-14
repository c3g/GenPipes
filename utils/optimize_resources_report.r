#!/usr/bin/env Rscript
#
# This script is made to create a summary of a job-output folder by creating graphs and a table.
# Graphs and table are used to compare the amount of resources asked and the amount used.
#
# IN : job_output folder path (str)
# OUT : folder containing graphs (.pdf) and a table containing all values (.csv)
#

library(dplyr, warn.conflicts = FALSE)  #avoid conflict message for duplicate functions
library(kimisc)                         #for seconds_to_hms function
library(lubridate)

# Library for plot
library(ggplot2)

#library for RMarkDown
library(knitr)
library(markdown)


# ask for job_output path
#cat("Job_output path: ");
#job_output_path <- readLines("stdin",n=1);

job_output_path = "~/Documents/local/projet/optimize_resources_report/job_output" 	#à modifier, permet de pas avoir à rentrer le nom du fichier à chaque fois
#job_output_path = "/scratch/matteol/genpipes_test/job_output"

folder_path_to_o_file_list <- function(job_output_path){
	#IN : folder path as character
	#OUT : list containing sub list for each step type containing .o file path

	#list all folders path in job_output folder as list
	job_output_list_folder <- (list.dirs(path = job_output_path,
										recursive = FALSE,
										full.names = TRUE))
	#folder path list
	list_path = list()

	for (path in job_output_list_folder){

		list_path <- append(list_path,list.files(path,
	                           pattern = "\\.o$",
	                           recursive = TRUE,
	                           full.names = TRUE))
	}

	#dataframe containing path and associated type of step
	df_path <<- data.frame(	path = as.character(),
								file_step_name = as.character())

	#create a dataframe containing all path as keys and associated name step as values
	for (path in list_path){
		file_step_name <- sapply(strsplit(x = path, split = "/"), tail, 1) %>% 	#element after last "/"
							strsplit(split = "\\.") %>%							#element before first "."
								sapply(head,1)

		#add new line with informations in dataframe
		df_path[nrow(df_path) + 1,] = list(path, file_step_name)
	}

	#list of steps
	unique_step_name <- unique(df_path$file_step_name)

	#list containing sub list for each folder
	o_files_list = list()

	#for each step, a list is created and add to a global list
	for (step in unique_step_name){

		#new list containing all path 
		new_list <- list(df_path$path[which(df_path$file_step_name == step)])

		o_files_list <- append(o_files_list,new_list)

	}
	return (o_files_list)
}

parsed_folder <- function(folder_path_list){
	#IN : list of .o file path
	#OUT : dataset containing wanted information from these files

	#récup path list
	file_path_list <- folder_path_to_o_file_list(folder_path_list)

	#creating df for parsed informations
	Info_df <<- data.frame(	JobName = as.character(),
							WaitingTime = as.character(),
							RunTime = as.character(),
							TimeLimit = as.character(),
							NumCPUs = as.integer(),
							Memory_Efficiency = as.numeric(),
							Memory_Request = as.numeric())
  
	DF_plot <<- data.frame(	JobName = as.character(),
	                        WaitingTime = as.character(),
	                        RunTime = as.character(),
	                        TimeLimit = as.character(),
	                        NumCPUs = as.integer(),
	                        Memory_Efficiency = as.numeric(),
	                        Memory_Request = as.numeric())
	
	for (i in file_path_list){

		#Create temporary dataframe 
		Info_df_temp <- data.frame(	JobName = as.character(),
									WaitingTime = as.character(),
									RunTime = as.character(),
									TimeLimit = as.character(),
									NumCPUs = as.integer(),
									Memory_Efficiency = as.numeric(),
									Memory_Request = as.numeric())

		for (j in i){

			#file reading
			FileInput = readLines(j)

			# Replace line feeds with spaces
			FileInput_Space <- gsub(pattern = "\\n", replacement = " ", x = FileInput)

			#split file to extract value
			FileInput_List <- strsplit(x = FileInput_Space, split = " ")

			#research StepName
			JobName <- research_Element(FileInput_List, "JobName")
			JobName <- strsplit(x = JobName, split = "\\.")[[1]][1]

			#research EligibleTime
			EligibleTime <- String_to_Date(research_Element(FileInput_List, "EligibleTime"))

			#research StartTime
			StartTime <- String_to_Date(research_Element(FileInput_List, "StartTime"))

			#calculation WaitingTime
			WaitingTime <- as.character(difftime(StartTime, EligibleTime, units = "mins"))

			#research RunTime
			RunTime <- research_Element(FileInput_List, "RunTime")
      
			#research TimeLimit
			TimeLimit <- (research_Element(FileInput_List, "TimeLimit"))

			#research NumCPUs
			NumCPUs <- research_Element(FileInput_List, "NumCPUs")

			#seff for memory request and
			#research JobId
			JobId <- research_Element(FileInput_List, "JobId")

			#seff command give memory efficiency information (and more)
			#Memory_Efficiency
			# seff_resp <- system2("seff", args = JobId, stdout = TRUE)
			# Pos_eli_time <- grep("Memory Efficiency",seff_resp)
			# #return(Pos_eli_time)
			# print(seff_resp)
			# print(Pos_eli_time)
			# Memory_Efficiency <- strsplit(x = seff_resp[Pos_eli_time], split = " ")[[1]][3]
			# Memory_Efficiency <- strsplit(x = Memory_Efficiency, split = "%")[[1]][1]

			# Memory_Efficiency <- as.numeric(as.character(Memory_Efficiency))

			# # #Memory_Request
			# Memory_Request <- strsplit(x = Pos_eli_time, split = " ")[[1]][5]

			Memory_Efficiency <- NA
			Memory_Request <- NA

			#fill df_info with new informations and rename columns
			df_info <- data.frame(JobName, WaitingTime, RunTime, TimeLimit, NumCPUs, Memory_Efficiency, Memory_Request)
			names(df_info) <- c("JobName", "WaitingTime", "RunTime", "TimeLimit", "NumCPUs", "Memory_Efficiency", "Memory_Request")

			Info_df_temp <- rbind(Info_df_temp, df_info)

		}

		#change TimeLimit format
		Info_df_temp$TimeLimit <- day_into_hours(Info_df_temp$TimeLimit)
		
		#Specfific dataframe for ploting
		DF_plot <- rbind(DF_plot, Info_df_temp)
		
		#compute average values if there is more than one .o file per folder
		if (nrow(Info_df_temp) > 1){

			#re-create the dataframe with the average values
			Info_df_temp <- data.frame(JobName,
									   mean_value(Info_df_temp$WaitingTime),
									   mean_value(Info_df_temp$RunTime),
									   mean_value(Info_df_temp$TimeLimit),
									   mean_value(as.integer(Info_df_temp$NumCPUs)),
									   mean_value(Info_df_temp$Memory_Efficiency),
									   mean_value(Info_df_temp$Memory_Request))

			#rename columns
			names(Info_df_temp) <- c("JobName", "WaitingTime", "RunTime", "TimeLimit", "NumCPUs", "Memory_Efficiency", "Memory_Request")
		}
		
		#put Info_df_temp at the end of the main dataframe containing all informations
		Info_df <- rbind(Info_df, Info_df_temp)
	
	}

	 #Change WaitingTime format
	 Info_df$WaitingTime <- ChangeWaitingTimeFormat(Info_df$WaitingTime)
	 DF_plot$WaitingTime <- ChangeWaitingTimeFormat(DF_plot$WaitingTime)
# 	 Info_df$WaitingTime <- round(as.double(Info_df$WaitingTime), digits = 2) %>% 
#  	 						gsub(pattern = "\\.", replacement = ":")
# 
# 	 for (i in 1:length(Info_df$WaitingTime)){
# 	 	time <- as.character(Info_df$WaitingTime[i])
# 	 	min <- strsplit(x = time, split = ":")[[1]][1]
# 	 	sec <- strsplit(x = time, split = ":")[[1]][2]
# 
# 	 	# transforme minutes into hours and min
# 	 	h <- as.numeric(min) %/% 60
# 	 	h <- add_0_time(h)
# 
# 	 	min <- as.numeric(min) %% 60
# 	 	min <- add_0_time(min)
# 
# 	 	#specific verification for sec because of
# 	 	if (is.na(sec)){
# 	 		sec <- "00"
# 	 	}else{
# 	 		sec <- as.numeric(sec)
# 	 		sec <- add_0_time(sec)
# 	 	}
# 
# 	 	#re-create the full WaitingTime value
# 	 	Info_df$WaitingTime[i] <- paste(c(h, min, sec), collapse=":")

	#}

	#return complete dataframe
	#return (Info_df, DF_plot)
	return(list(Info_df, DF_plot))
}

ChangeWaitingTimeFormat <- function(df){
  #IN : WaitingTime column as dataframe with random format
  #OUT : WaitingTime column with H:M:S format
  df <- round(as.double(df), digits = 2) %>% 
    gsub(pattern = "\\.", replacement = ":")
  
  for (i in 1:length(df)){
    
    time <- as.character(df[i])
    min <- strsplit(x = time, split = ":")[[1]][1]
    sec <- strsplit(x = time, split = ":")[[1]][2]
    
    # transforme minutes into hours and min
    h <- as.numeric(min) %/% 60
    h <- add_0_time(h)
    
    min <- as.numeric(min) %% 60
    min <- add_0_time(min)
    
    #specific verification for sec because of
    if (is.na(sec)){
      sec <- "00"
    }else{
      sec <- as.numeric(sec)
      sec <- add_0_time(sec)
    }
    
    #re-create the full WaitingTime value
    df[i] <- paste(c(h, min, sec), collapse=":")
  }
  return(df)
}

day_into_hours <- function(df){
	#Convert day into hours TimeLimit
	for (i in 1:length(df)){

		time <- as.character(df[i])
		day_in_hours <- as.integer(strsplit(x = time, split = "-")[[1]][1]) * 24
		
		#specific case where there is less than a day of TimeLimit
		if (is.na(day_in_hours)){
			day_in_hours <- 0 						#less than one day
			hminsec <- time
		}else{
			hminsec <- strsplit(x = time, split = "-")[[1]][2] 		#one day or more
		}

		h <- strsplit(x = hminsec, split = ":")[[1]][1]
		min <- strsplit(x = hminsec, split = ":")[[1]][2]
		sec <- strsplit(x = hminsec, split = ":")[[1]][3]
		day_and_hours <- as.numeric(day_in_hours) + as.numeric(h)
		day_and_hours <- add_0_time(day_and_hours)

		# re-create the full TimeLimit value
		df[i] <- paste(c(day_and_hours, min, sec), collapse=":")
	}
	return(df)
}

add_0_time <- function(number){
		#print("number")
		#print(number)
		if (is.na(number)){
		  number <- 0
		}

		if (floor(log10(number)) == 0 | number == 0){
			number <- paste(c(0,number), collapse="")
		}

		return(number)
}

research_Element <- function(Input_file, Researched_element){
	#IN : pre-treated file : \n replaced by space then space splited into list, researched element (character)
	#OUT : wanted value associated with researched_element given (character)
	out <- tryCatch(
		{
			#research element 
			Pos_eli_time <- grep(Researched_element,Input_file)
			Pos1 <- as.integer(strsplit(x = as.character(Pos_eli_time), split = " ")[[2]]) #position of list containing element

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

mean_value <- function(df){
	#IN : dataframe column containing numerical values or h:m:s values with character type
	#OUT : mean value of the dataframe column with character type

	#case where df contain character values
	if (is.character(df)){
	  #print(df)
	  #print(typeof(df))
		df <- df %>%
				hms()%>%
				period_to_seconds()%>%
				mean()%>%
		    seconds_to_period()%>%
		    seconds.to.hms()

	# df already contain numerical values
	}else{
		df <- mean(df)
	}
	
	return (as.character(df))
}

Keep_hour <- function(df, df2 = NA, unite = NA){
  #Keep only hours or minutes or seconds depending on the value
  #IN : dataframe containing H:M:S values
  #OUT : dataframe containing same hour value with different format depending on the value.
  #     df2 format will be the same as df format and decision based on df values
  
  #at least one hour ?
  if (!hour(hms(max(df))) == 0){
    df <- hour(hms(df))
    if (!is.na(max(df2))) { df2 <- hour(hms(df2)) }
    unite <- "hour"
  }
  #at least one minute ?
  else if (!minute(hms(max(df))) == 0){
    df <- hour(hms(df))*60 +  minute(hms(df))
    if (!is.na(max(df2))) { df2 <- hour(hms(df2))*60 +  minute(hms(df2)) }
    unite <- "minute"
  }
  #values are in seconds or value equal 0
  else{
    df <- hour(hms(df))*3600 +  minute(hms(df))*60 + second(hms(df))
    if (!is.na(max(df2))) { df2 <- hour(hms(df2))*3600 +  minute(hms(df2))*60 + second(hms(df2)) }
    unite <- "second"
  }
  
  if (is.na(max(df2))){
    return(list(df, unite))
  }
  
  return(list(df, df2, unite))

}


############# TO DO
### TO DO ###

#graph memory
  #verif modules beluga
  #test beluga
#création RMarkdown

#memory avec seff en dessous, à lancer depuis beluga donc push depuis local et pull depuis beluga
	# ga nom de fichier --> add
	# gc -m "message"
	# gp

	#gst pour avoir state et vérifier état

	# push faisable depuis sublime merge = plus facile

	# gl pour pull depuis ~/apps/genpipes pour tout à jour

############## RESULT ###########################################################

#function call with complete list of .o file
result <- parsed_folder(job_output_path)

Info_df <- result[[1]]
DF_plot <- result[[2]]

#' The following dataframe contain average values parsed in .o files contained in the job_output file gived.
print(Info_df)

#' The following dataframe contain all values found and will be used to create next plots
print(DF_plot)

############## Register file ###################################################
##Register dataframe as CSV
#Create name with date
actual_date_time <- strsplit(x= as.character(Sys.time()), split = " ")[[1]][1:2] %>%
                paste(collapse ="_")
file_name <- paste(c("job_output_analyse", actual_date_time), collapse ="_")
complete_path_name <- paste(c(job_output_path,file_name), collapse ="/")
complete_path_name_csv <- paste(c(complete_path_name,"csv"), collapse =".")

#Write csv in job_output folder
write.csv(Info_df, complete_path_name_csv, row.names=TRUE)

 
############## Custom specfific dataFrame for ploting ##########################
#Change WaitingTime, RunTime and TimeLimit format

Waiting_unite <- Keep_hour(DF_plot$WaitingTime)[[2]]
DF_plot$WaitingTime <- Keep_hour(DF_plot$WaitingTime)[[1]]

Run_Limit <- Keep_hour(DF_plot$RunTime, DF_plot$TimeLimit)

DF_plot$RunTime <- Run_Limit[[1]]
DF_plot$TimeLimit <- Run_Limit[[2]]
Time_unite <- Run_Limit[[3]]

#computing percentage amount for RunTime compared to TimeLimit
DF_plot$RunTime_Efficiency <- round((DF_plot$RunTime / DF_plot$TimeLimit) *100, 3)

#compute maximal Efficiency for each step
DF_max_Eff <- DF_plot[c(1, 8)] %>% group_by(JobName) %>% top_n(1, RunTime_Efficiency)

#Merge DF_max_Eff with DF_plot
DF_plot <- merge(DF_plot, DF_max_Eff, by = "JobName")

############## WaitingTime Plot ################################################

p_WaintingTime <- ggplot(DF_plot, aes(x=as.factor(JobName), y=WaitingTime, label= as.numeric(WaitingTime))) + 
              theme(panel.background = element_rect(fill = 'white', color = 'grey'), 
                    panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) +
              geom_boxplot(fill="slateblue", alpha=0.2) + 
              coord_flip() +
              geom_text(hjust=-0.5, vjust=-0.5) +
              ylab(paste(c("WaitingTime (", Waiting_unite, ")"), collapse ="")) +
              xlab("Step name") +
              ggtitle("WaitingTime for each step (EligibleTime to StartTime)")

############## RunTime vs RunTime_Efficiency Plot ##############################


p_RunTime <- ggplot(DF_plot, aes(x=as.factor(JobName))) +
  
              theme(panel.background = element_rect(fill = 'white', color = 'grey'), 
                    panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),
                    axis.title.x.top = element_text(color = "red", size=13),
                    axis.title.x.bottom = element_text(color = "blue", size=13),
                    axis.title.y = element_text(size=13),
                    ) +
              
              coord_flip() +
              
              xlab("Step name") +
              
              ylab(paste(c("RunTime (", Time_unite, ")"), collapse ="")) +
              
              geom_boxplot( aes(y=RunTime),
                            alpha=0.2,
                            color="blue",
                            fill="#69b3a2",
                            ) + 
            
              geom_point( aes(y=RunTime_Efficiency.y),
                            color="red",
                            alpha=0.7) +
            
              ggtitle("RunTime and RunTime_Efficiency") +
              
              geom_text(y =  as.numeric(DF_plot$RunTime_Efficiency.y),
                        label = as.numeric(DF_plot$RunTime_Efficiency.y),
                        color="red",
                        size=3,
                        nudge_x = -0.5, nudge_y = -0.5,
                        check_overlap = TRUE) +
              
              # geom_label(
              #   y =  as.numeric(DF_plot$RunTime),
              #   label = as.numeric(DF_plot$RunTime),
              #   nudge_x = 0.5, nudge_y = 0.5
              # ) +
            
              # geom_label(
              #   #data = data.frame(DF_plot$RunTime_Efficiency),
              #   y =  DF_plot$RunTime_Efficiency,
              #   label = DF_plot$RunTime_Efficiency,
              #   nudge_x = 0.5, nudge_y = 0.5
              # ) +
            
              scale_y_continuous(
                
                # Features of the first axis
                name = paste(c("RunTime (", Time_unite, ")"), collapse =""),
                
                # Add a second axis and specify its features
                #sec.axis = sec_axis(~.*coeff, name="Second Axis")
                sec.axis = sec_axis(~./ max(DF_plot$RunTime),
                                    name="RunTime_Efficiency (percentage)")
              )




############## Memory_Efficiency vs Memory_Request #############################

p_Memory <- ggplot(DF_plot, aes(x=as.factor(JobName))) +
  
              theme(panel.background = element_rect(fill = 'white', color = 'grey'), 
                    panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),
                    axis.title.x.top = element_text(color = "red", size=13),
                    axis.title.x.bottom = element_text(color = "blue", size=13),
                    axis.title.y = element_text(size=13),
              ) +
              
              coord_flip() +
              
              xlab("Step name") +
              
              ylab(paste(c("Memory_Request (", Time_unite, ")"), collapse ="")) +
              
              geom_boxplot( aes(y=Memory_Request),
                            alpha=0.2,
                            color="blue",
                            fill="#69b3a2",
              ) + 
              
              geom_point( aes(y=Memory_Efficiency),
                          color="red",
                          alpha=0.7) +
              
              ggtitle("Memory_Request and RunTime_Efficiency") +
              
              geom_text(y =  as.numeric(DF_plot$Memory_Efficiency),
                        label = as.numeric(DF_plot$Memory_Efficiency),
                        color="red",
                        size=3,
                        nudge_x = -0.5, nudge_y = -0.5,
                        check_overlap = TRUE) +
              
              # geom_label(
              #   y =  as.numeric(DF_plot$Memory_Request),
              #   label = as.numeric(DF_plot$Memory_Request),
              #   nudge_x = 0.5, nudge_y = 0.5
              # ) +
              
              # geom_label(
              #   #data = data.frame(DF_plot$Memory_Efficiency),
              #   y =  DF_plot$Memory_Efficiency,
              #   label = DF_plot$Memory_Efficiency,
            #   nudge_x = 0.5, nudge_y = 0.5
            # ) +
            
            scale_y_continuous(
              
              # Features of the first axis
              name = paste(c("Memory_Request (", Time_unite, ")"), collapse =""),
              
              # Add a second axis and specify its features
              #sec.axis = sec_axis(~.*coeff, name="Second Axis")
              sec.axis = sec_axis(~./ max(DF_plot$Memory_Request),
                                  name="Memory_Efficiency (between 0 and 1)")
            )




############## R MarkDown #######################################################
#toutes les infos (csv + plots + rapides explications de ce qui est montré)
#peut-être intéractif, couleurs changent en fonction des valeurs des plots
# 
#rmarkdown::render("~/Documents/local/apps/genpipes/utils/optimize_resources_report.r")

# Actual_path <- rstudioapi::getActiveDocumentContext()$path #get actual absolute path
# rmarkdown::render(input = Actual_path,
#                   output_dir = "/Users/mleguen/Documents/local/apps/genpipes/utils/",
#                   clean = TRUE
#                   )







