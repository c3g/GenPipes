#!/usr/bin/env Rscript
#
# This script is made to create a summary of a job-output folder by creating graphs and a table.
# Graphs and table are used to compare the amount of resources asked and the amount used.
#
# IN : job_output folder path (str)
# OUT : folder containing graphs (.pdf) and a table containing all values (.csv)
#
suppressMessages(library(dplyr, warn.conflicts = FALSE))  #avoid conflict message for duplicate functions
suppressMessages(library(kimisc))                         #for seconds_to_hms function
suppressMessages(library(lubridate))
suppressMessages(library(optparse))

# Library for plot
suppressMessages(library(ggplot2))

#library for RMarkDown
suppressMessages(library(knitr))
suppressMessages(library(markdown))


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
  #IN :   list of .o file path
  #OUT :  dataset containing wanted information from these files
  
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
  
  #data.frame creation for plots, this dataframe is made to keep all values
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
      seff_resp <- system2("seff", args = JobId, stdout = TRUE)
      Pos_eli_time <- grep("Memory Efficiency",seff_resp)
      # # # #return(Pos_eli_time)
      # # # print(seff_resp)
      # # # print(Pos_eli_time)
      Memory_Efficiency <- strsplit(x = seff_resp[Pos_eli_time], split = " ")[[1]][3]
      Memory_Efficiency <- strsplit(x = Memory_Efficiency, split = "%")[[1]][1]
      #  
      Memory_Efficiency <- as.numeric(as.character(Memory_Efficiency))
      #  
      # # #Memory_Request
      Memory_Request <- strsplit(x = seff_resp[Pos_eli_time], split = " ")[[1]][5]
      
      #Memory_Efficiency <- NA
      #Memory_Request <- NA
      
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
                                 mean(as.numeric(Info_df_temp$Memory_Request)))
      
      #rename columns
      names(Info_df_temp) <- c("JobName", "WaitingTime", "RunTime", "TimeLimit", "NumCPUs", "Memory_Efficiency", "Memory_Request")
    }
    
    #put Info_df_temp at the end of the main dataframe containing all informations
    Info_df <- rbind(Info_df, Info_df_temp)
    
  }
  
  #Change WaitingTime format
  Info_df$WaitingTime <- ChangeWaitingTimeFormat(Info_df$WaitingTime)
  DF_plot$WaitingTime <- ChangeWaitingTimeFormat(DF_plot$WaitingTime)
  
  
  #return complete dataframe
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
  #IN :   data.frame object containing D-H:M:S format values
  #OUT :  data.frame object containing H:M:S format values where days are turned into hours 
  
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
  #IN :   a number object
  #OUT :  same number with a 0 in front of it if it's smaller than 10
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
      print("The following element asn't been found : ")
      print(Researched_element)
      
      return(NA)
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
    df <- df %>%
      hms()%>%
      period_to_seconds()%>%
      mean() %>%
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



job_output_path = "~/Documents/local/projet/optimize_resources_report/job_output" 	#à modifier, permet de pas avoir à rentrer le nom du fichier à chaque fois
#job_output_path = "/scratch/matteol/genpipes_test/job_output"

##### OPTPARSE OPTIONS #####

option_list = list(
  make_option(c("-i", "--in_path"), 
              action = "store",
	      type = "character",
              default=NA, 
              help="Job_Output path to analyse"),

  make_option(c("-o", "--out_path"),
  	      action = "store",
              type = "character",
              default=NA, 
              help="Folder where report is registered"),

  make_option(c("-n", "--name"),
  	      action = "store",
              type = "character",
              default="Optimize_ressource_report", 
              help="Name of the report file [default %default]"),
  
  # make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
  #             help="Make the program not be verbose."),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Should the program print extra stuff out? [default %default]")  
)

opt = parse_args(OptionParser(option_list=option_list))






#output path not specified
if(is.na(opt$out_path)){
  opt$out_path <- opt$in_path
}

#main part of script
if(file.exists(opt$in_path) & file.exists(opt$out_path)) {
  
  # Options values
  job_output_path = opt$in_path
  output_path = opt$out_path
  
  
  
  ############## RESULT ###########################################################
  
  #function call with complete list of .o file
  result <- parsed_folder(job_output_path)
  
  Info_df <- result[[1]]
  DF_plot <- result[[2]]
  
  # The following dataframe contain average values parsed in .o files contained in the job_output file gived.
  #print(Info_df)
  
  # The following dataframe contain all values found and will be used to create next plots
  #print(DF_plot)
  
  ############## Register file ###################################################
  ##Register dataframe as CSV
  #Create name with date
  actual_date_time <- strsplit(x= as.character(Sys.time()), split = " ")[[1]][1:2] %>%
                      paste(collapse ="T") %>% 
                      gsub(pattern = ":", replacement = ".")
  file_name <- paste(c(opt$name, actual_date_time), collapse ="_")
  complete_path_name <- paste(c(output_path, file_name), collapse ="/")
  complete_path_name_csv <- paste(c(complete_path_name,"csv"), collapse =".")

  #Write csv in job_output folder
  write.csv(Info_df, complete_path_name_csv, row.names=FALSE)
  ############## Custom specfific dataFrame for ploting ##########################
  #Change WaitingTime, RunTime and TimeLimit format
  
  Waiting_unite <- Keep_hour(DF_plot$WaitingTime)[[2]]
  DF_plot$WaitingTime <- Keep_hour(DF_plot$WaitingTime)[[1]]
  
  Run_Limit <- Keep_hour(DF_plot$RunTime, DF_plot$TimeLimit)
  
  DF_plot$RunTime <- Run_Limit[[1]]
  DF_plot$TimeLimit <- Run_Limit[[2]]     #RunTime and TimeLimit have the same format depending on RunTime values
  Time_unite <- Run_Limit[[3]]            #unite associated with RunTime and TimeLimit columns for plot
  
  #computing percentage amount for RunTime compared to TimeLimit
  DF_plot$RunTime_Efficiency <- round((DF_plot$RunTime / DF_plot$TimeLimit) *100, 1)
  
  #compute maximal Efficiency for each step
  DF_max_Eff <- DF_plot[c(1, 8)] %>% group_by(JobName) %>% top_n(1, RunTime_Efficiency)
  
  #Merge DF_max_Eff with DF_plot
  DF_plot <- merge(DF_plot, DF_max_Eff, by = "JobName")
  
  #Change Memory_Request column type to numeric
  DF_plot$Memory_Request <- as.numeric(DF_plot$Memory_Request)
  
  ############## WaitingTime Plot ################################################

  p_WaintingTime <- ggplot(DF_plot, aes(x=as.factor(JobName), y=WaitingTime, label= as.numeric(WaitingTime))) + 
    theme(panel.background = element_rect(fill = 'white', color = 'grey'), 
          panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) +
    
    #boxplot
    geom_boxplot(color="blue", 
                 fill="#69b3a2",
                 alpha=0.2) +     
    
    scale_x_discrete(expand = c(0.05, 0)) +
    
    #invert x and y axes
    coord_flip() + 
    
    #label on highest values
    geom_text(data = . %>% group_by(JobName) %>% filter(WaitingTime == max(WaitingTime)), 
              nudge_y = 0.45,
              check_overlap = TRUE) +
    
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
    
    scale_x_discrete(expand = c(0.05, 0)) +
    
    coord_flip() +
    
    xlab("Step name") +
    
    ylab(paste(c("RunTime (", Time_unite, ")"), collapse ="")) +
    
    geom_boxplot( aes(y=RunTime),
                  alpha=0.2,
                  color="blue",
                  fill="#69b3a2",
    ) + 
    
    #Change RunTime_Efficiency.y values and axes to keep them between 0 and 100%
    geom_point( aes(y= round(RunTime_Efficiency.y *100 / max(RunTime),1)),
                color="red",
                alpha=0.5) +
    
    ggtitle("RunTime and RunTime_Efficiency") +
    
    #scale_x_discrete(expand = c(0.05, 0)) + 
    
    geom_text(y =  as.numeric(round(DF_plot$RunTime_Efficiency.y *100 / max(DF_plot$RunTime),1)),
              label = as.numeric(round(DF_plot$RunTime_Efficiency.y *100 / max(DF_plot$RunTime),1)),
              color="red",
              size=2.5,
              nudge_x = 0.5, 
              nudge_y = 0,
              check_overlap = TRUE) +
    
    scale_y_continuous(
      
      # Features of the first axis
      name = paste(c("RunTime (", Time_unite, ")"), collapse =""),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~. *10 / max(DF_plot$RunTime),
                          name="RunTime_Efficiency (percentage)")
    )
  
  
  
  ############## Memory_Efficiency vs Memory_Request #############################
  
  # print(DF_plot)
  # print(typeof(DF_plot$Memory_Efficiency))
  # print(typeof(DF_plot$Memory_Request))




  p_Memory <- ggplot(DF_plot, aes(x=as.factor(JobName))) +                               #CHANGER DF_MEMORY PAR DF_plot

    theme(panel.background = element_rect(fill = 'white', color = 'grey'),
          panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),
          axis.title.x.top = element_text(color = "red", size=13),
          axis.title.x.bottom = element_text(color = "blue", size=13),
          axis.title.y = element_text(size=13),
    ) +

    scale_x_discrete(expand = c(0.05, 0)) +

    coord_flip() +

    xlab("Step name") +

    geom_boxplot( aes(y=DF_plot$Memory_Request),
                  alpha=0.1,
                  color="blue",
                  fill="#69b3a2",
    ) +

    geom_point( aes(y=round(Memory_Efficiency * 100 / max(DF_plot$Memory_Request),1)),
                color="red",
                alpha=0.5) +

    ggtitle("Memory_Request and RunTime_Efficiency") +



    geom_text( y = round(as.numeric(DF_plot$Memory_Efficiency * 100 / max(DF_plot$Memory_Request)),1),
               label = round(as.numeric(DF_plot$Memory_Efficiency * 100 / max(DF_plot$Memory_Request)),1),
               color="red",
               size=3,
               nudge_x = 0.5,
               check_overlap = TRUE) +

    scale_y_continuous(

      # Features of the first axis
      name = "Memory_Request (GB)" ,

      # Add a second axis and specify its features
      sec.axis = sec_axis(~. ,
                          name="Memory_Efficiency (%)")
    )
  
  ############## PDF Result #######################################################
  #create path and file name
  pdf_file_name <- paste(c(opt$name,"plots", actual_date_time), collapse ="_")
  pdf_complete_path_name <- paste(c(output_path, pdf_file_name), collapse ="/")            
  pdf_complete_path_name <- paste(c(pdf_complete_path_name,"pdf"), collapse =".")
  
  #create pdf file from plots
  pdf(pdf_complete_path_name)
  print(p_WaintingTime)     # Page 1
  print(p_RunTime)          # Page 2 
  print(p_Memory)           # Page 3
  dev.off()
  
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
  
  #verbose option
  if (opt$verbose){
    cat("\n")  
    #Total steps
    cat(paste(nrow(DF_plot),"files parsed." ,"\n", sep = " "))
    
    #Total differents steps
    cat(paste(nrow(Info_df),"differents steps in total among all files." ,"\n", sep = " "))
    
    #Absolute path if the one gived is relative
    #cat(paste("", opt$in_path,"\n", sep = " "))
    cat("\n")
  }
  
} 

#Input folder problem
if  (!file.exists(opt$in_path)){
  cat(paste("The specified input path doesn't exist :", opt$in_path,"\n", sep = " "), file=stderr()) # print error messages to stderr
}

#Output folder problem
if  (!file.exists(opt$out_path)){
  cat(paste("The specified output path doesn't exist :", opt$out_path,"\n", sep = " "), file=stderr()) # print error messages to stderr
}






 









