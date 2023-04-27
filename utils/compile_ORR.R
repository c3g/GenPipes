#!/usr/bin/env Rscript
#
# This script is made to 
#
#

suppressMessages(library(optparse))
suppressMessages(library(dplyr))


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
              help="Folder where report is registered. If any folder is specified, the output folder will be the same as input folder."),
  
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
  if(opt$verbose){
    cat("\n")
    cat("Output folder will be the same as output folder.")
  }
  
  opt$out_path <- opt$in_path
}

#main part of script
if(file.exists(opt$in_path) & file.exists(opt$out_path)) {
  
  # Options values
  job_output_path = opt$in_path
  output_path = opt$out_path
  
  #creation of datetime element
  actual_date_time <- strsplit(x= as.character(Sys.time()), split = " ")[[1]][1:2] %>%
    paste(collapse ="T") %>% 
    gsub(pattern = ":", replacement = ".")
  
  #creation of file name
  file_name <- paste(c(opt$name,"plots", actual_date_time), collapse ="_")
  file_name_html <- paste(c(file_name,"html"), collapse =".")
  
  #call to Rmd file
  rmarkdown::render(input = "optimize_resource_report.Rmd",
                    params = list(
                      input = opt$in_path,             #joboutput path
                      output = opt$out_path,           #report output path
                      name = opt$name,                 #name of document (if needed)
                      verbose = opt$verbose            #verbose option
                      ),
                    output_dir = opt$out_path,
                    output_file = file_name_html
                    )
  
} 

#Input folder problem
if  (!file.exists(opt$in_path)){
  cat(paste("The specified input path doesn't exist :", opt$in_path,"\n", sep = " "), file=stderr()) # print error messages to stderr
}

#Output folder problem
if  (!file.exists(opt$out_path)){
  cat(paste("The specified output path doesn't exist :", opt$out_path,"\n", sep = " "), file=stderr()) # print error messages to stderr
}















