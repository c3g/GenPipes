#!/usr/bin/env python

### Edouard Henrion (2017/07/18)
### job2json

import os
import sys
import getopt
import re
import json
import subprocess
import datetime

def getarg(argument):
    step_name = ""
    job_name = ""
    job_log = ""
    job_done = ""
    json_files = ""
    status = True

    optli,arg = getopt.getopt(argument[1:], "s:j:l:d:o:f:h", ['step_name', 'job_name', 'job_log', 'job_done', 'json_outfiles', 'status', 'help'])

    if len(optli) == 0 :
        usage()
        sys.exit("Error : No argument given")

    for option, value in optli:
        if option in ("-s", "--step_name"):
            if str(value) == "" :
                sys.exit("Error - step_name (-s, --step_name) not provided...\n" + str(value))
            else :
                step_name = str(value)
        if option in ("-j", "--job_name"):
            if str(value) == "" :
                sys.exit("Error - job_name (-n, --job_name) not provided...\n")
            else :
                job_name = str(value)
        if option in ("-l", "--job_log"):
            if str(value) == "" :
                sys.exit("Error - job_log (-l, --job_log) not provided...\n")
            else :
                job_log = str(value)
        if option in ("-d", "--job_done"):
            if str(value) == "" :
                sys.exit("Error - job_done (-b, --job_done) not provided...\n")
            else :
                job_done = str(value)
        if option in ("-o", "--json_outfiles"):
            if str(value) == "" :
                sys.exit("Error - json_outfiles (-j, --json_outfiles) not provided...\n")
            else :
                json_files = str(value)
        if option in ("-f", "--status"):
            status = str(value)
        if option in ("-h", "--help"):
            usage()
            sys.exit()

    return step_name, job_name, job_log, job_done, json_files, status

def usage():
    print "\n-------------------------------------------------------------------------------------------"
    print "job2json.py will append a JSON section describing a pipeline job that has just finished"
    print "to a JSON file which was pre-generated when the pipeline was launched."
    print "This script is usually launched automatically before and after each pipeline job."
    print "This program was written by Edouard HENRION"
    print "For more information, contact: edouard.henrion@computationalgenomics.ca"
    print "-------------------------------------------------------------------------------------------\n"
    print "USAGE : job2json.py [option] "
    print "       -s    --step_name     : name of the step of the current job"
    print "       -j    --job_name      : name of the current job"
    print "       -l    --job_log       : name of the log file for the current job"
    print "       -d    --job_done      : name of the done file for the current job"
    print "       -o    --json_outfiles : comma-separated list of names of json files which need to be appended affected by the current job"
    print "       -f    --status        : boolean value to indicate if the job has failed (False/0) or succeeded (True/1) - Default : True"
    print "       -h    --help          : this help \n"

def main():
    #print "command line used :\n" + " ".join(sys.argv[:])

    step_name, job_name, job_log, job_done, json_files, status = getarg(sys.argv)

    for jfile in json_files.split(","):
        with open(jfile, 'r') as json_file:
            current_json = json.load(json_file)

        # Make sure the job_log file is not in absolute path anymore
        if current_json['pipeline']['general_information']['analysis_folder']:
            job_log = re.sub(current_json['pipeline']['general_information']['analysis_folder'], "", job_log)

        step_found = False
        job_found = False

        # Find the current step which should already exist in the json object (already created by bfx/jsonator.py)
        for jstep in current_json['pipeline']['step']:
            if jstep['name'] == step_name:
                step_found = True

                # Find the current job which should already exist in the json object (already created by bfx/jsonator.py)
                for jjob in jstep['job']:
                    if jjob['name'] == job_name:
                        job_found = True
                        if  status == "running":
                            jjob['job_start_date'] = re.sub("\.\d+$", "", str(datetime.datetime.now()))
                            jjob['status'] = "running"
                        else:
                            jjob['log_file'] = job_log
                            if status == "0":
                                jjob['status'] = "success"
                                jjob['done_file'] = job_done
                            else:
                                jjob['status'] = "error"
                            jjob['job_end_date'] = re.sub("\.\d+$", "", str(datetime.datetime.now()))

                # If job does not exists already, raise an exception
                if not job_found :
                    sys.exit("Error : job " + job_name + ", within step " + step_name + ", was not found in " + " json_file...")

        # If step does not exists already, raise an exception
        if not step_found :
            sys.exit("Error : step " + step_name + " was not found in " + " json_file...")

        # Print to file
        with open(jfile, 'w') as out_json:
            json.dump(current_json, out_json, indent=4)

main()
