#!/usr/bin/env python

### Edouard Henrion (2017/07/18)
### job2json

import os
import sys
import argparse
import getopt
import re
import json
import subprocess
import datetime
import time
from uuid import uuid4

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# MUGQIC Modules
from core.config import *

# def getarg(argument):
#     step_name = ""
#     job_name = ""
#     job_log = ""
#     job_done = ""
#     json_files = ""
#     config_files = []
#     user = ""
#     status = True
#
#     options, _ = getopt.getopt(argument[1:], "s:j:l:d:o:c:u:f:h", ['step_name', 'job_name', 'job_log', 'job_done', 'json_outfiles', 'config', 'user', 'status', 'help'])
#
#     if len(options) == 0:
#         usage()
#         sys.exit("Error : No argument given")
#
#     for option, value in options:
#         if option in ("-s", "--step_name"):
#             if str(value) == "" :
#                 sys.exit("Error - step_name (-s, --step_name) not provided...\n" + str(value))
#             else :
#                 step_name = str(value)
#         if option in ("-j", "--job_name"):
#             if str(value) == "" :
#                 sys.exit("Error - job_name (-n, --job_name) not provided...\n")
#             else :
#                 job_name = str(value)
#         if option in ("-l", "--job_log"):
#             if str(value) == "" :
#                 sys.exit("Error - job_log (-l, --job_log) not provided...\n")
#             else :
#                 job_log = str(value)
#         if option in ("-d", "--job_done"):
#             if str(value) == "" :
#                 sys.exit("Error - job_done (-b, --job_done) not provided...\n")
#             else :
#                 job_done = str(value)
#         if option in ("-c", "--config"):
#             if str(value) == "" :
#                 sys.exit("Error - config_files (-c, --config) not provided...\n")
#             else :
#                 config_files = str(value).split(',')
#         if option in ("-u", "--user"):
#             if str(value) == "" :
#                 sys.exit("Error - user (-u, --user) not provided...\n")
#             else :
#                 user = str(value)
#         if option in ("-o", "--json_outfiles"):
#             if str(value) == "" :
#                 sys.exit("Error - json_outfiles (-j, --json_outfiles) not provided...\n")
#             else :
#                 json_files = str(value)
#         if option in ("-f", "--status"):
#             status = str(value)
#         if option in ("-h", "--help"):
#             usage()
#             sys.exit()
#
#     return step_name, job_name, job_log, job_done, json_files, config_files, user, status
#
# def usage():
#     print "\n-------------------------------------------------------------------------------------------"
#     print "job2json.py will append a JSON section describing a pipeline job that has just finished"
#     print "to a JSON file which was pre-generated when the pipeline was launched."
#     print "This script is usually launched automatically before and after each pipeline job."
#     print "This program was written by Edouard HENRION"
#     print "For more information, contact: edouard.henrion@computationalgenomics.ca"
#     print "-------------------------------------------------------------------------------------------\n"
#     print "USAGE : job2json.py [option] "
#     print "       -s    --step_name     : name of the step of the current job"
#     print "       -j    --job_name      : name of the current job"
#     print "       -l    --job_log       : name of the log file for the current job"
#     print "       -d    --job_done      : name of the done file for the current job"
#     print "       -o    --json_outfiles : comma-separated list of names of json files which need to be appended affected by the current job"
#     print "       -f    --status        : boolean value to indicate if the job has failed (False/0) or succeeded (True/1) - Default : True"
#     print "       -h    --help          : this help \n"

def main():
    #print "command line used :\n" + " ".join(sys.argv[:])

    parser = argparse.ArgumentParser()

    # Options for job2json.py

    parser.add_argument("-u", "--user", help="name of the user processing pipeline")
    parser.add_argument("-s", "--steps", help="name of the step of the current job")
    parser.add_argument("-j", "--job_name", help="name of current job")
    parser.add_argument("-l", "--job_log", help="name of the log file for the current job")
    parser.add_argument("-d", "--job_done", help="name of the done file for the current job")
    parser.add_argument("-o", "--json_outfiles", help="output directory (default: current)", nargs="+", type=file)
    #parser.add_argument("-o", "--json_outfiles", help="output directory (default: current)")
    parser.add_argument("-c", "--config",
                        help="config INI-style list of files; config parameters are overwritten based on files order",
                        nargs="+", type=file)
    parser.add_argument("-f", "--status",
                        help="boolean value to indicate if the job has failed (False/0) or succeeded (True/1) - Default : True")
                        #action="store_true")

    args = parser.parse_args()
    args.logLevel = "INFO"

    #step_name, job_name, job_log, job_done, json_files, config_files, user, status = getarg(sys.argv)

    global_config_parser.parse_files(args.config)

    job_log = args.job_log
    
    for jfile in args.json_outfiles:
        #wait_for_lock(jfile)
        #lock(jfile)
        
        #with open(jfile, 'r') as json_file:
        current_json = json.load(jfile)

        # Make sure the job_log file is not in absolute path anymore
        if current_json['pipeline']['general_information']['analysis_folder']:
            job_log = re.sub(current_json['pipeline']['general_information']['analysis_folder'], "", job_log)

        step_found = False
        job_found = False

        # Find the current step which should already exist in the json object (already created by bfx/jsonator.py)
        for jstep in current_json['pipeline']['step']:
            if jstep['name'] == args.steps:
                step_found = True

                # Find the current job which should already exist in the json object (already created by bfx/jsonator.py)
                for jjob in jstep['job']:
                    if jjob['name'] == args.job_name:
                        job_found = True
                        if  args.status == "running":
                            jjob['job_start_date'] = re.sub("\.\d+$", "", str(datetime.datetime.now()))
                            jjob['status'] = "running"
                        else:
                            jjob['log_file'] = job_log
                            if args.status == "0":
                                jjob['status'] = "success"
                                jjob['done_file'] = args.job_done
                            else:
                                jjob['status'] = "error"
                            jjob['job_end_date'] = re.sub("\.\d+$", "", str(datetime.datetime.now()))

                # If job does not exists already, raise an exception
                if not job_found :
                    sys.exit("Error : job " + args.job_name + ", within step " + args.steps + ", was not found in " + " json_file...")

        # If step does not exists already, raise an exception
        if not step_found :
            sys.exit("Error : step " + args.steps + " was not found in " + " json_file...")

        # Print to file
        with open(jfile, 'w') as out_json:
            json.dump(current_json, out_json, indent=4)

        # Print a copy of it for the monitoring interface
        portal_output_dir = global_config_parser.param('DEFAULT', 'portal_output_dir', required=False, param_type='dirpath')
        if portal_output_dir != '':
            with open(os.path.join(portal_output_dir, args.user + '.' + uuid4().get_hex() + '.json'), 'w') as out_json:
                json.dump(current_json, out_json, indent=4)
                
def wait_for_lock(filepath):
    while os.path.isfile(filepath + '.lock'):
        time.sleep(1)

def lock(filepath):
    with open(filepath + '.lock', 'w') as file:
        file.write('')

def unlock(filepath):
    try:
        os.remove(filepath + '.lock')
    except:
        pass


main()
