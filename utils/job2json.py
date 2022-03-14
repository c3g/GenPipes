#!/usr/bin/env python

### Edouard Henrion (2017/07/18)
### job2json

import os
import errno
import sys
import getopt
import re
import json
import datetime
import time
import random
import shutil
import signal
import uuid

from uuid import uuid4

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# MUGQIC Modules
from core.config import *

def getarg(argument):
    step_name = ""
    job_name = ""
    job_log = ""
    job_done = ""
    json_files = ""
    config_files = []
    user = ""
    status = True

    options, _ = getopt.getopt(argument[1:], "s:j:l:d:o:c:u:f:h", ['step_name', 'job_name', 'job_log', 'job_done', 'json_outfiles', 'config', 'user', 'status', 'help'])

    if len(options) == 0:
        usage()
        sys.exit("Error : No argument given")

    for option, value in options:
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
        if option in ("-c", "--config"):
            if str(value) == "" :
                sys.exit("Error - config_files (-c, --config) not provided...\n")
            else :
                config_files = str(value).split(',')
        if option in ("-u", "--user"):
            if str(value) == "" :
                sys.exit("Error - user (-u, --user) not provided...\n")
            else :
                user = str(value)
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

    return step_name, job_name, job_log, job_done, json_files, config_files, user, status

def usage():
    print("\n-------------------------------------------------------------------------------------------")
    print("job2json.py will append a JSON section describing a pipeline job that has just finished")
    print("to a JSON file which was pre-generated when the pipeline was launched.")
    print("This script is usually launched automatically before and after each pipeline job.")
    print("This program was written by Edouard HENRION")
    print("For more information, contact: edouard.henrion@computationalgenomics.ca")
    print("-------------------------------------------------------------------------------------------\n")
    print("USAGE : job2json.py [option] ")
    print("       -s    --step_name     : name of the step of the current job")
    print("       -j    --job_name      : name of the current job")
    print("       -l    --job_log       : name of the log file for the current job")
    print("       -d    --job_done      : name of the done file for the current job")
    print("       -o    --json_outfiles : comma-separated list of names of json files which need to be appended by the current job")
    print("       -f    --status        : boolean value to indicate if the job has failed (False/0) or succeeded (True/1) - Default : True")
    print("       -h    --help          : this help \n")

def main():

    step_name, job_name, job_log, job_done, json_files, config_files, user, status = getarg(sys.argv)

    #print config_files
    global_config_parser.parse_files(config_files)


    for jfile in json_files.split(","):

        # finally (unlock) will execute even if exceptions occur
        try:

            # Make sure the the jfile is unlock if process receive SIGTERM too (not python exception)
            def sigterm_handler(_signo, _stack_frame):
                unlock(jfile)
                sys.exit(0)
            signal.signal(signal.SIGTERM, sigterm_handler)

            # First lock the file to avoid multiple and synchronous writing attemps
            lock(jfile)

            with open(jfile, 'r') as json_file:
                current_json = json.load(json_file)
            json_file.close()

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

                    # If job does not exist already, raise an exception
                    if not job_found :
                        sys.exit("Error : job " + job_name + ", within step " + step_name + ", was not found in " + " json_file...")

            # If step does not exist already, raise an exception
            if not step_found :
                sys.exit("Error : step " + step_name + " was not found in " + " json_file...")

            # Print to file
            with open(jfile, 'w') as out_json:
                json.dump(current_json, out_json, indent=4)

            out_json.close()

            # Print a copy of it for the monitoring interface
            portal_output_dir = global_config_parser.param('DEFAULT', 'portal_output_dir', required=False, param_type='dirpath')
            if portal_output_dir != '':
                with open(os.path.join(portal_output_dir, user + '.' + current_json['sample_name'] + '.' + str(uuid4().hex) + '.json'), 'w') as out_json:
                    json.dump(current_json, out_json, indent=4)
        finally:
            # Finally unlock the file
            unlock(jfile)

def lock(filepath):
    unlocked = True
    while unlocked:
        try:
            os.makedirs(filepath + '.lock')
        except OSError as exception:
            if exception.errno == errno.EEXIST and os.path.isdir(filepath + '.lock'):
                # The lock folder already exists, we need to wait for it to be deleted
                sleep_time = random.randint(1, 100)
                time.sleep(sleep_time)
                pass
            else:
                # An unexpected error has occured : let's stop the program and raise the error"
                raise exception
        else:
            # The lock folder was successfully created !"
            unlocked = False

def unlock(filepath):
    shutil.rmtree(filepath + '.lock', ignore_errors=True)

if __name__ == '__main__':
    main()
