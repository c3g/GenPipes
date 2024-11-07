#!/usr/bin/env python

### Edouard Henrion (2017/07/18)
### Paul Stretenowich (2024/09/27)
### job2json

import argparse
import os
import errno
import sys
import re
import json
import datetime
import time
import random
import shutil
import signal

from uuid import uuid4

# GenPips Modules
from genpipes.core.config import global_conf

def main(args=None):
    """
    Main function
    """
    if args is None:
        parser = argparse.ArgumentParser(prog='job2json.py', description="Appends a JSON section describing a pipeline job that has just finished to a JSON file which was pre-generated when the pipeline was launched. This script is usually launched automatically before and after each pipeline job. /!\\ This version is for project tracking database only.")
        parser.add_argument('-s', '--step_name', required=True, help="name of the step of the current job")
        parser.add_argument('-j','--job_name', required=True, help="name of the current job")
        parser.add_argument('-l','--job_log', required=True, help="name of the log file for the current job")
        parser.add_argument('-d','--job_done', required=True, help="name of the done file for the current job")
        parser.add_argument('-o','--json_outfile', required=True, help="comma-separated list of names of json files which need to be appended by the current job")
        parser.add_argument('-f','--status', required=True, help="job status comming from GenPipes. Either running when job starts or $GenPipes_STATE")
        parser.add_argument('-u','--user', required=True, help="user to use GenPipes $USER")
        args = parser.parse_args()

    global_conf.parse_files(args.config_files)

    for jfile in args.json_files.split(","):
        # finally (unlock) will execute even if exceptions occur
        try:

            # Make sure the jfile is unlock if process receive SIGTERM too (not python exception)
            def sigterm_handler(_signo, _stack_frame):
                unlock(jfile)
                sys.exit(0)
            signal.signal(signal.SIGTERM, sigterm_handler)

            # First lock the file to avoid multiple and synchronous writing attemps
            lock(jfile)

            with open(jfile, 'r') as json_file:
                current_json = json.load(json_file)
            json_file.close()

            if not current_json:
                os.remove(jfile)

            else:
                # Make sure the job_log file is not in absolute path anymore
                if current_json['pipeline']['general_information']['analysis_folder']:
                    job_log = re.sub(current_json['pipeline']['general_information']['analysis_folder'], "", args.job_log)
                else:
                    job_log = args.job_log

                step_found = False
                job_found = False

                # Find the current step which should already exist in the json object (already created by bfx/jsonator.py)
                for jstep in current_json['pipeline']['step']:
                    if jstep['name'] == args.step_name:
                        step_found = True

                        # Find the current job which should already exist in the json object (already created by bfx/jsonator.py)
                        for jjob in jstep['job']:
                            if jjob['name'] == args.job_name:
                                job_found = True
                                if  args.status == "running":
                                    jjob['job_start_date'] = re.sub(r"\.\d+$", "", str(datetime.datetime.now()))
                                    jjob['status'] = "running"
                                else:
                                    jjob['log_file'] = job_log
                                    if args.status == "0":
                                        jjob['status'] = "success"
                                        jjob['done_file'] = args.job_done
                                    else:
                                        jjob['status'] = "error"
                                    jjob['job_end_date'] = re.sub(r"\.\d+$", "", str(datetime.datetime.now()))

                        # If job does not exist already, raise an exception
                        if not job_found :
                            sys.exit(f"Error : job {args.job_name}, within step {args.step_name}, was not found in json_file...")

                # If step does not exist already, raise an exception
                if not step_found :
                    sys.exit(f"Error : step {args.step_name} was not found in json_file...")

                # Let's do 5 attempts to write the file (because sometimes we weirly end up with malformed JSON files...)
                count = 5
                while count:
                    # Print to file
                    with open(jfile, 'w') as out_json:
                        json.dump(current_json, out_json, indent=4)
                    out_json.close()

                    # Test opening the written file
                    try:
                        with open(jfile, 'r') as json_file:
                            current_json_hash = json.load(json_file)
                        if current_json_hash:
                            # Print a copy of the JSON file for the monitoring interface
                            portal_output_dir = global_conf.global_get('DEFAULT', 'portal_output_dir', required=False, param_type='dirpath')
                            if portal_output_dir != '':
                                shutil.copy(jfile, os.path.join(portal_output_dir, f"{args.user}.{current_json['sample_name']}.{str(uuid4().hex)}.json"))
                            count = 0
                        else:
                            count -= 1
                    except json.decoder.JSONDecodeError:
                        count -= 1

            # Print a copy of it for the monitoring interface
            portal_output_dir = global_conf.global_get('DEFAULT', 'portal_output_dir', required=False, param_type='dirpath')
            if portal_output_dir != '':
                with open(os.path.join(portal_output_dir, f"{args.user}.{current_json['sample_name']}.{str(uuid4().hex)}.json"), 'w') as out_json:
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
            else:
                # An unexpected error has occured : let's stop the program and raise the error"
                raise exception
        else:
            # The lock folder was successfully created !"
            unlocked = False

def unlock(filepath):
    """
    Unlocking filepath by removing the .lock file
    """
    shutil.rmtree(filepath + '.lock', ignore_errors=True)

if __name__ == '__main__':
    main()
