#!/usr/bin/env python

### Paul Stretenowich (2023/06/05)
### job2json_project_tracking

import argparse
import os
import errno
import sys
import json
import time
import random
import shutil
import signal

from datetime import datetime

def main():

    parser = argparse.ArgumentParser(prog='job2json_project_tracking.py', description="Appends a JSON section describing a pipeline job that has just finished to a JSON file which was pre-generated when the pipeline was launched. This script is usually launched automatically before and after each pipeline job. /!\\ This version is for project trasking database only.")
    parser.add_argument('-s', '--sample_names', required=True, help="comma-separated list of names of the samples of the current job")
    parser.add_argument('-r', '--readset_names', required=True, help="comma-separated list of names of the readsets of the current job")
    parser.add_argument('-j','--job_name', required=True, help="name of the current job")
    parser.add_argument('-m', '--metrics', required=False, help="comma-separated list of metrics of the current job: name=value,name=value,... With <name> = metric name; <value> = metric value")
    # parser.add_argument('-l','--job_log', required=True, help="name of the log file for the current job")
    # parser.add_argument('-d','--job_done', required=True, help="name of the done file for the current job")
    # parser.add_argument('-c','--config', required=True, help="")
    # parser.add_argument('-u','--user', required=True, help="name of user running GenPipes")
    parser.add_argument('-o','--json_outfile', required=True, help="name of json output file")
    parser.add_argument('-f', '--status', required=False, help="status of job")
    args = parser.parse_args()

    # step_name, job_name, job_log, job_done, json_files, config_files, user, status = getarg(sys.argv)
    # job_log = args.job_log
    sample_list = args.sample_names.split(",")
    readset_list = args.readset_names.split(",")
    if args.metrics:
        metrics_list = args.metrics.split(",")
    else:
        metrics_list = []

    # finally (unlock) will execute even if exceptions occur
    try:
        # Make sure the args.json_outfile is unlock if process receive SIGTERM too (not python exception)
        def sigterm_handler(_signo, _stack_frame):
            unlock(args.json_outfile)
            sys.exit(0)
        signal.signal(signal.SIGTERM, sigterm_handler)

        # First lock the file to avoid multiple and synchronous writing attemps
        lock(args.json_outfile)

        with open(args.json_outfile, 'r') as json_file:
            current_json = json.load(json_file)

        for sample in current_json['sample']:
            if sample['sample_name'] in sample_list:
                for readset in sample['readset']:
                    if readset['readset_name'] in readset_list:
                        for job in readset['job']:
                            if job['job_name'] == args.job_name:
                                if args.status == "RUNNING":
                                    job['job_start'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                                    job['job_status'] = args.status
                                else:
                                    if args.status == "0":
                                        job['job_status'] = "COMPLETED"
                                    else:
                                        job['job_status'] = "FAILED"
                                    job['job_stop'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                                for metric in metrics_list:
                                    metric_split = metric.split("=")
                                    metric_name = metric_split[0]
                                    metric_value = metric_split[1]
                                    try:
                                        job['metric'].append({
                                            'metric_name': metric_name,
                                            'metric_value': metric_value
                                            })
                                    except KeyError:
                                        job['metric'] = [{
                                            'metric_name': metric_name,
                                            'metric_value': metric_value
                                            }]

        # Print to file
        with open(args.json_outfile, 'w') as out_json:
            json.dump(current_json, out_json, indent=4)

    finally:
        # Finally unlock the file
        unlock(args.json_outfile)

def lock(filepath):
    """
    Locking filepath by creating a .lock file
    """
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
