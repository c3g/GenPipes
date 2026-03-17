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
import hashlib
import logging

from datetime import datetime

log = logging.getLogger(__name__)

def main(args=None):
    """
    Main function
    """
    if args is None:
        parser = argparse.ArgumentParser(prog='job2json_project_tracking.py', description="Appends a JSON section describing a pipeline job that has just finished to a JSON file which was pre-generated when the pipeline was launched. This script is usually launched automatically before and after each pipeline job. /!\\ This version is for project tracking database only.")
        parser.add_argument('-s', '--sample_names', required=True, help="comma-separated list of names of the samples of the current job")
        parser.add_argument('-r', '--readset_names', required=True, help="comma-separated list of names of the readsets of the current job")
        parser.add_argument('-j','--job_name', required=True, help="name of the current job")
        parser.add_argument('-m', '--metrics', required=False, help="comma-separated list of metrics of the current job: name=value,name=value,... With <name> = metric name; <value> = metric value")
        parser.add_argument('-o','--json_outfile', required=True, help="name of json output file")
        parser.add_argument('-f', '--status', required=False, help="status of job")
        args = parser.parse_args()

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
                                for fobj in job.get('file', []):
                                    uri = fobj.get("location_uri")
                                    local_path = uri_to_local_path(uri)

                                    if local_path and os.path.exists(local_path) and os.path.isfile(local_path):
                                        md5 = compute_md5sum(local_path)
                                    else:
                                        md5 = None

                                    fobj["file_md5sum"] = md5
                        
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

def compute_md5sum(filepath, block_size=8 * 1024 * 1024):
    """
    Compute MD5 hex digest of filepath
    """
    try:
        with open(filepath, 'rb') as fh:
            h = hashlib.md5()
            for chunk in iter(lambda: fh.read(block_size), b''):
                h.update(chunk)
            return h.hexdigest()
    except Exception as exc:
        log.debug ("compute_md5sum: could not read '%s': %s", filepath, exc)
        return None
    
def uri_to_local_path(uri):
    """
    Convert location_uri to local path (removing name of cluster)
    """
    if uri is None:
        return None
    parts = uri.split("://", 1)
    if len(parts) == 2:
        return parts[1]   
    return uri

if __name__ == '__main__':
    main()
