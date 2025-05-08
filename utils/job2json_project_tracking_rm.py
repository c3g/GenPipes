#!/usr/bin/env python3

### Paul Stretenowich (2025/04/29)
### job2json

import argparse
import os
import errno
import sys
import re
import json
import time
import random
import shutil
import signal

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# MUGQIC Modules
# from core.config import *

def main():

    parser = argparse.ArgumentParser(prog='job2json_project_tracking_rm.py', description="Removes a JSON section to a JSON file which was pre-generated when the pipeline was launched. This script is usually launched in the context of cleaning GenPipes output files. /!\\ This version is for project tracking database only.")
    parser.add_argument('-f', '--file_regex', required=True, help="regex to match the location_uri to be removed")
    parser.add_argument('-o', '--json_outfile', required=True, help="name of json output file")
    args = parser.parse_args()

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

        files_to_remove = []
        for sample in current_json['sample']:
            for readset in sample['readset']:
                for job in readset['job']:
                    if 'file' in job:
                        for file in job['file']:
                            if re.search(args.file_regex, file['location_uri'], re.UNICODE):
                                files_to_remove.append((job, file))
        for job, file in files_to_remove:
            job['file'].remove(file)
            if len(job['file']) == 0:
                job.pop('file')
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
