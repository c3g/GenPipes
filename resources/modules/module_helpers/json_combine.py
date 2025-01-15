#!/usr/bin/env python
'''
Combines all the JSONs within a certain directory.

Traverses the folders inside the given folder, and extracts the 
JSON data.

Takes 2 command line arguments.

python json_combine.py <PATH_TO_SOFT_STACK> <PATH_TO_STORE_JSON>
'''
import os
import json
import logging
import argparse
import requests

log = logging.getLogger(__name__)

def read_json(filepath):
    return json.loads(read_file(filepath))

def read_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    return content

def red(text):
    return '\x1b[31m%s\x1b[39m' % text

def yellow(text):
    return '\x1b[33m%s\x1b[39m' % text

def get_versions(
    soft,
    dir_
    ):

    files = os.listdir(dir_)
    files = [file for file in files if file[0] != '.']
    versions = []
    for version in files:
        version
        versions.append(version)
    versions.sort(key=lambda s: [int(u) for u in s.split('.') if u.isdigit()], reverse=True)
    if len(versions) > 3:
      return versions[:3] + ['...']
    else:
      return versions[:3]

def add_to_jsons(
    dict_
    ):

    for json_p in software_paths:
        load_json = None
        json_f = os.path.join(json_p, '.metadata.json')
        tmp_json_f = os.path.join(tmp_stack, os.path.basename(json_p), '.metadata.json')
        if os.path.isfile(json_f):
            log.info("%s found" % json_f)
            load_json = read_json(json_f)
            sw_name = os.path.split(json_p)[1]
            sw_name = sw_name.lower()
            sw_name = sw_name.replace('.json', '')
            sw_name = sw_name.replace('NOTFOUND_', '')
            sw_name = sw_name.replace('VERIFY_', '')
            for elem in dict_:
                if elem.lower() == sw_name:
                    log.info("%s found" % sw_name)
                    load_json['versions'] = dict_[elem]
                    # Create the current soft folder in the temporary stack
                    try:
                        os.makedirs(os.path.dirname(tmp_json_f), exist_ok = True)
                        log.info('Directory %s created successfully' % os.path.dirname(tmp_json_f))
                    except OSError as error:
                        print(error)
                        log.error(red('Directory %s cannot be created. Skipping...' % os.path.dirname(tmp_json_f)))
                    with open(tmp_json_f, 'w') as f:
                        json.dump(load_json, f, indent=6)
        else:
            log.warning(yellow('Missing metadata : %s' % json_f))

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Script to combine all .metadata jsons in the C3G software stack.')
    parser.add_argument('-p', '--path', type=str, help='Path to the software stack', required=True)
    parser.add_argument('-o', '--output_file', type=str, help='Path of the output JSON file. Default: None i.e. prints in stdout')
    parser.add_argument('-l', '--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    soft_stack = args.path
    json_output = args.output_file

    # Setting of a temporary/writable stack folder
    if json_output:
        tmp_stack = os.path.join(os.path.dirname(json_output), '.tmp')
    else:
        tmp_stack = os.path.join('/tmp', '.tmp')

    list_of_software = [soft for soft in os.listdir(soft_stack) if os.path.isdir(os.path.join(soft_stack, soft))]
    software_paths = [os.path.join(soft_stack, soft) for soft in list_of_software]
    tmp_software_paths = [os.path.join(tmp_stack, soft) for soft in list_of_software]

    soft_list = {}
    for index in range(0, len(software_paths)):
        extracted_versions = get_versions(list_of_software[index], software_paths[index])
        curr_soft = list_of_software[index]
        soft_list[curr_soft] = extracted_versions

    add_to_jsons(soft_list)

    arr_ = []
    for path in sorted(tmp_software_paths, key=str.casefold):
        json_p = os.path.join(path, '.metadata.json')
        if os.path.isfile(json_p):
            with open(json_p, 'r') as f:
                arr_.append(json.load(f))

    if json_output:
        with open(json_output, 'w') as f:
            json.dump(arr_, f, indent=6)
    else:
        print(json.dumps(arr_))
