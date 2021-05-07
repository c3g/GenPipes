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

logger = logging.getLogger(__name__)

def read_json(filepath):
    return json.loads(read_file(filepath))

def read_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    return content

def bold(text):
    return '\x1b[1m%s\x1b[21m' % text

def red(text):
    return '\x1b[31m%s\x1b[39m' % text

def yellow(text):
    return '\x1b[33m%s\x1b[39m' % text

def get_versions(
    soft,
    dir_
    ):

    files = os.listdir(dir_)
    files = [file for file in files if file[0] is not '.']
    versions = []
    for version in files:
        version 
        versions.append(version)
    return versions

def add_to_jsons(
    dict_
    ):

    for json_p in software_paths:
        load_json = None
        json_f = os.path.join(json_p, '.metadata.json')
        tmp_json_f = os.path.join(tmp_stack, os.path.basename(json_p), '.metadata.json')
        if os.path.isfile(json_f):
            load_json = read_json(json_f)
            sw_name = os.path.split(json_p)[1]
            sw_name = sw_name.lower()
            sw_name = sw_name.replace('.json', '')
            sw_name = sw_name.replace('NOTFOUND_', '')
            sw_name = sw_name.replace('VERIFY_', '')
            for elem in dict_:
                if elem.lower() == sw_name:
                    logger.info("%s found" % sw_name)
                    load_json['versions'] = dict_[elem]
                    # Create the current soft folder in the temporary stack
                    try:
                        os.makedirs(os.path.dirname(tmp_json_f), exist_ok = True)
                        logger.info('Directory %s created successfully' % os.path.dirname(tmp_json_f))
                    except OSError as error:
                        print(error)
                        logger.error(red('Directory %s cannot be created. Skipping...' % os.path.dirname(tmp_json_f)))
                    with open(tmp_json_f, 'w') as f:
                        json.dump(load_json, f, indent=6)
        else:
            logger.warning(yellow('Missing metadata : %s' % json_f))

def send_file(
    filepath,
    url
    ):

    try:
        data = read_json(filepath)
    except Exception as error:
        print(error)
        logger.error(red('Failed to read file : %s ' % filepath))
        return

    logger.info("Sending file %s to %s..." % (filepath, url))

    try:
        response = requests.post(url, json=data)
        result = response.json()
    except Exception as error:
        print(error)
        logger.error(red('Got error while sending file : %s ' % filepath))
        return

    if response.status_code == 200 and result.get('ok') is True:
        logger.info('File %s sent successfully to %s' % (filepath, url))
    else:
        logger.error(red('Request failed %d ' % response.status_code) + ('[%s] %s: %s : %s' % (bold(url), filepath, response.reason, response.text)))
        return

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Script to combine all .metadata jsons in the C3G software stack.')
    parser.add_argument('-p', '--path', type=str, help='Path to the software stack', required=True)
    parser.add_argument('-o', '--output_file', type=str, help='Path of the output JSON file', required=True)
    parser.add_argument('-u', '--url', type=str, help='URL where send the JSON files to', required=True)
    parser.add_argument('-l', '--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    soft_stack = args.path
    json_output = args.output_file
    url = args.url

    # Setting of a temporary/writable stack folder
    tmp_stack = os.path.join(os.path.dirname(json_output), '.tmp')

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
    for path in tmp_software_paths:
        json_p = os.path.join(path, '.metadata.json')
        if os.path.isfile(json_p):
            with open(json_p, 'r') as f:
                arr_.append(json.load(f))

    with open(json_output, 'w') as f:
        json.dump(arr_, f, indent=6)

    send_file(
        json_output,
        url
    )
