#!/usr/bin/env python
from __future__ import print_function
import os
import os.path
import sys
import json
import time
import getopt
import shutil
import requests
from deepdiff import DeepDiff

def main():

    options = get_arguments()

    if not os.path.isdir(options.cache_folder):
        try:
            os.makedirs(options.cache_folder)
        except Exception as e:
            print(e)
            print(red('Couldnt create cache folder at %s' % options.cache_folder))
            sys.exit(1)

    print('Checking for JSON files in %s' % options.watch_folder)

    if options.update_interval is None:
        run(options)
        sys.exit()

    while True:
        time.sleep(options.update_interval)
        run(options)


def run(options):
    files = [f for f in os.listdir(options.watch_folder) if f.endswith('.json')]

    print('Found %i files' % len(files))

    details = []
    for i, filename in enumerate(files):
        filepath = os.path.join(options.watch_folder, filename)
        sample_name = None
        filename_parts = filename.split('.')

        if len(filename_parts) < 4:
            # Old filename format: $USER.[UUID].json: we need to read the content to know the sample_name
            if not os.path.isfile(filepath):
                print(red('  File %s does not exist anymore. Skipping' % filename))
                continue

            try:
                content = read_file(filepath)
            except Exception as e:
                print(e)
                print(red('  Failed to read file "%s". Skipping' % filename))
                continue

            try:
                data = json.loads(content)
            except Exception as e:
                print(e)
                print(red('  Failed to parse JSON "%s". Skipping' % filename))
                continue

            sample_name = data['sample_name']
        else:
            # New filename format: $USER.[sample_name].[UUID].json
            sample_name = '.'.join(filename_parts[1:-2])

        print('  Read %s (%i/%i)' % (filename, i + 1, len(files)))
        details.append({
            "filepath": filepath,
            "sample_name": sample_name
        })

    print('Read %i files' % len(details))

    details_by_sample = group_by(details, lambda detail: detail['sample_name'])
    for sample_name in details_by_sample:
        details = details_by_sample[sample_name]
        details.sort(key=lambda detail: os.path.getmtime(detail['filepath']))
        details.reverse()

    print('Found %i samples' % len(details_by_sample.keys()))

    for sample_name in details_by_sample:
        send_files(options, sample_name, details_by_sample[sample_name])

def send_files(options, sample_name, details):

    cache_filepath = os.path.join(options.cache_folder, '%s.json' % sample_name)

    # First detail is of the newest file
    detail = details[0]
    filepath = detail['filepath']

    try:
        data = read_json(filepath)
    except Exception as e:
        print(e)
        print(red('Failed to read file "%s": ' % filepath))
        return

    previous_data = None
    if os.path.isfile(cache_filepath):
        previous_data = read_json(cache_filepath)

    if not os.path.isfile(filepath):
        print('File %s not on disk anymore. Skipping sample.' % filepath)
        return

    url = None
    username = filepath.split('/')[-1].split('.')[0]
    should_send = True

    if not previous_data:
        url = options.url + '/api/samples/external-update/' + username
    else:
        url = options.url + '/api/samples/external-update-diff/' + username
        operations = get_diff(previous_data, data)

        if len(operations) == 0:
            print(yellow('No difference for file %s. No request made.' % filepath))
            should_send = False
        else:
            data = {
                'sample_name': data['sample_name'],
                'operations': operations,
            }

    if should_send:
        try:
            response = requests.post(url, json=data)
            result   = response.json()
        except Exception as e:
            print(red('Got error while sending file. Skipping.'))
            print(e)
            print(options)
            print('Url: ' + url)
            print('Filepath: ' + detail['filepath'])
            #  if response:
            #     print('Response:\n----------\n' + response.text + '\n----------')
            return

        if response.status_code == 200 and result.get('ok') is True:
            print('Sent %s (deleting %i files)' % (filepath, len(details)))
        else:
            print(red('Request failed %d ' % response.status_code) +
                ('[%s] %s: %s : %s' % (bold(url), filepath, response.reason, response.text)))
            return

    shutil.copy(filepath, cache_filepath)
    for detail in details:
        os.remove(detail['filepath'])

def read_json(filepath):
    return json.loads(read_file(filepath))

def read_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    return content

def write_file(filename, content):
    with open(filename, 'w') as file:
        file.write(content)

def usage():
    print("watch_portal_folder.py - watches a folder for JSON files and sends them to and external URL.")
    print("")
    print(bold("Usage:") + " watch_portal_folder.py [options] ")
    print("")
    print(bold("Options:"))
    print("    -w, --watch     - folder to watch")
    print("    -c, --cache     - folder for caching JSONs")
    print("    -u, --url       - URL to send the JSON files to")
    print("    -i, --interval  - folder check interval (in seconds) (default: None - doesnt watch)")
    print("    -h, --help      - display this message")

def get_arguments():
    options = dotdict({})
    options.watch_folder    = './buffer'
    options.cache_folder    = '/tmp/watch_portal_folder'
    options.url             = 'http://localhost:3000'
    options.update_interval = None

    optli, arg = getopt.getopt(sys.argv[1:], 'w:c:u:i:h', ['watch=', 'cache=', 'url=', 'interval=', 'help'])

    if len(optli) == 0:
        usage()
        exit('Error: No arguments given')

    for option, value in optli:
        if option in ('-w', '--watch'):
            if str(value) == '':
                exit('Error: --watch folder not provided\n')
            else:
                options.watch_folder = str(value)
        if option in ('-c', '--cache'):
            if str(value) == '':
                exit('Error: --cache folder not provided\n')
            else:
                options.cache_folder = str(value)
        if option in ('-u', '--url'):
            if str(value) == '':
                exit('Error: --url not provided\n')
            else:
                options.url = str(value)
        if option in ('-i', '--interval'):
            if int(value) == '':
                exit('Error: --interval not provided\n')
            else:
                options.update_interval = int(value)
        if option in ('-h', '--help'):
            usage()
            exit()

    return options

def bold(text):
    return '\x1b[1m%s\x1b[21m' % text

def red(text):
    return '\x1b[31m%s\x1b[39m' % text

def yellow(text):
    return '\x1b[33m%s\x1b[39m' % text

def exit(message=None):
    sys.exit(red(message) if message else None)

class dotdict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def group_by(seq, key=lambda x: x):
    result = {}
    for value in seq:
        result.setdefault(key(value), []).append(value)
    return result

def get_diff(old, new):
    operations = []

    result = DeepDiff(old, new)

    #  { 'type_changes': { 'root[2]': { 'new_type': <class 'str'>,
    #                                  'new_value': '2',
    #                                  'old_type': <class 'int'>,
    #                                  'old_value': 2}}}
    type_changes = result.get('type_changes', {})
    for key in type_changes:
        operations.append({'op': 'replace', 'path': get_path_from_key(key), 'value': type_changes[key]['new_value']})

    #  {'values_changed': {'root[2]': {'new_value': 4, 'old_value': 2}}}
    values_changed = result.get('values_changed', {})
    for key in values_changed:
        operations.append({'op': 'replace', 'path': get_path_from_key(key), 'value': values_changed[key]['new_value']})

    #  {'dictionary_item_added': {'root[5]', 'root[6]'}}
    dictionary_item_added = result.get('dictionary_item_added', {})
    for key in dictionary_item_added:
        operations.append({'op': 'add', 'path': get_path_from_key(key), 'value': get_value_from_key(new, key)})

    #  {'dictionary_item_removed': {'root[4]'}}
    dictionary_item_removed = result.get('dictionary_item_removed', {})
    for key in dictionary_item_removed:
        operations.append({'op': 'delete', 'path': get_path_from_key(key)})

    #  {'iterable_item_removed': {"root[4]['b'][2]": 3, "root[4]['b'][3]": 4}}
    iterable_item_removed = result.get('iterable_item_removed', {})
    for key in iterable_item_removed:
        path = get_path_from_key(key)
        index = path[-1]
        path = path[:-1]
        operations.append({'op': 'remove', 'path': path, 'index': index})

    #  {'iterable_item_added': {"root[4]['b'][3]": 3}}
    iterable_item_added = result.get('iterable_item_added', {})
    for key in iterable_item_added:
        path = get_path_from_key(key)
        index = path[-1]
        path = path[:-1]
        value = iterable_item_added[key]
        operations.append({'op': 'insert', 'path': path, 'index': index, 'value': value})

    return operations

def get_path_from_key(key):
    return map(lambda x: eval(x), key[5:-1].split(']['))

def get_value_from_key(root, key):
    return eval(key)


# Benchmarking helpers

def get_current_milliseconds():
    return int(round(time.time() * 1000))

labels = {}

def start_timer(label):
    labels[label] = get_current_milliseconds()

def end_timer(label):
    print('%s: %i' % (label, get_current_milliseconds() - labels[label]))


if __name__ == '__main__':
    main()
