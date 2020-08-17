'''
Combines all the JSONs within a certain directory.

Traverses the folders inside the given folder, and extracts the 
JSON data.

Takes 2 command line arguments.

python json_combine.py <PATH_TO_SOFT_STACK> <PATH_TO_STORE_JSON>
'''
import os
import sys
import json


soft_stack = sys.argv[1]
if len(sys.argv) == 3:
    json_folder = sys.argv[2]
else:
    json_folder = os.getcwd()

# soft_stack = '/lustre03/project/6007512/yatharth/mugqic/'
# json_folder = '/lustre03/project/6007512/yatharth/C3G_2/json_folder'

directories_ = os.listdir(soft_stack)
list_of_software = []

for elem in directories_:
    if os.path.isdir(os.path.join(soft_stack, elem)):
        list_of_software.append(elem)

software_paths = [os.path.join(soft_stack, soft) for soft in list_of_software]

def getVersions(soft, dir_):
    files = os.listdir(dir_)
    files = [file for file in files if file[0] is not '.']
    versions = []
    for version in files:
        version 
        versions.append(version)
    return versions

softList = {}
for index in range(0, len(software_paths)):
    extractedVersions = getVersions(list_of_software[index], software_paths[index])
    currSoft = list_of_software[index]
    softList[currSoft] = extractedVersions

def addToJson(dict_):
    for json_p in software_paths:
        load_json = None
        json_f = os.path.join(json_p, '.metadata.json')
        if os.path.isfile(json_f):
            with open(json_f, 'r') as f:
                load_json = json.load(f)
            sw_name = os.path.split(json_p)[1]
            sw_name = sw_name.lower()
            sw_name = sw_name.replace('.json', '')
            sw_name = sw_name.replace('NOTFOUND_', '')
            sw_name = sw_name.replace('VERIFY_', '')
            for elem in dict_:
                if elem.lower() == sw_name:
                    print("{} found".format(sw_name))
                    load_json['VERSIONS_AVAILABLE'] = dict_[elem]
                    with open(json_f, 'w') as f:
                        json.dump(load_json, f, indent=6)

addToJson(softList)
                
arr_ = []
for path in software_paths:
    json_p = os.path.join(path, '.metadata.json')
    if os.path.isfile(json_p):
        with open(json_p, 'r') as f:
            arr_.append(json.load(f))

with open(os.path.join(json_folder, 'combined_json.json'), 'w') as f:
    json.dump(arr_, f, indent=6)
