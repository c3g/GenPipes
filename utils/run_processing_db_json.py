#!/usr/bin/env python

import os
import re
import sys
import json
import logging
import argparse

def arg_parse():
    """
    Argument parser.

    Returns argparse parsed object.
    """
    parser = argparse.ArgumentParser(
        prog=os.path.basename(__file__),
        description="Aggregate database information in json format following \
                     successful run processing",
        )
    parser.add_argument(
        "json",
        nargs="+",
        help="Run validation json file(s)"
        )
    parser.add_argument(
        "-p", "--platform",
        help="Operation platform (Default: Abacus)",
        default="abacus",
        )
    parser.add_argument(
        "-v", "--verbose",
        help="Increase verbose level (Default: Warning)",
        action='count',
        default=0,
        )
    parser.add_argument(
        "-q", "--quiet",
        help="Decrease verbose level (Levels: Info, Warning)",
        action='count',
        default=0,
        )
    return parser.parse_args()

def read_run_jsons(json_list):
    """
    Convert the list of jsons files supplied via arguments into a dictionary
    of json parsed objects.

    Args:
        json_list   list
            Paths of json files to load in the dictionary. These should be
            run_validation json obtained from GenPipes' run_processing.

    Returns:
        Dictionary of the json, with content available to call with keys. 
    """
    print("As list:", json_list)
    read_dict = {}
    for json_file in json_list:
        read_dict[json_file] = json.load(open(json_file,
                            mode="r"))
    print("As keys:", read_dict.keys())
    return read_dict

class database_json(object):

    def __init__(self, args, project):
        print("Init", self)
        self._json_dict =  read_run_jsons(args.json)
        self._args = args
        self._project = project

    @property
    def input_files(self):
        return self._json_dict.keys()

    @property
    def operation_platform(self):
        platform = self._args.platform
        if not re.search(self._args.platform, os.uname()[1]):
            LOGGER.warning("""Operation platform supplied is \
\"{platform}\", but running from a different hostname \"{host}\". The \
specified platform will be used but it is recommended to run on the same \
platform.""".format(
                    platform=platform,
                    host=os.uname()[1])
                )
        return platform

    @property
    def project_fms_id(self):
        return self._project[0]

    @property
    def project_name(self):
        return self._project[1]

    @property
    def run_fms_id(self):
        run_obj_id = []
        for run_validation_json in self._json_dict.keys():
            if self._json_dict[run_validation_json]["run_obj_id"] not in run_obj_id:
                run_obj_id.append(self._json_dict[run_validation_json]["run_obj_id"])
        if len(run_obj_id) == 1:
            return run_obj_id
        else:
            LOGGER.error("Multiple runs detected")
            raise Exception("Single run at the time")

    @property
    def run_name(self): #TODO incomplete value And need to turn the checks in a function
        run_name = []
        for run_validation_json in self._json_dict.keys():
            if self._json_dict[run_validation_json]["run"] not in run_name:
                run_name.append(self._json_dict[run_validation_json]["run"])
        if len(run_name) == 1:
            return run_name
        else:
            LOGGER.error("Multiple run names detected")
            raise Exception("Single run at the time")

    @property
    def run_instrument(self): #TODO Need to turn the checks in a function
        run_instrument = []
        for run_validation_json in self._json_dict.keys():
            if self._json_dict[run_validation_json]["seqtype"] not in run_instrument:
                run_instrument.append(self._json_dict[run_validation_json]["seqtype"])
        if len(run_instrument) == 1:
            return run_instrument
        else:
            LOGGER.error("Multiple run instrument detected")
            raise Exception("Single run at the time")

def main():
    args = arg_parse()
    global LOGGER 
    logging.basicConfig(level=30-args.verbose*10+args.quiet*10)
    LOGGER = logging.getLogger(__name__)
    json_dict = read_run_jsons(args.json)
    # Need to loop through the values that must be unique by database entry:
    #   project based info: project_id, project_name
    #   run based info:     run_name, run_date, etc
    projects = {}
    for json_path, json_content in json_dict.items():
        for sample, content in json_content["readsets"].items():
            if content["project_obj_id"] not in projects.keys() \
                    or content["project_name"] not in projects.values():
                projects[content["project_obj_id"]] = content["project_name"]
                LOGGER.info(
                    "Project ID:{obj_id} \"{project}\" identified.".format(
                        obj_id=content["project_obj_id"],
                        project=content["project_name"],
                        ),
                    )
    print(projects)
    databases = {}
    for project in projects.items():
        databases[project] = database_json(args, project)
        databases[project].operation_platform
        print("TEST CURRENT:",
              [a for a in dir(databases[project]) if not a.startswith('_')],
              databases[project].run_name,
              )
    print("End of main()")


if __name__ == "__main__":
    main()
