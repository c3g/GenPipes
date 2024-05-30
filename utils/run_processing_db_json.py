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
    read_dict = {}
    for json_file in json_list:
        read_dict[json_file] = json.load(open(json_file,
                            mode="r"))
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
        if self.is_same_value_between_jsons("run_obj_id"):
            return self._json_dict[sorted(self._json_dict)[0]]["run_obj_id"]

    @property
    def run_name(self):
        if self.is_same_value_between_jsons("run"):
            return "".join([
                self._json_dict[sorted(self._json_dict)[0]]["run"],
                "-",
                self._json_dict[sorted(self._json_dict)[0]]["seqtype"],
                ])

    @property
    def run_instrument(self):
        if self.is_same_value_between_jsons("seqtype"):
            return self._json_dict[sorted(self._json_dict)[0]]["seqtype"]

    def properties(self):
        class_items = self.__class__.__dict__.items()
        return dict((k, getattr(self, k)) 
                    for k, v in class_items 
                    if isinstance(v, property))

    def is_same_value_between_jsons(self, run_validation_field):
        values = []
        for run_validation_json in self._json_dict.keys():
            if self._json_dict[run_validation_json][run_validation_field] not in values:
                values.append(self._json_dict[run_validation_json][run_validation_field])
        if len(values) == 1:
            return True
        else:
            LOGGER.error("Multiple {run_validation_field} detected".format(
                run_validation_field=run_validation_field,
                ))
            raise Exception("Conflict in input values for \
\"{run_validation_field}\". Confirm that the input files are \
from a single run".format(
                run_validation_field=run_validation_field,
                ))

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
    databases = {}
    for project in projects.items():
        databases[project] = database_json(args, project)
        databases[project].operation_platform
        print("TEST CURRENT:",
              #[a for a in dir(databases[project]) if not a.startswith('_')],
              databases[project].properties(),
              databases[project].run_name
              )
    print("End of main()")


if __name__ == "__main__":
    main()
