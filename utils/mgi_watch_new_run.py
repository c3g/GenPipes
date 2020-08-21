#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import argparse
import logging
import gspread
import os
import textwrap

from spreadsheet_utils import print_sample_sheet, parse_google_sheet

logger = logging.getLogger(__name__)

columns = {}

def get_flowcell_from_run(
    columns,
    run
    ):
    for x in range(len(columns['RUN_ID'])):
        if columns['RUN_ID'][x] == run:
            if columns['Flowcell_ID'][x]:
                return columns['Flowcell_ID'][x]
    # If no flowcell was found
    return None

def get_project_from_run(
    columns,
    run
    ):
    for x in range(len(columns['RUN_ID'])):
        if columns['RUN_ID'][x] == run:
            if columns['Project_ID'][x]:
                return columns['Project_ID'][x]

def get_lanes_from_run(
    columns,
    run
    ):
    lanes = []
    for x in range(len(columns['RUN_ID'])):
        if columns['RUN_ID'][x] == run:
            if columns['Lane'][x] and columns['Lane'][x] not in lanes:
                lanes.append(columns['Lane'][x])
    return lanes

def compare_runs(columns):
    header = "RUN_ID"
    if os.path.isfile("/tmp/mgi_runs.ref"):
        new_list = filter(None, sorted(set(columns['RUN_ID'])))
#        ref_list = []
        ref_list = filter(None, open("/tmp/mgi_runs.ref", "r").read().split("\n"))

        if new_list != ref_list:
            print "New runs have been added to the Run Management spreadsheet... building the new samples sheets now"
            runs_added = list(set(new_list) - set(ref_list))
            print runs_added
            for run in runs_added:
                flowcell = get_flowcell_from_run(
                    columns,
                    run
                )
                if flowcell:
                    print "Processing run " + run
                    project = get_project_from_run(
                        columns,
                        run
                    )
                    print project
                    lanes = get_lanes_from_run(
                        columns,
                        run
                    )

                    for lane in lanes:
                        print "    Lane " + lane
                        print_sample_sheet(
                            columns,
                            project,
                            run,
                            lane
                        )
            # replace referece run list by the current run list to set it as the reference for next watch round
            print_runs(columns)
        else:
            print "No new run detected..."

    else:
        # print current run list and set it as the reference for next watch round
        print_runs(columns)
        ref_list = open("/tmp/mgi_runs.ref", "r").read().split("\n")
        print ref_list

def print_runs(columns, path="/tmp/mgi_runs.ref"):
    f = open(path, "wb+")
    for i in filter(None, sorted(set(columns['RUN_ID']))):
        f.write(i + "\n")
    f.close

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--spreadsheet_name', help='Name of the MGI run management Google spreadsheet to parse', required=True, dest="spreadsheet_name")
    parser.add_argument('-a', '--authentication_file', help="JSON authentication file used to connect the spreadsheets", type=file, required=True, dest="json_file")
    parser.add_argument('-o', '--outfile', help="Output file (csv format)", required=False)
    parser.add_argument('--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    MGI_spreadsheet_name = args.spreadsheet_name
    outfile = args.outfile
    authentication_file = args.json_file.name

    dict_of_columns = parse_google_sheet(
        MGI_spreadsheet_name,
        authentication_file
    )

    compare_runs(
        dict_of_columns
    )
