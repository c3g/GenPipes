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
import csv
import textwrap

logger = logging.getLogger(__name__)

def parse_google_sheet():
    gc = gspread.service_account()

    wks = gc.open("COVID 19 MGI Run Management ").sheet1

    # As the labels in the header of the sample sheet differs from the ones in the header of the google spread sheet
    # we work with indices
    list_of_lists = wks.get_all_values()

    # Use the second line of the google sheet as the HEADER,
    # strip trailing space characters and replace new-line charaters by spaces
    header = [x.strip() for x in list_of_lists[1]]
    header = [x.replace("\n", " ") for x in header]

    # Now build a dict of lists to store the columns
    columns = {}
    for i in range(len(header)):
        h = header[i].encode('utf8')
        if h == "Library Plate Barcode-Well (Library ID)":
            h = "Library"
        if h == "Index Name":
            h = "Index"
        h = h.replace(" ", "_")
        columns[h] = []
        header[i] = h

    # Actual data starts at line 3 
    for row in list_of_lists[2:]:
        for h, v in zip(header, row):
            columns[h].append(v.encode('utf8'))
    # Add the Readset column
    columns["Readset"] = [sample + "_" + library for sample, library in zip(columns["Sample_Name"], columns["Library"])] 

    # Split the Suencer column in to 2 columns
    columns["SequencerID"] = [seq.split('-')[1] if "-" in seq else "" for seq in columns["Sequencer"]]
    columns["Sequencer"] = [seq.split('-')[0]  if "-" in seq else seq for seq in columns["Sequencer"]]

    return columns

def print_sheet(
    columns,
    project,
    run,
    lane,
    outfile=None,
    verbose=None
    ):

    if project not in columns['Project']:
        logger.error("Project " + str(project) + " was not found in " + columns['Project'])
    if run not in columns['RUN_ID']:
        logger.error("Run ID " + str(run) + " was not found in " + columns['Run_ID'])
    if str(lane) not in columns['Lane']:
        logger.error("Lane " + str(lane) + " was not found in " + columns['Lane'])

    if not outfile:
        outfile = str(project) + "." + str(run) + ".L0" + str(lane) + ".sample_sheet.csv"

    f = csv.writer(open(outfile, "wb+"))

    # Write CSV Header
    header = ["Sample_Name", "Readset", "Library", "Project", "Project_ID", "Protocol", "Index", "Pool_ID", "RUN_ID", "Flowcell_ID", "Lane", "Run_Date", "Sequencer", "SequencerID"]
    f.writerow(header)

    for x in range(len(columns['RUN_ID'])):
        if columns['Project'][x] == project and columns['RUN_ID'][x] == run and columns['Lane'][x] == str(lane):
            if columns['Pool_ID'][x] != "FAIL":
                f.writerow([columns[col][x] for col in header])
            else:
                logger.info("Skipping failed sample " + columns['Sample'][x])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--project', help="Project ID", required=True, dest="project_id")
    parser.add_argument('-r', '--run', help="Run ID", required=True, dest="run_id")
    parser.add_argument('-l', '--lane', help="Lane number", type=int, required=True, dest="lane_number")
    parser.add_argument('-o', '--outfile', help="Output file (csv format). Default : PROJECT_ID.RUN_ID.LANE_NUMBER.sample_sheet.csv", required=False)
    parser.add_argument('--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    project = args.project_id
    run = args.run_id
    lane = args.lane_number

    outfile = args.outfile

    dict_of_columns = parse_google_sheet()

    print_sheet(
        dict_of_columns,
        project,
        run,
        lane,
        outfile=outfile
    )

