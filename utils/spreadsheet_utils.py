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

# Python Standard Modules
import argparse
import logging
import gspread
import csv
import os
import textwrap

logger = logging.getLogger(__name__)

def parse_google_sheet(sheet_name, authentication_file):
    gc = gspread.service_account(authentication_file)

    wks = gc.open(sheet_name).sheet1

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
            columns[h].append(v.encode('utf8').lstrip().rstrip().replace(" ", "_"))
    # Add the Readset column
    columns["Readset"] = [sample + "_" + library for sample, library in zip(columns["Sample_Name"], columns["Library"])]

    # Split the Suencer column in to 2 columns
    columns["SequencerID"] = [seq.split('-')[1] if "-" in seq else "" for seq in columns["Sequencer"]]
    columns["Sequencer"] = [seq.split('-')[0]  if "-" in seq else seq for seq in columns["Sequencer"]]

    return columns

def print_sample_sheet(
    columns,
    project,
    run,
    lane,
    output_dir=None,
    group=True,
    verbose=None
    ):

    if project not in columns['Project_ID']:
        logger.error("Project ID" + str(project) + " was not found in the spreadsheet...")
    if run not in columns['RUN_ID']:
        logger.error("Run ID " + str(run) + " was not found in the spreadsheet...")
    if lane and str(lane) not in columns['Lane']:
        logger.error("Lane " + str(lane) + " was not found in the spreadsheet...")

    process_dir = os.path.join(str(project) , str(run))
    if group:
        outfile = str(project) + "." + str(run) + ".sample_sheet.csv"
    elif lane:
        outfile = str(project) + "." + str(run) + ".L0" + str(lane) + ".sample_sheet.csv"
    else:
        print "Error :  no lane... no group..."
        exit(10)

    if output_dir :
        process_dir = os.path.join(output_dir, process_dir)
    if not os.path.exists(process_dir):
        os.makedirs(process_dir)
    outfile = os.path.join(process_dir, outfile)

    f = csv.writer(open(outfile, "wb+"))

    # Write CSV Header
    header = ["Sample_Name", "Readset", "Library", "Project", "Project_ID", "Protocol", "Library_Source", "Index", "Pool_ID", "RUN_ID", "Flowcell_ID", "Lane", "Run_Date", "Sequencer", "SequencerID"]
    f.writerow(header)

    for x in range(len(columns['RUN_ID'])):
        if columns['Project_ID'][x] == project and columns['RUN_ID'][x] == run and (group or (lane and columns['Lane'][x] == str(lane))):
            if not columns['Library_Source'][x]:
                logger.error("Missing Library Source in run '" + run + "' (sample '" + columns['Sample'][x] + "') Skipping...")
                break
            if not columns['Flowcell_ID'][x]:
                logger.error("Missing Flowcell ID in run '" + run + "' (sample '" + columns['Sample'][x] + "') Skipping...")
                break
            if columns['Pool_ID'][x] != "FAIL":
                f.writerow([columns[col][x] for col in header])
            else:
                logger.info("Skipping failed sample " + columns['Sample'][x])

