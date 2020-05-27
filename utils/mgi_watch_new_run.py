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

logger = logging.getLogger(__name__)

columns = {}

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

def compare_runs(columns):
    header = "RUN_ID"
    if os.path.isfile("/tmp/mgi_runs.ref"):
        new_list = filter(None, sorted(set(columns['RUN_ID'])))
        ref_list = open("/tmp/mgi_runs.ref", "r").read().split("\n")
        print new_list
        print ref_list
        if new_list != ref_list:
            # replace referece run list by the current run list to set it as the reference for next watch round
            print "Lists differ !!!"
            print_runs()

        #os.remove()
    else:
        # print current run list and set it as the reference for next watch round
        print_runs()
        ref_list = open("/tmp/mgi_runs.ref", "r").read().split("\n")
        print ref_list

def print_runs(path="/tmp/mgi_runs.ref"):
    f = open(path, "wb+")
    for i in sorted(set(columns['RUN_ID'])):
        f.write(i)
    f.close

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', help="Output file (csv format)", required=False)
    parser.add_argument('--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    outfile = args.outfile

    dict_of_columns = parse_google_sheet()

    compare_runs(
        dict_of_columns
    )
