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

from spreadsheet_utils import parse_google_sheet, print_sample_sheet

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--project', help="Project ID", required=True, dest="project_id")
    parser.add_argument('-r', '--run', help="Run ID", required=True, dest="run_id")
    parser.add_argument('-l', '--lane', help="Lane number", type=int, required=True, dest="lane_number")
    parser.add_argument('-o', '--outfile', help="Output file (csv format). Default : PROJECT_ID.RUN_ID.LANE_NUMBER.sample_sheet.csv", required=False)
    parser.add_argument('-a', '--authentication_file', help="JSON authentication file used to connect the spreadsheets", required=True, type=file, dest="json_file")
    parser.add_argument('--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    project = args.project_id
    run = args.run_id
    lane = args.lane_number

    outfile = args.outfile
    authentication_file = args.json_file.name

    MGI_spreadsheet_name = "ALL MGI Run Management"
    dict_of_columns = parse_google_sheet(
        MGI_spreadsheet_name,
        authentication_file
    )

    print_sample_sheet(
        dict_of_columns,
        project,
        run,
        lane,
        outfile=outfile
    )

