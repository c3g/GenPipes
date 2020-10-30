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
import subprocess
import csv

from spreadsheet_utils import print_sample_sheet, parse_google_sheet

logger = logging.getLogger(__name__)

columns = {}

def get_sequencer_from_run(
    columns,
    run
    ):
    for x in range(len(columns['RUN_ID'])):
        if columns['RUN_ID'][x] == run:
            if columns['Sequencer'][x]:
                seq = columns['Sequencer'][x]
                if seq == "01":
                    return "seq1/R2130400190016"
                elif seq == "02":
                    return "seq2/R2130400190018"
    # If no sequencer was found
    return None

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

def compare_runs(
    columns,
    outdir,
    process_dir,
    genpipes_repo,
    mgi_runs_file
    ):

    header = "RUN_ID"
    if os.path.isfile(mgi_runs_file):
        new_list = filter(None, sorted(set(columns['RUN_ID'])))
#        ref_list = []
        ref_list = filter(None, open(mgi_runs_file, "r").read().split("\n"))

        if new_list != ref_list:
            print "New runs have been added to the Run Management spreadsheet... building the new samples sheets now"
            runs_added = list(set(new_list) - set(ref_list))
            print runs_added
            for run in runs_added:
                sequencer_path = get_sequencer_from_run(
                    columns,
                    run
                )
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
                        print "        Generating sample sheet..."
                        print_sample_sheet(
                            columns,
                            project,
                            run,
                            lane,
                            outdir
                        )
                        print "        Generating GenPipes script..."
                        print_genpipes_scripts(
                            outdir,
                            process_dir,
                            project,
                            flowcell,
                            run,
                            lane,
                            sequencer_path,
                            genpipes_repo
                        )
            # replace referece run list by the current run list to set it as the reference for next watch round
            print_runs(columns)

        else:
            print "No new run detected..."

    else:
        # print current run list and set it as the reference for next watch round
        print_runs(
            columns,
            mgi_runs_file
        )
        ref_list = open(mgi_runs_file, "r").read().split("\n")
        print ref_list

def print_runs(
    columns,
    path
    ):

    f = open(path, "wb+")
    for i in filter(None, sorted(set(columns['RUN_ID']))):
        f.write(i + "\n")
    f.close

def print_genpipes_scripts(
    samplesheet_dir,
    process_dir,
    project,
    flowcell,
    run,
    lane,
    sequencer_path,
    genpipes_repo
    ):

    # Check if --raw-fastq' flag should be used int the pipeline call
    sample_sheet = os.path.join(samplesheet_dir, project, run, "L0"+lane, project+"."+run+".L0"+lane+".sample_sheet.csv")
    sample_sheet_rows = [row for row in csv.DictReader(open(sample_sheet, 'rb'), delimiter=',')]
    print len(sample_sheet_rows)
    print sample_sheet_rows[0]['Index']
    if len(sample_sheet_rows) == 1 and "MGI" in sample_sheet_rows[0]['Index']:
        raw_fastq = True
    else: 
        raw_fastq = False

    if not os.path.exists(process_dir+"/../genpipes_scripts/"+project+"/"+run):
        os.makedirs(process_dir+"/../genpipes_scripts/"+project+"/"+run)
    genpipes_script = open(process_dir+"/../genpipes_scripts/"+project+"/"+run+"/"+project+"."+run+".L0"+lane+".genpipes_script.sh", 'wb+') 
    genpipes_script.write("""\
mkdir -p {process_dir}/{project}/{run}/L0{lane} && \\
python {genpipes_repo}/pipelines/mgi_run_processing/mgi_run_processing.py \\
  -c {genpipes_repo}/pipelines/mgi_run_processing/mgi_run_processing.base.ini $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \\
  --no-json -l debug {raw_fastq}\\
  -d /nb/Research/MGISeq/{sequencer_path}/{flowcell} \\
  -r {samplesheet_dir}/{project}/{run}/L0{lane}/{project}.{run}.L0{lane}.sample_sheet.csv \\
  --lane {lane} \\
  -o {process_dir}/{project}/{run}/L0{lane} \\
  > {process_dir}/{project}/{run}/{project}.{run}.L0{lane}.sh \\
  2> {process_dir}/{project}/{run}/{project}.{run}.L0{lane}.trace.log""".format(
        process_dir=process_dir,
        samplesheet_dir=samplesheet_dir,
        genpipes_repo=genpipes_repo,
        raw_fastq="--raw-fastq " if raw_fastq else "",
        project=project,
        flowcell=flowcell,
        run=run,
        lane=lane,
        sequencer_path=sequencer_path
    ))
    genpipes_script.close()
    subprocess.call("bash "+process_dir+"/../genpipes_scripts/"+project+"/"+run+"/"+project+"."+run+".L0"+lane+".genpipes_script.sh", shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--spreadsheet_name', help='Name of the MGI run management Google spreadsheet to parse', required=False, dest="spreadsheet_name", default='ALL MGI Run Management')
    parser.add_argument('-a', '--authentication_file', help="JSON authentication file used to connect the spreadsheets", type=file, required=True, dest="json_file")
    parser.add_argument('-o', '--sample_sheet_outdir', help="Path where the sample sheets will be written in their respective project/run/lane folder", required=False, dest="outdir", default='/nb/Research/processingmgiscratch/sample_sheets')
    parser.add_argument('-p', '--processing_dir', help="Path where the MGI run processging will happen", required=False, dest="process_dir", default='/nb/Research/processingmgiscratch/processing')
    parser.add_argument('-r', '--genpipes_repo', help="Path of the GenPipes repository to use for generating the run processing pipeline scripts", required=True, dest="genpipes_repo")
    parser.add_argument('--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    MGI_spreadsheet_name = args.spreadsheet_name
    outdir = args.outdir
    process_dir = args.process_dir
    authentication_file = args.json_file.name
    genpipes_repo = args.genpipes_repo

    dict_of_columns = parse_google_sheet(
        MGI_spreadsheet_name,
        authentication_file
    )

    mgi_runs_file = os.path.join(os.path.dirname(authentication_file), ".mgi_runs.ref")
    compare_runs(
        dict_of_columns,
        outdir,
        process_dir,
        genpipes_repo,
        mgi_runs_file
    )
