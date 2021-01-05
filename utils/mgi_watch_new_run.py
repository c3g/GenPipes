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

def get_date_from_run(
    columns,
    run
    ):
    for x in range(len(columns['RUN_ID'])):
        if columns['RUN_ID'][x] == run:
            if columns['Run_Date'][x]:
                return columns['Run_Date'][x]

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

def new_naming_convention(
    run_date
    ):

    is_new_naming = True

    [year, month, day] = run_date.split("-")
    if int(year) == 2020:
        if int(month) == 12:
            if int(day) < 1:
                is_new_naming = False
        elif int(month) < 12:
            is_new_naming = False

    return is_new_naming


def compare_runs(
    columns,
    outdir,
    genpipes_scr_dir,
    process_dir,
    mgi_runs_file,
    run_id,
    is_demultiplexed,
    extra_options
    ):

    header = "RUN_ID"
    if os.path.isfile(mgi_runs_file):
        new_list = filter(None, sorted(set(columns['RUN_ID'])))
#        ref_list = []
        ref_list = filter(None, open(mgi_runs_file, "r").read().split("\n"))

        if run_id or new_list != ref_list:
            if run_id:
                runs_added = [run_id]
            else:
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
                print flowcell
                if flowcell:
                    run_date = get_date_from_run(
                        columns,
                        run
                    )
                    if new_naming_convention(run_date):
                        run_folder_basename = flowcell + "_" + run
                    else:
                        run_folder_basename = flowcell

                    print "Processing run " + run
                    project = get_project_from_run(
                        columns,
                        run
                    )
                    print project
                    print "    Generating sample sheet..."
                    print_sample_sheet(
                        columns,
                        project,
                        run,
                        None,
                        outdir,
                        group=True
                    )
                    print "    Generating GenPipes script..."
                    print_genpipes_scripts(
                        outdir,
                        process_dir,
                        genpipes_scr_dir,
                        project,
                        run,
                        None,
                        sequencer_path,
                        run_folder_basename,
                        is_demultiplexed,
                        extra_options
                    )

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
                            outdir,
                            group=False
                        )
                        print "        Generating GenPipes script..."
                        print_genpipes_scripts(
                            outdir,
                            process_dir,
                            genpipes_scr_dir,
                            project,
                            run,
                            lane,
                            sequencer_path,
                            run_folder_basename,
                            is_demultiplexed,
                            extra_options
                        )
            # replace referece run list by the current run list to set it as the reference for next watch round
            print_runs(
                columns,
                mgi_runs_file
            )

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
    genpipes_scr_dir,
    project,
    run,
    lane,
    sequencer_path,
    run_folder_basename,
    is_demultiplexed,
    extra_options
    ):

    if lane:
        sample_sheet = os.path.join(samplesheet_dir, project, run, project+"."+run+".L0"+lane+".sample_sheet.csv")
    else:
        sample_sheet = os.path.join(samplesheet_dir, project, run, project+"."+run+".sample_sheet.csv")

    sample_sheet_rows = [row for row in csv.DictReader(open(sample_sheet, 'rb'), delimiter=',')]
    print len(sample_sheet_rows)
    print sample_sheet_rows[0]['Index']

    if not os.path.exists(os.path.join(genpipes_scr_dir, project, run)):
        os.makedirs(os.path.join(genpipes_scr_dir, project, run))

    if lane:
        genpipes_script = open(os.path.join(genpipes_scr_dir, project, run, project + "." + run + ".L0" + lane + ".genpipes_script.sh"), 'wb+')
    else:
        genpipes_script = open(os.path.join(genpipes_scr_dir, project, run, project + "." + run + ".genpipes_script.sh"), 'wb+')

    genpipes_script.write("""\
module load mugqic/python/2.7.14 mugqic_dev/genpipes/3.1.6 && \\
mkdir -p {process_dir}/{process_dir_suffix} && \\
python $MUGQIC_PIPELINES_HOME/pipelines/mgi_run_processing/mgi_run_processing.py \\
  -c $MUGQIC_PIPELINES_HOME/pipelines/mgi_run_processing/mgi_run_processing.base.ini $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \\
  --no-json -l debug {demux_fastq}{lane}{extra_options} \\
  -d /nb/Research/MGISeq/{sequencer_path}/{run_folder_basename} \\
  -r {samplesheet_dir}/{process_dir_suffix}/{outfile_prefix}.sample_sheet.csv \\
  -o {process_dir}/{process_dir_suffix} \\
  > {process_dir}/{project}/{run}/{outfile_prefix}.sh \\
  2> {process_dir}/{project}/{run}/{outfile_prefix}.trace.log""".format(
        process_dir=process_dir,
        process_dir_suffix=project+"/"+run,
        outfile_prefix=project+"."+run+".L0"+lane if lane else project+"."+run,
        samplesheet_dir=samplesheet_dir,
        demux_fastq="--demux-fastq " if is_demultiplexed else "",
        project=project,
        run_folder_basename=run_folder_basename,
        run=run,
        lane="--lane "+lane+" " if lane else "",
        sequencer_path=sequencer_path,
        extra_options=extra_options
    ))
    genpipes_script.close()

    if lane:
        subprocess.call("bash " + os.path.join(genpipes_scr_dir, project, run, project + "." + run + ".L0" + lane + ".genpipes_script.sh"), shell=True)
    else:
        subprocess.call("bash " + os.path.join(genpipes_scr_dir, project, run, project + "." + run + ".genpipes_script.sh"), shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--spreadsheet_name', help='Name of the MGI run management Google spreadsheet to parse', required=False, dest="spreadsheet_name", default='ALL MGI Run Management')
    parser.add_argument('-a', '--authentication_file', help="JSON authentication file used to connect the spreadsheets", type=file, required=True, dest="json_file")
    parser.add_argument('-o', '--sample_sheet_outdir', help="Path where the sample sheets will be written in their respective project/run/lane subfolder", required=False, dest="outdir", default='/nb/Research/processingmgiscratch/sample_sheets')
    parser.add_argument('-g', '--genpipes_scripts_outdir', help="Path where the genpipes scripts will be written and executed in their respective project/run/lane subfolder", required=False, dest="genpipes_scr_dir", default='/nb/Research/processingmgiscratch/genpipes_scripts')
    parser.add_argument('-p', '--processing_dir', help="Path where the MGI run processging will happen", required=False, dest="process_dir", default='/nb/Research/processingmgiscratch/processing')
    parser.add_argument('-r', '--run', help="RUN ID : sample sheets and genpipes_scripts will only be created for the specified RUN ID", required=False, dest="run_id")
    parser.add_argument('-d', '--demultiplexed-fastqs', help="Fastqs of the run are expected to be demultiplexed : will prepare genpipes scripts without --raw-fastq flag", action="store_true", required=False, dest="is_demultiplexed")
    parser.add_argument('-e', '--extra-options', help="Extra parameters to include in the genpipes command - mostly used for tests and corner cases...", required=False, dest="extra_options", default="")
    parser.add_argument('--loglevel', help="Standard Python log level", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    MGI_spreadsheet_name = args.spreadsheet_name
    outdir = args.outdir
    genpipes_scr_dir = args.genpipes_scr_dir
    process_dir = args.process_dir
    authentication_file = args.json_file.name
    run_id = args.run_id
    if args.is_demultiplexed:
        is_demultiplexed = True
    else:
        is_demultiplexed = False
    extra_options = args.extra_options

    dict_of_columns = parse_google_sheet(
        MGI_spreadsheet_name,
        authentication_file
    )

    mgi_runs_file = os.path.join(os.path.dirname(authentication_file), ".mgi_runs.ref")
    compare_runs(
        dict_of_columns,
        outdir,
        genpipes_scr_dir,
        process_dir,
        mgi_runs_file,
        run_id,
        is_demultiplexed,
        extra_options
    )
