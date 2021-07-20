#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import csv
import glob
import http.client
import logging
import os
import re

log = logging.getLogger(__name__)

def get_nanuq_file(nanuq_auth_file, nanuq_url, nanuq_file):

    if os.path.exists(nanuq_file):
        log.warning("File " + nanuq_file + " already exists! Skipping...")
    else:
        https_connection = http.client.HTTPSConnection("genomequebec.mcgill.ca")
        https_connection.set_debuglevel(1)
        headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/plain"}

        log.info("Fetching Nanuq file from server...")

        with open(nanuq_auth_file) as auth_file:
            https_connection.request("POST", nanuq_url, auth_file, headers)

        http_response = https_connection.getresponse()

        log.info("HTTP Response Status: " + str(http_response.status))
        log.info("HTTP Response Reason: " + http_response.reason)

        with open(nanuq_file, 'w') as file:
            file.write(http_response.read())

        https_connection.close()


def get_nanuq_readset_file(nanuq_auth_file, nanuq_project_id, nanuq_readset_file, seq_type):
    get_nanuq_file(nanuq_auth_file, "/nanuqMPS/csv/technology/" + seq_type + "/project/" + nanuq_project_id + "/filename/" + nanuq_readset_file, nanuq_readset_file)


def get_nanuq_bed_file(nanuq_auth_file, bed_file):
    get_nanuq_file(nanuq_auth_file, "/nanuqLimsCgi/targetRegion/downloadBed.cgi?bedName=" + bed_file, bed_file)


def create_readsets(nanuq_readset_file, seq_type, mugqic_pipelines_readset_file="readsets.tsv", args_nanuq_auth_file=None, args_run=None):
    if seq_type != "NovaSeq":
        # Lowercase the first seq_type character
        lcfirst_seq_type = seq_type[0].lower() + seq_type[1:]
    else:
        lcfirst_seq_type = seq_type

    nanuq_readset_root_directory = "/lb/robot/" + lcfirst_seq_type + "Sequencer/" + lcfirst_seq_type + "Runs"
    raw_reads_directory = "raw_reads"
    symlinks = []
    mugqic_pipelines_readset_csv_rows = []
    bed_files = set()  # Only used for HiSeq/MiSeq sequencing type

    # Parse Nanuq readset file and list symlinks to be created
    log.info("Parse Nanuq readset file " + nanuq_readset_file + " ...")
    nanuq_readset_csv = csv.DictReader(open(nanuq_readset_file, 'rb'), delimiter=',', quotechar='"')

    for line in nanuq_readset_csv:
        if line['Status'] and line['Status'] == "Data is valid":
            mugqic_pipelines_readset_csv_row = {}

            if (not args_run) or (line['Run'] in args_run):
                log.info("Parsing Nanuq run " + line['Run'] + " ...")

                if seq_type == "Pacbio":
                    mugqic_pipelines_readset_csv_row['Sample'] = line['Sample Group'] if 'Sample Group' in line and line['Sample Group'] != "" else line['Name']
                    mugqic_pipelines_readset_csv_row['Readset'] = ".".join([line['Name'], line['Library Barcode'], line['Run'], line['Well']])

                    nanuq_vs_mugqic_pipelines_readset_keys = [
                        ['Run', 'Run'],
                        ['Well', 'Smartcell'],
                        ['Collection Protocol', 'Protocol']
                    ]
                    formats = ['BAS', 'BAX']

                    fieldnames = ['Sample', 'Readset'] + [key[1] for key in nanuq_vs_mugqic_pipelines_readset_keys] + ['NbBasePairs', 'EstimatedGenomeSize'] + formats

                    nb_basepairs = re.search("^\([^/]*/[^/]*/(.*)\)$", line['Longest Subreads (count mean bp)'])
                    mugqic_pipelines_readset_csv_row['NbBasePairs'] = re.sub(",", "", nb_basepairs.group(1))
                    mugqic_pipelines_readset_csv_row['EstimatedGenomeSize'] = line['Genome size from Nanuq reception (Mb)']
                    log.warning('EstimatedGenomeSize: ' + line['Genome size from Nanuq reception (Mb)'])

                    if mugqic_pipelines_readset_csv_row['EstimatedGenomeSize'] != '':
                        mugqic_pipelines_readset_csv_row['EstimatedGenomeSize'] = int(float(mugqic_pipelines_readset_csv_row['EstimatedGenomeSize'])*1000*1000)
                    if line.get('Results Directory', None):
                        nanuq_readset_prefix = os.path.normpath(os.path.join(nanuq_readset_root_directory, line['Results Directory'], line['Movie name']))
                        for format in formats:
                            nanuq_readset_paths = sorted(glob.glob(nanuq_readset_prefix + "*." + format.lower() + ".h5"))
                            mugqic_pipelines_readset_paths = [os.path.join(raw_reads_directory, line['Name'], os.path.basename(nanuq_readset_path)) for nanuq_readset_path in nanuq_readset_paths]
                            mugqic_pipelines_readset_csv_row[format] = ",".join(mugqic_pipelines_readset_paths)
                            for nanuq_readset_path, mugqic_pipelines_readset_path in zip(nanuq_readset_paths, mugqic_pipelines_readset_paths):
                                symlinks.append([nanuq_readset_path, mugqic_pipelines_readset_path])

                else:  # seq_type = HiSeq or MiSeq or NovaSeq or iSeq
                    mugqic_pipelines_readset_csv_row['Sample'] = line['Sample Group'] if 'Sample Group' in line and line['Sample Group'] != "" else line['Name']
                    mugqic_pipelines_readset_csv_row['Readset'] = ".".join([line['Name'], line['Library Barcode'], line['Run'], line['Region']])

                    nanuq_vs_mugqic_pipelines_readset_keys = [
                        ['Library Barcode', 'Library'],
                        ['Run Type', 'RunType'],
                        ['Run', 'Run'],
                        ['Region', 'Lane'],
                        ['Adaptor Read 1 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)', 'Adapter1'],
                        ['Adaptor Read 2 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)', 'Adapter2'],
                        ['Forward Primer Sequence', 'Primer1'],
                        ['Reverse Primer Sequence', 'Primer2'],
                        ['Quality Offset', 'QualityOffset'],
                        ['BED Files', 'BED']
                    ]
                    formats = ['FASTQ1', 'FASTQ2', 'BAM']

                    fieldnames = ['Sample', 'Readset'] + [key[1] for key in nanuq_vs_mugqic_pipelines_readset_keys] + formats

                    for format in formats:
                        if line.get(format, None):
                            nanuq_readset_path = os.path.normpath(os.path.join(nanuq_readset_root_directory, line[format]))
                            if os.path.isfile(nanuq_readset_path) and (format != 'FASTQ2' or line['Run Type'] != 'SINGLE_END'):  # Ignore FASTQ2 value if any, in case of SINGLE_END readset
                                if format == 'BAM':
                                    mugqic_pipelines_readset_basename = mugqic_pipelines_readset_csv_row['Readset'] + ".bam"
                                else:  # format = FASTQ1 or FASTQ2
                                    if line['Run Type'] == 'PAIRED_END':
                                        mugqic_pipelines_readset_basename = mugqic_pipelines_readset_csv_row['Readset'] + ".pair" + format[-1] + ".fastq.gz"
                                    else:  # format = FASTQ1 and Run Type = SINGLE_END
                                        mugqic_pipelines_readset_basename = mugqic_pipelines_readset_csv_row['Readset'] + ".single.fastq.gz"

                                mugqic_pipelines_readset_path = os.path.join(raw_reads_directory, mugqic_pipelines_readset_csv_row['Sample'], mugqic_pipelines_readset_basename)
                                symlinks.append([nanuq_readset_path, mugqic_pipelines_readset_path])
                                mugqic_pipelines_readset_csv_row[format] = mugqic_pipelines_readset_path

                                # Add BAM index to symlinks if it exists, log a warning otherwise
                                if format == 'BAM':
                                    nanuq_readset_index_path = re.sub("\.bam$", ".bai", nanuq_readset_path)
                                    if os.path.isfile(nanuq_readset_index_path):
                                        symlinks.append([nanuq_readset_index_path, re.sub("\.bam$", ".bai", mugqic_pipelines_readset_path)])
                                    else:
                                        log.warning("Nanuq readset index path " + nanuq_readset_index_path + " is invalid!")
                            else:
                                raise Exception("Error: Nanuq readset path " + nanuq_readset_path + " is invalid!")

                    # BED files
                    # Filter empty strings returned by split with string ";" separator
                    for bed_file in filter(None, line['BED Files'].split(';')):
                        # Retrieve BED file if not previously done
                        if bed_file not in bed_files:
                            if args_nanuq_auth_file:
                                get_nanuq_bed_file(args_nanuq_auth_file.name, bed_file)
                            else:
                                log.warning("Nanuq authentication file missing: skipping retrieval of " + bed_file + "...")
                            bed_files.add(bed_file)

                for nanuq_key, mugqic_pipelines_key in nanuq_vs_mugqic_pipelines_readset_keys:
                    value = line.get(nanuq_key, None)
                    if value:
                        mugqic_pipelines_readset_csv_row[mugqic_pipelines_key] = value

                mugqic_pipelines_readset_csv_rows.append(mugqic_pipelines_readset_csv_row)

            else:
                log.debug(str(line) + " line data is not in requested runs... skipping")
        else:
            log.warning(str(line) + " line data is not valid... skipping")

    # Create symbolic links and parent directories if necessary
    for target, link_name in symlinks:
        if os.path.islink(link_name):
            log.warning("Symlink " + link_name + " already exists! Skipping...")
        else:
            link_directory = os.path.dirname(link_name)
            if not os.path.isdir(link_directory):
                os.makedirs(link_directory)
            log.info("Creating symlink " + link_name + " ...")
            os.symlink(target, link_name)
            log.info("Symlink " + link_name + " created successfully.")

    # Write mugqic_pipelines readset file if necessary
    if os.path.exists(mugqic_pipelines_readset_file):
        log.warning("File " + mugqic_pipelines_readset_file + " already exists! Skipping...")
    else:
        mugqic_pipelines_readset_csv = csv.DictWriter(open(mugqic_pipelines_readset_file, 'wb'), fieldnames=fieldnames, delimiter='\t')
        try:
            mugqic_pipelines_readset_csv.writeheader()
        except AttributeError as ae:
            log.error("This script requires at minimum Python 2.7/3.2.")
        mugqic_pipelines_readset_csv.writerows(mugqic_pipelines_readset_csv_rows)

#-------------------------------------------------------------------------------
# Main script

# Parse options
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-p", "--nanuq-project-id", help="Nanuq project ID used to fetch readset file from server (incompatible with --nanuq-readset-file)")
group.add_argument("-r", "--nanuq-readset-file", help="Nanuq readset file to use instead of fetching it from server (incompatible with --nanuq-project-id)", type=file)
parser.add_argument("-rn", "--run", help="Nanuq run ID to fetch only readset procesed in specified run(s). Multiple runs can be provided, should they be comma-separated")
parser.add_argument("-s", "--seq-type", help="Sequencing type (default: HiSeq)", choices=["HiSeq", "NovaSeq", "MiSeq", "iSeq", "Pacbio"], default="HiSeq")
parser.add_argument("-a", "--nanuq-auth-file", help="Nanuq authentication file containing your Nanuq username and password e.g. $HOME/.nanuqAuth.txt\nTo create it:\n$ echo -n \"user=<USERNAME>&password=<PASSWORD>\" > $HOME/.nanuqAuth.txt ; chmod u+r,go-rwx $HOME/.nanuqAuth.txt\nNote '-n' option since trailing newline is not allowed at the end of the file.", type=file)
parser.add_argument("-nl", "--no-links", help="Do not create raw_reads directory and symlinks (default: false)", action="store_true")
parser.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")

args = parser.parse_args()

logging.basicConfig(level=getattr(logging, args.log.upper()))

if args.nanuq_project_id:
    if not args.nanuq_auth_file:
        raise Exception("Error: missing Nanuq authentication file (use -a or --nanuq-auth-file)!")

    nanuq_readset_file = "project.nanuq." + args.nanuq_project_id + ".csv"
    mugqic_pipelines_readset_file = "readsets." + args.seq_type + "." + args.nanuq_project_id + ".tsv"

    get_nanuq_readset_file(args.nanuq_auth_file.name, args.nanuq_project_id, nanuq_readset_file, args.seq_type)

elif args.nanuq_readset_file:
    nanuq_readset_file = args.nanuq_readset_file.name
    mugqic_pipelines_readset_file = "readsets."+ args.seq_type +".tsv"

if not args.no_links:
    create_readsets(nanuq_readset_file, args.seq_type, mugqic_pipelines_readset_file, args.nanuq_auth_file, args.run)
