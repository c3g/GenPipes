#!/usr/bin/env python

import argparse
import csv
import glob
import logging
import os
import re
import subprocess

log = logging.getLogger(__name__)

def get_nanuq_readset_file(nanuq_auth_file, nanuq_project_id, nanuq_readset_file, seq_type):

    command = "wget --no-cookies --post-file " + nanuq_auth_file + " https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/" + seq_type + "/project/" + nanuq_project_id + "/filename/" + nanuq_readset_file

    log.info("Fetching Nanuq readset file from server...")
    log.info(command)

    subprocess.check_call(["bash", "-c", command])

def get_nanuq_bed_files(nanuq_auth_file, bed_files):
    for bed in bed_files:
        command = "wget --no-cookies --post-file " + nanuq_auth_file + " https://genomequebec.mcgill.ca/nanuqLimsCgi/targetRegion/downloadBed.cgi?bedName=" + bed + " -O " + bed
        log.info("Fetching Nanuq BED file from server...")
        log.info(command)

        subprocess.check_call(["bash", "-c", command])

def create_readsets(nanuq_auth_file, nanuq_readset_file, seq_type, mugqic_pipeline_readset_file="readsets.tsv"):
    # Lowercase the first seq_type character
    lcfirst_seq_type = seq_type[0].lower() + seq_type[1:]

    nanuq_readset_root_directory = "/lb/robot/" + lcfirst_seq_type + "Sequencer/" + lcfirst_seq_type + "Runs"
    raw_reads_directory = "raw_reads"
    symlinks = []
    mugqic_pipeline_readset_csv_rows = []

    # Parse Nanuq readset file and list symlinks to be created
    log.info("Parse Nanuq readset file " + nanuq_readset_file + " ...")
    nanuq_readset_csv = csv.DictReader(open(nanuq_readset_file, 'rb'), delimiter=',', quotechar='"')
    bed_files = set();
    for line in nanuq_readset_csv:
        if line['Status'] and line['Status'] == "Data is valid":
            mugqic_pipeline_readset_csv_row = {}

            if seq_type == "Pacbio":
                nanuq_vs_mugqic_pipeline_readset_keys = [
                    ['Name', 'Sample'],
                    ['Filename Prefix', 'Readset'],
                    ['Run', 'Run'],
                    ['Well', 'Smartcell'],
                    ['Collection Protocol', 'Protocol']
                ]
                formats = ['BAS', 'BAX']

                fieldnames = [key[1] for key in nanuq_vs_mugqic_pipeline_readset_keys] + ['NbBasePairs', 'EstimatedGenomeSize'] + formats

                nb_basepairs = re.search("^\([^/]*/[^/]*/(.*)\)$", line['Longest Subreads (count mean bp)'])
                mugqic_pipeline_readset_csv_row['NbBasePairs'] = re.sub(",", "", nb_basepairs.group(1))

                if line.get('Results Directory', None):
                    nanuq_readset_prefix = os.path.normpath(os.path.join(nanuq_readset_root_directory, line['Results Directory'], line['Movie name']))
                    for format in formats:
                        nanuq_readset_paths = sorted(glob.glob(nanuq_readset_prefix + "*." + format.lower() + ".h5"))
                        mugqic_pipeline_readset_paths = [os.path.join(raw_reads_directory, line['Name'], os.path.basename(nanuq_readset_path)) for nanuq_readset_path in nanuq_readset_paths]
                        mugqic_pipeline_readset_csv_row[format] = ",".join(mugqic_pipeline_readset_paths)
                        for nanuq_readset_path, mugqic_pipeline_readset_path in zip(nanuq_readset_paths, mugqic_pipeline_readset_paths):
                            symlinks.append([nanuq_readset_path, mugqic_pipeline_readset_path])

            else:  # seq_type = HiSeq or MiSeq
                nanuq_vs_mugqic_pipeline_readset_keys = [
                    ['Name', 'Sample'],
                    ['Filename Prefix', 'Readset'],
                    ['Library Barcode', 'Library'],
                    ['Run Type', 'RunType'],
                    ['Run', 'Run'],
                    ['Region', 'Lane'],
                    ['Adaptor Read 1 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)', 'Adaptor1'],
                    ['Adaptor Read 2 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)', 'Adaptor2'],
                    ['Quality Offset', 'QualityOffset'],
                    ['BED Files', 'BED']
                ]
                formats = ['FASTQ1', 'FASTQ2', 'BAM']

                fieldnames = [key[1] for key in nanuq_vs_mugqic_pipeline_readset_keys] + formats
                # Filter empty strings returned by split with string ";" separator
                available_beds = filter(None, line['BED Files'].split(';'))
                for bed in available_beds:
                    bed_files.add(bed)

                for format in formats:
                    if line.get(format, None):
                        nanuq_readset_path = os.path.normpath(os.path.join(nanuq_readset_root_directory, line[format]))
                        if os.path.isfile(nanuq_readset_path):
                            mugqic_pipeline_readset_path = os.path.join(raw_reads_directory, line['Name'], os.path.basename(nanuq_readset_path))
                            symlinks.append([nanuq_readset_path, mugqic_pipeline_readset_path])
                            mugqic_pipeline_readset_csv_row[format] = mugqic_pipeline_readset_path
                        else:
                            raise Exception("Error: Nanuq readset path " + nanuq_readset_path + " is invalid!")

            for nanuq_key, mugqic_pipeline_key in nanuq_vs_mugqic_pipeline_readset_keys:
                value = line.get(nanuq_key, None)
                if value:
                    mugqic_pipeline_readset_csv_row[mugqic_pipeline_key] = value

            mugqic_pipeline_readset_csv_rows.append(mugqic_pipeline_readset_csv_row)
        else:
            log.warning(str(line) + " line data is not valid... skipping")

    # Get the bed files once
    get_nanuq_bed_files(nanuq_auth_file, bed_files)

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

    mugqic_pipeline_readset_csv = csv.DictWriter(open(mugqic_pipeline_readset_file, 'wb'), fieldnames=fieldnames, delimiter='\t')
    mugqic_pipeline_readset_csv.writeheader()
    mugqic_pipeline_readset_csv.writerows(mugqic_pipeline_readset_csv_rows)

#-------------------------------------------------------------------------------
# Main script

# Parse options
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-p", "--nanuq-project-id", help="Nanuq project ID used to fetch readset file from server (incompatible with --nanuq-readset-file)")
group.add_argument("-r", "--nanuq-readset-file", help="Nanuq readset file to use instead of fetching it from server (incompatible with --nanuq-project-id)", type=file)

parser.add_argument("-s", "--seq-type", help="Sequencing type (default: HiSeq)", choices=["HiSeq", "MiSeq", "Pacbio"], default="HiSeq")
parser.add_argument("-a", "--nanuq-auth-file", help="Nanuq authentication file containing your Nanuq username and password e.g. $HOME/.nanuqAuth.txt\nTo create it:\n$ echo -n \"user=<USERNAME>&password=<PASSWORD>\" > $HOME/.nanuqAuth.txt ; chmod u+r,go-rwx $HOME/.nanuqAuth.txt\nNote '-n' option since trailing newline is not allowed at the end of the file.", type=file)
parser.add_argument("-nl", "--no-links", help="Do not create raw_reads directory and symlinks (default: false)", action="store_true")
parser.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")

args = parser.parse_args()

logging.basicConfig(level=getattr(logging, args.log.upper()))

if args.nanuq_project_id:
    if not args.nanuq_auth_file:
        raise Exception("Error: missing Nanuq authentication file (use -a or --nanuq-auth-file)!")

    nanuq_readset_file = "project.nanuq." + args.nanuq_project_id + ".csv"
    mugqic_pipeline_readset_file = "readsets." + args.nanuq_project_id + ".tsv"

    get_nanuq_readset_file(args.nanuq_auth_file.name, args.nanuq_project_id, nanuq_readset_file, args.seq_type)

elif args.nanuq_readset_file:
    nanuq_readset_file = args.nanuq_readset_file.name
    mugqic_pipeline_readset_file = "readsets.tsv"

if not args.no_links:
    create_readsets(args.nanuq_auth_file.name, nanuq_readset_file, args.seq_type, mugqic_pipeline_readset_file)
