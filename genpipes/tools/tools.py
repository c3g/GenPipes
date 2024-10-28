################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes Pipelines.
#
# GenPipes Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import os
import subprocess
import sys

from . import log_report
from . import job2json
from . import job2json_project_tracking

def add_subcommands(parser):
    """
    To be called by the main genpipes parser to add the tools subcommands.
    """
    tools_subparsers = parser.add_subparsers(title='tools subcommands', description='valid tools subcommands', help='additional tools help')

    # Create the parser for the "chunk_genpipes" subcommand
    parser_chunk_genpipes = tools_subparsers.add_parser('chunk_genpipes', help='Chunk Genpipes submitting script in multiple sub scripts with n jobs in them.')
    parser_chunk_genpipes.add_argument('-n', type=int, help='Maximum number of jobs in chunk')
    parser_chunk_genpipes.add_argument('genpipes_script', help="Genpipes output script.")
    parser_chunk_genpipes.add_argument('output_folder', help="Folder where to store chunks.")
    parser_chunk_genpipes.set_defaults(func=run_chunk_genpipes)

    # Create the parser for the "csvToreadset" subcommand
    parser_csvtoreadset = tools_subparsers.add_parser('csvToreadset', help='Take the csv project file downloadable from nanuq and creates the readset output tsv file needed to run the GenPipes.')
    parser_csvtoreadset.add_argument('input_csv', help="It can be downloaded from nanuq's Readset page (not Libraries page).")
    parser_csvtoreadset.add_argument('output_name', help="Name of the output readset file.")
    parser_csvtoreadset.add_argument('input_type', help="Type of input. Either fastq for bam.")
    parser_csvtoreadset.add_argument('data_path', help="Relative path to the data from the project folder.")
    parser_csvtoreadset.set_defaults(func=run_csvtoreadset)

    # Create the parser for the "job2json" subcommand
    parser_job2json = tools_subparsers.add_parser('job2json', help='Appends a JSON section describing a pipeline job that has just finished to a JSON file which was pre-generated when the pipeline was launched. This script is usually launched automatically before and after each pipeline job. /!\\ This version is for project trasking database only.')
    parser_job2json.add_argument('-s', '--step_name', required=True, help="name of the step of the current job")
    parser_job2json.add_argument('-j','--job_name', required=True, help="name of the current job")
    parser_job2json.add_argument('-l','--job_log', required=True, help="name of the log file for the current job")
    parser_job2json.add_argument('-d','--job_done', required=True, help="name of the done file for the current job")
    parser_job2json.add_argument('-o','--json_outfile', required=True, help="comma-separated list of names of json files which need to be appended by the current job")
    parser_job2json.add_argument('-f','--status', required=True, help="job status comming from GenPipes. Either running when job starts or $GenPipes_STATE")
    parser_job2json.add_argument('-u','--user', required=True, help="user to use GenPipes $USER")
    parser_job2json.set_defaults(func=run_job2json)

    # Create the parser for the "job2json_project_tracking" subcommand
    parser_job2json_project_tracking = tools_subparsers.add_parser('job2json_project_tracking', help='Appends a JSON section describing a pipeline job that has just finished to a JSON file which was pre-generated when the pipeline was launched. This script is usually launched automatically before and after each pipeline job. /!\\ This version is for project tracking database only.')
    parser_job2json_project_tracking.add_argument('-s', '--sample_names', required=True, help="comma-separated list of names of the samples of the current job")
    parser_job2json_project_tracking.add_argument('-r', '--readset_names', required=True, help="comma-separated list of names of the readsets of the current job")
    parser_job2json_project_tracking.add_argument('-j','--job_name', required=True, help="name of the current job")
    parser_job2json_project_tracking.add_argument('-m', '--metrics', required=False, help="comma-separated list of metrics of the current job: name=value,name=value,... With <name> = metric name; <value> = metric value")
    parser_job2json_project_tracking.add_argument('-o','--json_outfile', required=True, help="name of json output file")
    parser_job2json_project_tracking.add_argument('-f', '--status', required=False, help="status of job")
    parser_job2json_project_tracking.set_defaults(func=run_job2json_project_tracking)

    # Create the parser for the "log_report" subcommand
    parser_log_report = tools_subparsers.add_parser('log_report', help='Generate the log report for a GenPipes run')
    parser_log_report.add_argument('job_list_path', help="Path to a GenPipes job list")
    parser_log_report.add_argument('--remote', '-r', help="Remote HPC where the job was run", choices=['beluga', 'cedar', 'narval', 'abacus'], default=None)
    # As using sys.argv to check remote value making sure it's only in the context of a "genpipes tools log_report call"
    if any((remote in sys.argv for remote in ['--remote', '-r'])) and "tools" in sys.argv and "log_report" in sys.argv:
        if 'abacus' in sys.argv:
            parser_log_report.add_argument('--memtime', '-m', help="Output also memtime values if present in job output files", action='store_true', default=False)
            parser_log_report.add_argument('--success', '-s', help="Show successful jobs only", action='store_true', default=False)
            parser_log_report.add_argument('--nosuccess', '--nos', help="Show unsuccessful jobs only i.e. failed or uncompleted jobs", action='store_true', default=False)
        elif any((cluster in sys.argv for cluster in ['beluga', 'cedar', 'narval'])):
            parser_log_report.add_argument('--loglevel', help="Standard Python log level. Default: ERROR", choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"], default='ERROR')
            parser_log_report.add_argument('--tsv', help="Output to tsv file")
            parser_log_report.add_argument('--quiet', '-q', help="No report printed to terminal", action='store_true', default=False)
    parser_log_report.set_defaults(func=run_log_report)

    # Create the parser for the "submit_genpipes" subcommand
    parser_submit_genpipes = tools_subparsers.add_parser('submit_genpipes', help='Control the number of jobs in scheduler queue, and resubmit jobs if then fail at submit time.')
    parser_submit_genpipes.add_argument('-n', type=int, help='Maximum number of job in slurm queue. Default: 500', default=500)
    parser_submit_genpipes.add_argument('-s', type=int, help='Number of second to sleep when queue is full. Default: 120', default=120)
    parser_submit_genpipes.add_argument('-S', help='Scheduler running on the cluster. Default: slurm', choices=["slurm", "pbs"], default="slurm")
    parser_submit_genpipes.add_argument('-l', type=int, help='Will retry N time(s) to resubmit a chunk if error occurs. Default: 10', default=10)
    parser_submit_genpipes.add_argument('chunk_folder', help="The output folder from the chunk_genpipes.sh script.")
    parser_submit_genpipes.set_defaults(func=run_submit_genpipes)

    # Create the parser for the "get_wrapper" subcommand
    parser_get_wrapper = tools_subparsers.add_parser('get_wrapper', help='Get Genpipes In A Container image.')
    parser_get_wrapper.set_defaults(func=run_get_wrapper)

def run_chunk_genpipes(args):
    """
    Run the chunk_genpipes.sh script with the given arguments.
    """
    chunk_genpipes = os.path.join(os.path.dirname(__file__), 'chunk_genpipes.sh')
    cmd = [chunk_genpipes]
    if args.n:
        cmd += ['-n', str(args.n)]
    cmd += [args.genpipes_script]
    cmd += [args.output_folder]
    subprocess.run(cmd, check=False)

def run_csvtoreadset(args):
    """
    Run the csvToreadset.R script with the given arguments.
    """
    csvtoreadset = os.path.join(os.path.dirname(__file__), 'csvToreadset.R')
    cmd = ['Rscript', csvtoreadset, args.input_csv, args.output_name, args.input_type, args.data_path]
    subprocess.run(cmd, check=False)

def run_job2json(args):
    """
    Run the job2json.py script with the given arguments.
    """
    job2json.main(args)

def run_job2json_project_tracking(args):
    """
    Run the job2json_project_tracking.py script with the given arguments.
    """
    job2json_project_tracking.main(args)

def run_log_report(args):
    """
    Run the log_report.pl script with the given arguments.
    """
    if args.remote == 'abacus':
        # Call log_report.pl
        log_report_pl = os.path.join(os.path.dirname(__file__), 'log_report.pl')
        cmd = [log_report_pl, args.job_list_path]
        if args.memtime:
            cmd += ['--memtime']
        if args.success:
            cmd += ['--success']
        if args.nosuccess:
            cmd += ['--nosuccess']
        if args.tsv:
            cmd += ['--tsv', args.tsv]
        if args.quiet:
            cmd += ['--quiet']
        subprocess.run(cmd, check=False)
    else:
        # Call the main function of log_report.py
        log_report.main(args)

def run_submit_genpipes(args):
    """
    Run the submit_genpipes.sh script with the given arguments.
    """
    submit_genpipes = os.path.join(os.path.dirname(__file__), 'submit_genpipes.sh')
    cmd = [submit_genpipes]
    if args.n:
        cmd += ['-n', str(args.n)]
    if args.s:
        cmd += ['-s', str(args.s)]
    if args.S:
        cmd += ['-S', args.S]
    if args.l:
        cmd += ['-l', str(args.l)]
    cmd += [args.chunk_folder]
    subprocess.run(cmd, check=False)

def run_get_wrapper():
    """
    Run the get_wrapper.sh script.
    """
    get_wrapper = os.path.join(os.path.dirname(__file__), 'get_wrapper.sh')
    cmd = [get_wrapper]
    subprocess.run(cmd, check=False)
