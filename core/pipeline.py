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

# Check Python version
import sys

import utils.utils

if sys.version_info < (3,):
    raise SystemExit("Incompatible Python version: " + sys.version + "\nPython 3.X is required")

# Python Standard Modules
import argparse
import collections
import datetime
import hashlib
import logging
import os
import re
import subprocess
import textwrap
from uuid import uuid4

# MUGQIC Modules
from .config import _raise, SanitycheckError, global_config_parser
from .job import Job
from .scheduler import create_scheduler
from .step import Step

from bfx import jsonator

log = logging.getLogger(__name__)

class Pipeline(object):
    # Pipeline version
    version = None
    _argparser = None

    def __init__(self, config_files, sanity_check=False,
                 output_dir=None, job_scheduler=None, container=None, genpipes_file=None, no_json=False,
                 steps=None, report=False, clean=False, force=False):

        self._timestamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

        self.sanity_check = sanity_check
        if sanity_check:
            global_config_parser.sanity = True

        global_config_parser.parse_files(config_files)
        # Create a config trace from merged config file values
        with open(self.__class__.__name__ + ".config.trace.ini", 'w') as config_trace:
            config_trace.write(textwrap.dedent("""\
              # {self.__class__.__name__} Config Trace
              # Created on: {self.timestamp}
              # From:
              #   {config_files}
              # DO NOT EDIT THIS AUTOMATICALLY GENERATED FILE - edit the master config files

            """).format(config_files="\n#   ".join([config_file.name for config_file in config_files]), self=self))
            global_config_parser.write(config_trace)
            global_config_parser._filepath = os.path.abspath(config_trace.name)

        if output_dir:
            self._output_dir = output_dir
        else:
            self._output_dir = os.getcwd()

        if job_scheduler is not None:
            self._job_scheduler = create_scheduler(job_scheduler, config_files,
                                               container=container, genpipes_file=genpipes_file)
        # sorry about the double negative...
        self._json = not no_json

        step_counter = collections.Counter(self.step_list)
        duplicated_steps = [step.__name__ for step in step_counter if step_counter[step] > 1]
        if duplicated_steps:
            raise ValueError("Error: pipeline contains duplicated steps: " + ", "
                             .join(duplicated_steps) + "!")

        self._step_to_execute = None
        self.steps = steps
        if self.steps:
            if re.search("^\d+([,-]\d+)*$", steps):
                self._step_to_execute = [Step(self.step_list[i - 1]) for i in parse_range(steps)]
            else:
                raise Exception("Error: step range \"" + steps +
                    "\" is invalid (should match \d+([,-]\d+)*)!")
        else:
            log.warning("No step provided by the user => creating the entire pipeline\n")
            self._step_to_execute = [Step(step) for step in self.step_list]

        self._sample_list = []
        self._sample_paths = []

        # For job reporting, all jobs must be created first, no matter whether they are up to date or not
        if report:
            self._force_jobs = True
            self.create_jobs()
            self.report_jobs()
        # For job cleaning, all jobs must be created first, no matter whether they are up to date or not
        elif clean:
            self._force_jobs = True
            self.create_jobs()
            self.clean_jobs()
        elif sanity_check:
            self._force_jobs = force
            self.create_jobs()
        else:
            try:
                self._force_jobs = force
                self.create_jobs()
            except SanitycheckError as e:
                log.error("""
***The pipeline encounterered an error :
    {error}
***Please try running the pipeline in SANITY CHECK mode using the '--sanity-check' flag to check for more potential issues...""".format(
                   error=e
                   ))
                exit(1)
            
    @classmethod
    def process_help(cls, argv):
        """
        Genpipes help info
        Returns: the epilog

        """
        step_lst = []
        steps_doc = []
        for protocol, steps in cls.protocols(cls).items():
            step_lst.append('Protocol {}'.format(protocol))
            for i, step in enumerate(steps):
                step_lst.append('{} {}'.format(i, step.__name__))
                steps_doc.append('{} \n{}\n {}'.format(step.__name__, "-" * len(step.__name__),
                                                       textwrap.dedent(step.__doc__) if step.__doc__ else ""))
        epilog = "\n".join(step_lst)
        if '--help' in argv:
            epilog = epilog + "\n".join(steps_doc)
        return epilog

    @classmethod
    def argparser(cls, argparser):

            cls._argparser = argparser

            cls._argparser.add_argument("-c", "--config", help="config INI-style list of files; config parameters are "
                                                                "overwritten based on files order",
                                        nargs="+", type=argparse.FileType('r'), required=True)

            cls._argparser.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'")

            cls._argparser.add_argument("-o", "--output-dir", help="output directory (default: current)",
                                        default=os.getcwd())

            cls._argparser.add_argument("-j", "--job-scheduler", help="job scheduler type (default: slurm)",
                                        choices=["pbs", "batch", "daemon", "slurm"], default="slurm")

            cls._argparser.add_argument("-f", "--force", help="force creation of jobs even if up to date "
                                                              "(default: false)", action="store_true")

            cls._argparser.add_argument("--no-json", help="do not create JSON file per analysed sample to track the "
                                                          "analysis status (default: false i.e. JSON file will be "
                                                          "created)", action="store_true")

            cls._argparser.add_argument("--report", help="create 'pandoc' command to merge all job markdown report "
                                                         "files in the given step range into HTML, if they exist; if "
                                                         "--report is set, --job-scheduler, --force, --clean options "
                                                         "and job up-to-date status are ignored (default: false)",
                                        action="store_true")

            cls._argparser.add_argument("--clean", help="create 'rm' commands for all job removable files in the given"
                                                        " step range, if they exist; if --clean is set,"
                                                        " --job-scheduler, --force options and job up-to-date status "
                                                        "are ignored (default: false)", action="store_true")
            cls._argparser.add_argument("--container", nargs=2,
                                        help="Run inside a container providing a valid"
                                        "singularity image path", action=ValidateContainer,
                                        metavar=("{wrapper, singularity}",
                                                 "<IMAGE PATH>"))

            cls._argparser.add_argument("--genpipes_file", '-g', default=sys.stdout, type=argparse.FileType('w'),
                                        help="Command file output path. This is the command used to process "
                                             "the data, or said otherwise, this command will \"run the "
                                             "Genpipes pipeline\". Will be redirected to stdout if the "
                                             "option is not provided.")

            cls._argparser.add_argument("-l", "--log", help="log level (default: info)",
                                        choices=["debug", "info", "warning", "error", "critical"], default="info")

            cls._argparser.add_argument("--sanity-check", help="run the pipeline in `sanity check mode` to verify that"
                                                               " all the input files needed for the pipeline to run "
                                                               "are available on the system (default: false)",
                                        action="store_true")

            cls._argparser.add_argument('--wrap', type=str, help="Path to the genpipe cvmfs wrapper script.\n"
                                        "Default is genpipes/ressources/container/bin/container_wrapper.sh. This is "
                                        "a convenience options for using genpipes in a container",
                                        nargs='?')

            return cls._argparser

    # Pipeline command line arguments
    @property
    def args(self):
        return self._args

    @property
    def timestamp(self):
        return self._timestamp

    @property
    def output_dir(self):
        return self._output_dir

    @property
    def report_template_dir(self):
        return os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))
                            , "bfx", "report")

    @property
    def job_scheduler(self):
        return self._job_scheduler

    @property
    def force_jobs(self):
        return self._force_jobs

    @property
    def protocol(self):
        return self._protocol

    @protocol.setter
    def protocol(self, val):
        self._protocol = val


    @property
    def protocols(self):
        raise NotImplementedError("protocols method needs to be defined in the pipeline child class")

    @property
    def step_list(self):
        return self.protocols[self.protocol]

    @property
    def step_to_execute(self):
        return self._step_to_execute

    @property
    def sample_list(self):
        return self._sample_list

    @property
    def sample_paths(self):
        return self._sample_paths

    @property
    def jobs(self):
        jobs = []
        for step in self.step_to_execute:
            jobs.extend(step.jobs)
        return jobs

    @property
    def json(self):
        return self._json

    @classmethod
    def genpipes_version(cls):
        if cls.version is None:
            cls.version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                             "VERSION"), 'r').read().split('\n')[0]
        return cls.version

    # Given a list of lists of input files, return the first valid list of input files which can be found either in previous jobs output files or on file system.
    # Thus, a job with several candidate lists of input files can find out the first valid one.
    def select_input_files(self, candidate_input_files):
        log.debug("candidate_input_files: \n" + str(candidate_input_files))

        selected_input_files = []

        # Create a reversed copy to pop the candidates ordered by priority
        remaining_candidate_input_files = list(candidate_input_files)
        remaining_candidate_input_files.reverse()
        previous_jobs_output_files = set([output_file for job in self.jobs for output_file in job.output_files])

        while not selected_input_files and remaining_candidate_input_files:
            input_files = [_f for _f in remaining_candidate_input_files.pop() if _f]
            # Skip empty candidate input files
            if input_files:
                # dependency_jobs() checks if the current candidate input files is valid, otherwise raises an exception
                try:
                    job = Job(input_files=input_files)
                    job.output_dir = self.output_dir
                    self.dependency_jobs(job)
                    selected_input_files = input_files
                except Exception as e:
                    log.debug("Caught Exception for candidate input file: " +  ", ".join(input_files))
                    log.debug(e)

        if selected_input_files:
            log.debug("selected_input_files: " + ", ".join(input_files) + "\n")
            return selected_input_files
        else:
            _raise(SanitycheckError("Error: missing candidate input files: " + str(candidate_input_files) +
                " neither found in dependencies nor on file system!"))


    def dependency_jobs(self, current_job):
        dependency_jobs = []
        dependency_input_files = set()
        for step in self.step_to_execute:
            for step_job in step.jobs:
                # If current job input files intersect with step job output files, step job is a dependency
                shared_files = set(current_job.input_files).intersection(set(step_job.output_files))
                if shared_files:
                    dependency_jobs.append(step_job)
                    dependency_input_files.update(shared_files)

        # Check if job input files not found in dependencies are on file system
        missing_input_files = set()
        # Add current_job.output_files in case of "... && ..." command
        # where first command output becomes second command input
        for remaining_input_file in set(current_job.input_files).difference(dependency_input_files).difference(set(current_job.output_files)):
            # Use 'exists' instead of 'isfile' since input file can be a directory
            if not os.path.exists(current_job.abspath(remaining_input_file)):
                missing_input_files.add(remaining_input_file)
        if missing_input_files:
            raise Exception("Warning: missing input files for job " + current_job.name + ": " +
                ", ".join(missing_input_files) + " neither found in dependencies nor on file system!")

        return dependency_jobs

    def create_jobs(self):
        for step in self.step_to_execute:
            if self.sanity_check :
                log.warning("* Checking jobs for step " + step.name + "...")
            else :
                log.info("Create jobs for step " + step.name + "...")
            jobs = step.create_jobs()
            for job in jobs:
                # Job name is mandatory to create job .done file name
                if not job.name:
                    _raise(SanitycheckError("Error: job \"" + job.command + "\" has no name!"))

                log.debug("Job name: " + job.name)
                log.debug("Job input files:\n  " + "\n  ".join(job.input_files))
                log.debug("Job output files:\n  " + "\n  ".join(job.output_files) + "\n")

                # Job .done file name contains the command checksum.
                # Thus, if the command is modified, the job is not up-to-date anymore.
                job.done = os.path.join("job_output", step.name, job.name + "." + hashlib.md5(job.command_with_modules.encode('utf-8')).hexdigest() + ".mugqic.done")
                job.output_dir = self.output_dir
                job.dependency_jobs = self.dependency_jobs(job)

                if not self.force_jobs and job.is_up2date():
                    log.info("Job " + job.name + " up to date... skipping\n")
                else:
                    step.add_job(job)
                    if job.samples:
                        for sample in job.samples:
                            if sample not in self.sample_list:
                                self.sample_list.append(sample)

            if not self.sanity_check :
                log.info("Step " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "") + " created" + ("" if step.jobs else "... skipping") + "\n")

        # Now create the json dumps for all the samples if not already done
        if self.json:
            # Check if portal_output_dir is set from a valid environment variable
            self.portal_output_dir = global_config_parser.param('DEFAULT', 'portal_output_dir', required=False)
            log.info(self.portal_output_dir.startswith("$"))
            log.info(os.environ.get(re.search("^\$(.*)\/?", self.portal_output_dir).group(1)))
            if self.portal_output_dir.startswith("$") and (os.environ.get(re.search("^\$(.*)\/?", self.portal_output_dir).group(1)) is None or os.environ.get(re.search("^\$(.*)\/?", self.portal_output_dir).group(1)) == ""):
                if self.portal_output_dir == "$PORTAL_OUTPUT_DIR":
                    self.portal_output_dir = ""
                    log.info(" --> PORTAL_OUTPUT_DIR environment variable is not set... no JSON file will be generated during analysis...\n")
                    self._json = False
                else:
                    _raise(SanitycheckError("Environment variable \"" + re.search("^\$(.*)\/?", self.portal_output_dir).group(1) + "\" does not exist or is not valid!"))
            elif not os.path.isdir(os.path.expandvars(self.portal_output_dir)):
                _raise(SanitycheckError("Directory path \"" + self.portal_output_dir + "\" does not exist or is not a valid directory!"))
            else:
                for sample in self.sample_list:
                    self.sample_paths.append(jsonator.create(self, sample))

        log.info("TOTAL: " + str(len(self.jobs)) + " job" + ("s" if len(self.jobs) > 1 else "") + " created" + ("" if self.jobs else "... skipping") + "\n")

    def submit_jobs(self):
        self.job_scheduler.submit(self)

    def report_jobs(self, output_dir=None):
        if not output_dir:
            output_dir = self.output_dir  # Default to pipeline output directory
        report_files = []
        for job in self.jobs:
            # Retrieve absolute paths of report files
            for report_file in job.report_files:
                if report_file not in report_files:
                    if os.path.exists(os.path.join(output_dir, report_file)) :
                        report_files.append(report_file)
                    else:
                        log.warning("Report file: " + report_file + " not found!... skipping")
        if report_files:
            # Copy images and other HTML dependencies into report directory
            # Print pandoc command with all markdown report files and config/references sections at the end
            print("""\
#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

module load {module_pandoc}
cd {output_dir}
mkdir -p report
cp -r \\
  {self.report_template_dir}/css/ \\
  {self.report_template_dir}/images/ \\
  {self.report_template_dir}/_js/ \\
  {self.report_template_dir}/config_and_references.md \\
  report/
cp {config_file} report/config.ini
pandoc \\
  --toc \\
  --toc-depth=6 \\
  --template={self.report_template_dir}/template.html \\
  --css=css/style.css \\
  --metadata title="{title}" \\
  --citeproc \\
  <(pandoc --to=markdown \\
    --template {introduction} \\
    --variable date="{date}" \\
    {introduction} \\
  ) \\
  {report_files} \\
  report/config_and_references.md \\
  --output report/index.html""".format(
                module_pandoc=global_config_parser.param('report', 'module_pandoc'),
                output_dir=output_dir,
                config_file=global_config_parser.filepath,
                self=self,
                title=global_config_parser.param('report', 'title'),
                introduction=os.path.join(self.report_template_dir, self.__class__.__name__ + '.introduction.md'),
                date=datetime.datetime.now().strftime("%Y-%m-%d"),
                report_files=" \\\n  ".join(report_files)
            ))

    def clean_jobs(self):
        abspath_removable_files = []
        for job in self.jobs:
            # Retrieve absolute paths of removable files
            abspath_removable_files.extend([job.abspath(removable_file) for removable_file in job.removable_files])
        # Remove removable file duplicates but keep the order
        for removable_file in list(collections.OrderedDict.fromkeys(abspath_removable_files)):
            if os.path.exists(removable_file):
                print("rm -rf " + removable_file)

# Return a range list given a string.
# e.g. parse_range('1,3,5-12') returns [1, 3, 5, 6, 7, 8, 9, 10, 11, 12]
def parse_range(astr):
    result = set()
    for part in astr.split(','):
        x = part.split('-')
        result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)

class ValidateContainer(argparse.Action):

    VALID_TYPE = ('singularity', 'wrapper')

    def __call__(self, parser, args, values, option_string=None):
        c_type, container = values
        if c_type not in self.VALID_TYPE:
            raise ValueError('{} is not supported, choose from {}'.format(c_type, self.VALID_TYPE))
        Container = collections.namedtuple('container', 'type name')

        setattr(args, self.dest, Container(c_type, container))
