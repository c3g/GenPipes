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
if sys.version_info < (2,7):
    raise SystemExit("Incompatible Python version: " + sys.version + "\nPython 2.7 or higher is required")

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
from .config import config, _raise, SanitycheckError
from .job import Job
from .scheduler import create_scheduler
from .step import Step

from bfx import jsonator

log = logging.getLogger(__name__)

class Pipeline(object):

    def __init__(self):
        self._timestamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        self._args = self.argparser.parse_args()
        self._genpipes_version = subprocess.check_output("cat " + os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))), "VERSION"), shell=True)
            
        if self.protocol is None:
            step_list = self.steps
        elif self.args.help:
            step_list = []
            for i in range(0, len(self.protocol)):
                step_list = step_list + self.steps[i]
            tmp_set_list = {}
            step_list = [tmp_set_list.setdefault(step, step) for step in step_list if step not in tmp_set_list]
        else:
            pos = 0
            for i in range(0, len(self.protocol)):
                if self.protocol[i] == self.args.type:
                    pos = i
            step_list = self.steps[pos]

        if self.args.sanity_check:
            logging.basicConfig(stream=sys.stdout, level=logging.WARNING, format='%(message)s')
            log.warning("Sanity Check Report :")
        else:
            logging.basicConfig(level=getattr(logging, self.args.log.upper()))

        if self.args.help:
            print(textwrap.dedent("""\
              [TOC]

              {pipeline_doc}

              Usage
              -----
              ```
              #!text

              {help}
              ```

              {step_doc}
            """).format(
            pipeline_doc=textwrap.dedent(self.__doc__ or ""),
            help=self.argparser.format_help(),
            overview=self.__doc__ or "",
            #step_doc="\n".join([str(idx + 1) + "- " + step.__name__ + "\n" + "-" * len(str(idx + 1) + "- " + step.__name__) + (textwrap.dedent(step.__doc__) if step.__doc__ else "") for idx, step in enumerate(step_list)])
            step_doc="\n".join([step.__name__ + "\n" + "-" * len(step.__name__) + (textwrap.dedent(step.__doc__) if step.__doc__ else "") for step in step_list])
            ))
            self.argparser.exit()

        # Normal pipeline execution
        if self.args.config:
            if self.args.sanity_check: config.sanity = True
            config.parse_files(self.args.config)
        else:
            self.argparser.error("argument -c/--config is required!")

        # Create a config trace from merged config file values
        with open(self.__class__.__name__ + ".config.trace.ini", 'w') as config_trace:
            config_trace.write(textwrap.dedent("""\
              # {self.__class__.__name__} Config Trace
              # Created on: {self.timestamp}
              # From:
              #   {config_files}
              # DO NOT EDIT THIS AUTOMATICALLY GENERATED FILE - edit the master config files

            """).format(config_files="\n#   ".join([config_file.name for config_file in self.args.config]), self=self))
            config.write(config_trace)
            config._filepath = os.path.abspath(config_trace.name)

        self._output_dir = os.path.abspath(self.args.output_dir)
        self._scheduler = create_scheduler(self.args.job_scheduler, self.args.config, container=self.args.container,
                                           genpipes_file=self.args.genpipes_file)

        self._json = True
        if self.args.no_json:
            self._json = False

        step_counter = collections.Counter(step_list)
        duplicated_steps = [step.__name__ for step in step_counter if step_counter[step] > 1]
        if duplicated_steps:
            raise Exception("Error: pipeline contains duplicated steps: " + ", ".join(duplicated_steps) + "!")
        else:
            self._step_list = [Step(step) for step in step_list]

        if self.args.steps:
            if re.search("^\d+([,-]\d+)*$", self.args.steps):
                self._step_range = [self.step_list[i - 1] for i in parse_range(self.args.steps)]
            else:
                raise Exception("Error: step range \"" + self.args.steps +
                    "\" is invalid (should match \d+([,-]\d+)*)!")
        else:
#            self.argparser.error("argument -s/--steps is required!")
            log.warning("No step provided by the user => launching the entire pipeline\n")
            self._step_range = self.step_list
                

        self._sample_list = []
        self._sample_paths = []

        # For job reporting, all jobs must be created first, no matter whether they are up to date or not
        if self.args.report:
            self._force_jobs = True
            self.create_jobs()
            self.report_jobs()
        # For job cleaning, all jobs must be created first, no matter whether they are up to date or not
        elif self.args.clean:
            self._force_jobs = True
            self.create_jobs()
            self.clean_jobs()
        elif self.args.sanity_check:
            config.sanity = True
            self._force_jobs = self.args.force
            self.create_jobs()
        else:
            try:
                self._force_jobs = self.args.force
                self.create_jobs()
                self.submit_jobs()
            except SanitycheckError as e:
                log.error("""
***The pipeline encountered an error :
    {error}
***Please try running the pipeline in SANITY CHECK mode using the '--sanity-check' flag to check for more potential issues...""".format(
                   error=e
                   ))
                exit(1)
            
    # Pipeline command line arguments parser
    @property
    def argparser(self):
        if self.protocol is None:
            steps = "\n".join([str(idx + 1) + "- " + step.__name__ for idx, step in enumerate(self.steps)])
        else:
            steps = ""
            for i in range(0, len(self.protocol)):
                steps += "\n----\n"+self.protocol[i]+":\n"+"\n".join([str(idx + 1) + "- " + step.__name__ for idx, step in enumerate(self.steps[i])])

        if not hasattr(self, "_argparser"):
            epilog = textwrap.dedent("""\
                Steps:
                ------
                {steps}
            """).format(steps=steps)

            # Create ArgumentParser with numbered step list and description as epilog
            self._argparser = argparse.ArgumentParser(
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog,
                    conflict_handler='resolve')

            # Common options for all pipelines
            self._argparser.add_argument("--help", help="show detailed description of pipeline and steps",
                                         action="store_true")
            self._argparser.add_argument("-c", "--config", help="config INI-style list of files; config parameters are "
                                                                "overwritten based on files order", nargs="+", type=argparse.FileType('r'))
            self._argparser.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'")
            self._argparser.add_argument("-o", "--output-dir", help="output directory (default: current)",
                                         default=os.getcwd())
            self._argparser.add_argument("-j", "--job-scheduler", help="job scheduler type (default: slurm)",
                                         choices=["pbs", "batch", "daemon", "slurm"], default="slurm")
            self._argparser.add_argument("-f", "--force", help="force creation of jobs even if up to date "
                                                               "(default: false)", action="store_true")
            self._argparser.add_argument("--no-json", help="do not create JSON file per analysed sample to track the "
                                                           "analysis status (default: false i.e. JSON file will be "
                                                           "created)", action="store_true")
            self._argparser.add_argument("--report", help="create 'pandoc' command to merge all job markdown report "
                                                          "files in the given step range into HTML, if they exist; if "
                                                          "--report is set, --job-scheduler, --force, --clean options "
                                                          "and job up-to-date status are ignored (default: false)",
                                         action="store_true")
            self._argparser.add_argument("--clean", help="create 'rm' commands for all job removable files in the given"
                                                         " step range, if they exist; if --clean is set,"
                                                         " --job-scheduler, --force options and job up-to-date status "
                                                         "are ignored (default: false)", action="store_true")
            self._argparser.add_argument("-l", "--log", help="log level (default: info)",
                                         choices=["debug", "info", "warning", "error", "critical"], default="info")
            self._argparser.add_argument("--sanity-check", help="run the pipeline in `sanity check mode` to verify that"
                                                                " all the input files needed for the pipeline to run "
                                                                "are available on the system (default: false)",
                                         action="store_true")
            self._argparser.add_argument("--container", nargs=2,
                                         help="Run inside a container providing a valid "
                                         "singularity image path", action=ValidateContainer,
                                          metavar=("{wrapper, singularity}",
                                                   "<IMAGE PATH>"))
            self._argparser.add_argument("--genpipes_file", '-g', default=sys.stdout, type=argparse.FileType('w'),
                                         help="Command file output path. This is the command used to process "
                                              "the data, or said otherwise, this command will \"run the "
                                              "Genpipes pipeline\". Will be redirected to stdout if the "
                                              "option is not provided.")

        return self._argparser

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
        return os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))), "bfx", "report")

    @property
    def scheduler(self):
        return self._scheduler

    @property
    def force_jobs(self):
        return self._force_jobs

    @property
    def protocol(self):
        return self._protocol

    @property
    def steps(self):
        # Needs to be defined in pipeline child class
        raise NotImplementedError

    @property
    def step_list(self):
        return self._step_list

    @property
    def step_range(self):
        return self._step_range

    @property
    def sample_list(self):
        return self._sample_list

    @property
    def sample_paths(self):
        return self._sample_paths

    @property
    def jobs(self):
        jobs = []
        for step in self.step_range:
            jobs.extend(step.jobs)
        return jobs

    @property
    def json(self):
        return self._json

    @property
    def genpipes_version(self):
        return self._genpipes_version

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
                    log.debug("Caught Exception for candidate input file: " + ", ".join(input_files))
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
        for step in self.step_range:
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
        for step in self.step_range:
            if self.args.sanity_check:
                log.warn("* Checking jobs for step " + step.name + "...")
            else:
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

            if not self.args.sanity_check:
                log.info("Step " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "") + " created" + ("" if step.jobs else "... skipping") + "\n")

        # Now create the json dumps for all the samples if not already done
        if self.json:
            # Check if portal_output_dir is set from a valid environment variable
            self.portal_output_dir = config.param('DEFAULT', 'portal_output_dir', required=False)
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
        self.scheduler.submit(self)

        ## Print a copy of sample JSONs for the genpipes dashboard
#        if self.json and self.portal_output_dir != "":
#            copy_commands = []
#            for i, sample in enumerate(self.sample_list):
#                input_file = self.sample_paths[i]
#                output_file = os.path.join(self.portal_output_dir, '$USER.' + sample.name + '.' + uuid4().get_hex() + '.json')
#                copy_commands.append("cp \"{input_file}\" \"{output_file}\"".format(
#                    input_file=input_file, output_file=output_file))
#            print(textwrap.dedent("""
#                #------------------------------------------------------------------------------
#                # Print a copy of sample JSONs for the genpipes dashboard
#                #------------------------------------------------------------------------------
#                {copy_commands}
#            """).format(copy_commands='\n'.join(copy_commands)))

    def report_jobs(self, output_dir=None):
        if not output_dir:
            output_dir = self.output_dir  # Default to pipeline output directory
        report_files = []
        for job in self.jobs:
            # Retrieve absolute paths of report files
            for report_file in job.report_files:
                if report_file not in report_files:
                    if os.path.exists(os.path.join(output_dir, report_file)):
                        report_files.append(report_file)
                    else:
                        log.warn("Report file: " + report_file + " not found!... skipping")
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
                module_pandoc=config.param('report', 'module_pandoc'),
                output_dir=output_dir,
                config_file=config.filepath,
                self=self,
                title=config.param('report', 'title'),
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
