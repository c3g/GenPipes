#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import datetime
import hashlib
import logging
import os
import re

# MUGQIC Modules
from config import *
from scheduler import *
from step import *

log = logging.getLogger(__name__)

class Pipeline(object):
    def __init__(self):
        self._timestamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        self._args = self.argparser.parse_args()

        logging.basicConfig(level=getattr(logging, self.args.log.upper()))
        config.parse_files(self.args.config)

        # Create a config trace from merged config file values
        with open(self.__class__.__name__ + ".config.trace.ini", 'wb') as config_trace:        
            config_trace.write("""\
# {self.__class__.__name__} Config Trace
# Created on: {self.timestamp}
# From:
#   {config_files}
# DO NOT EDIT THIS AUTOMATICALLY GENERATED FILE - edit the master config files

""".format(config_files="\n#   ".join([config_file.name for config_file in self.args.config]), self=self))
            config.write(config_trace)
            config.filepath = os.path.abspath(config_trace.name)

        self._output_dir = os.path.abspath(self.args.output_dir)
        self._scheduler = create_scheduler(self.args.job_scheduler)
        self._force_jobs = self.args.force

        step_counter = collections.Counter(self.steps)
        duplicated_steps = [step.__name__ for step in step_counter if step_counter[step] > 1]
        if duplicated_steps:
            raise Exception("Error: pipeline contains duplicated steps: " + ", ".join(duplicated_steps) + "!")
        else:
            self._step_list = [Step(step) for step in self.steps]

        if re.search("^\d+([,-]\d+)*$", self.args.steps):
            self._step_range = [self.step_list[i - 1] for i in parse_range(self.args.steps)]
        else:
            raise Exception("Error: step range \"" + self.args.steps +
                "\" is invalid (should match \d+([,-]\d+)*)!")

        self.create_jobs()

    # Pipeline command line arguments parser
    @property
    def argparser(self):
        if not hasattr(self, "_argparser"):
            # Create ArgumentParser with numbered step list as epilog
            self._argparser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=self.name, epilog="Steps:\n" + "\n".join([str(idx + 1) + "- " + step.__name__ for idx, step in enumerate(self.steps)]))

            # Common options for all pipelines
            self._argparser.add_argument("-v", "--version", action="version", version=self.name, help="show the version information and exit")
            self._argparser.add_argument("-c", "--config", help="config INI-style list of files; config parameters are overwritten based on files order", nargs="+", type=file, required=True)
            self._argparser.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'", required=True)
            self._argparser.add_argument("-o", "--output-dir", help="output directory (default: current)", default=os.getcwd())
            self._argparser.add_argument("-j", "--job-scheduler", help="job scheduler type (default: pbs)", choices=["pbs", "batch"], default="pbs")
            self._argparser.add_argument("-f", "--force", help="force creation of jobs even if up to date (default: false)", action="store_true")
            self._argparser.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")

        return self._argparser

    @property
    def name(self):
        return self.__class__.__name__ + "-" + open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "VERSION"), 'r').read().split('\n')[0]

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
    def scheduler(self):
        return self._scheduler

    @property
    def force_jobs(self):
        return self._force_jobs

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
    def jobs(self):
        jobs = []
        for step in self.step_range:
            jobs.extend(step.jobs)
        return jobs

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
            raise Exception("Error: missing input files for job " + current_job.name + ": " +
                ", ".join(missing_input_files) + " neither found in dependencies nor on file system!")

        return dependency_jobs

    def create_jobs(self):
        for step in self.step_range:
            log.info("Create jobs for step " + step.name + "...")
            jobs = step.create_jobs()
            for job in jobs:
                # Job name is mandatory to create job .done file name
                if not job.name:
                    raise Exception("Error: job \"" + job.command + "\" has no name!")

                log.debug("Job name: " + job.name)
                log.debug("Job input files:\n  " + "\n  ".join(job.input_files))
                log.debug("Job output files:\n  " + "\n  ".join(job.output_files) + "\n")

                # Job .done file name contains the command checksum.
                # Thus, if the command is modified, the job is not up-to-date anymore.
                job.done = os.path.join("job_output", step.name, job.name + "." + hashlib.md5(job.command_with_modules).hexdigest() + ".mugqic.done")
                job.output_dir = self.output_dir
                job.dependency_jobs = self.dependency_jobs(job)
                if not self.force_jobs and job.is_up2date():
                    log.info("Job " + job.name + " up to date... skipping")
                else:
                    step.add_job(job)
            log.info("Step " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "") + " created" + ("" if step.jobs else "... skipping") + "\n")
        log.info("TOTAL: " + str(len(self.jobs)) + " job" + ("s" if len(self.jobs) > 1 else "") + " created" + ("" if self.jobs else "... skipping") + "\n")

    def submit_jobs(self):
        self.scheduler.submit(self)


# Return a range list given a string.
# e.g. parse_range('1,3,5-12') returns [1, 3, 5, 6, 7, 8, 9, 10, 11, 12]
def parse_range(astr):
    result = set()
    for part in astr.split(','):
        x = part.split('-')
        result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)
