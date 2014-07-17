#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import hashlib
import logging
import os
import re

# MUGQIC Modules
from config import *
from scheduler import *
from step import *

log = logging.getLogger(__name__)

class Pipeline:
    def __init__(self, args):
        self._output_dir = os.path.abspath(args.output_dir)
        self._scheduler = create_scheduler(args.job_scheduler)
        self._force_jobs = args.force

        step_counter = collections.Counter(self.steps)
        duplicated_steps = [step.__name__ for step in step_counter if step_counter[step] > 1]
        if duplicated_steps:
            raise Exception("Error: pipeline contains duplicated steps: " + ", ".join(duplicated_steps) + "!")
        else:
            self.step_list = [Step(step) for step in self.steps]

        if re.search("^\d+([,-]\d+)*$", args.steps):
            self._step_range = [self.step_list[i - 1] for i in parse_range(args.steps)]
        else:
            raise Exception("Error: step range \"" + args.steps +
                "\" is invalid (should match \d+([,-]\d+)*)!")

        self.create_jobs()

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
            if not os.path.isfile(current_job.abspath(remaining_input_file)):
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

                # Job .done file name contains the command checksum.
                # Thus, if the command is modified, the job is not up-to-date anymore.
                job.done = os.path.join("job_output", step.name, job.name + "." + hashlib.md5(job.command).hexdigest() + ".mugqic.done")
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

# Default command line argument parser holding common options shared by all pipelines
class PipelineArgumentParser(argparse.ArgumentParser):

    def __init__(self, steps):
        # Create ArgumentParser with numbered step list as epilog
        argparse.ArgumentParser.__init__(self, formatter_class=argparse.RawDescriptionHelpFormatter, epilog="Steps:\n" + "\n".join([str(idx + 1) + "- " + step.__name__ for idx, step in enumerate(steps)]))

        # Common options for all pipelines
        self.add_argument("-c", "--config", help="config INI-style file", nargs="+", type=file, required=True)
        self.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'", required=True)
        self.add_argument("-o", "--output-dir", help="output directory (default: current)", default=os.getcwd())
        self.add_argument("-j", "--job-scheduler", help="job scheduler type (default: torque)", choices=["torque", "batch", "daemon"], default="torque")
        self.add_argument("-f", "--force", help="force creation of jobs even if up to date (default: false)", action="store_true")
        self.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")

    def parse_args(self):
        args = argparse.ArgumentParser.parse_args(self)
        logging.basicConfig(level=getattr(logging, args.log.upper()))
        config.parse_files(args.config)
        return args
