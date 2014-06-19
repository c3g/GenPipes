#!/usr/bin/env python

# Python Standard Modules
import collections
import logging
import os

# MUGQIC Modules
from config import *

log = logging.getLogger(__name__)

class Job:

    def __init__(self, input_files, output_files, module_entries = []):
        # Remove undefined input/output files if any
        self._input_files = filter(None, input_files)
        self._output_files = filter(None, output_files)

        # Retrieve modules from config, removing duplicates but keeping the order
        self._modules = list(collections.OrderedDict.fromkeys([config.param(section, option) for section, option in module_entries]))

    def show(self):
        print("Job: input_files: " + \
            ", ".join(self.input_files))

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def input_files(self):
        return self._input_files

    @property
    def output_files(self):
        return self._output_files

    # Input/output files are first defined relative to the pipeline output directory
    def update_files(self, output_dir):
        self.input_files = [file if os.path.isabs(file) else os.path.join(output_dir, file) for file in self.input_files]
        self.output_files = [file if os.path.isabs(file) else os.path.join(output_dir, file) for file in self.output_files]

    @property
    def done(self):
        return self._done

    @property
    def dependency_jobs(self):
        return self._dependency_jobs

    @property
    def modules(self):
        return self._modules

    @property
    def command(self):
        return self._command

    @property
    def command_with_modules(self):
        command = self.command
        if self.modules:
            command = "module load " + " ".join(self.modules) + " && \\\n" + command
        return command

    def is_up2date(self):
        # Job is up to date by default
        is_job_up2date = True

        # If .done is missing, job is not up to date
        if not os.path.isfile(os.path.expandvars(self.done)):
            is_job_up2date = False

        # If any input file is missing, job is not up to date
        for input_file in self.input_files:
            if not os.path.isfile(os.path.expandvars(input_file)):
                is_job_up2date = False

        # If any output file file is missing, job is not up to date
        for output_file in self.output_files:
            if not os.path.isfile(os.path.expandvars(output_file)):
                is_job_up2date = False

        # Skip further tests if job is already out of date
        if is_job_up2date:
            # Retrieve latest input file modification time i.e. maximum stat mtime
            # Use 'echo' system command to expand environment variables in input file paths if any
            latest_input_time = max([os.stat(os.path.expandvars(input_file)).st_mtime for input_file in self.input_files])

            # Same with earliest output file modification time
            earliest_output_time = max([os.stat(os.path.expandvars(output_file)).st_mtime for output_file in self.output_files])
            is_job_up2date = earliest_output_time > latest_input_time

        return is_job_up2date


# Create a new job from a job list by merging their modules and commands with a specified separator
def group_jobs(jobs, separator):

    # At least 2 jobs are required
    if len(jobs) < 2:
        raise Exception("Error: group_jobs requires at least 2 jobs!")

    # Merge all input/output files and modules
    input_files = []
    output_files = []
    modules = []
    for job_item in jobs:
        input_files.extend(job_item.input_files)
        output_files.extend(job_item.output_files)
        modules.extend(job_item.modules)

    # Remove duplicates if any, keeping the order
    input_files = list(collections.OrderedDict.fromkeys([input_file for input_file in input_files]))
    output_files = list(collections.OrderedDict.fromkeys([output_file for output_file in output_files]))
    modules = list(collections.OrderedDict.fromkeys([module for module in modules]))

    job = Job(input_files, output_files);
    job.modules = modules

    # Merge commands with specified separator
    job.command = separator.join([job_item.command for job_item in jobs])

    return job;

# Create a new job by piping a list of jobs together
def pipe_jobs(jobs):
    return group_jobs(jobs, " | \\\n")

# Create a new job by concatenating a list of jobs together
def concat_jobs(jobs):
    return group_jobs(jobs, " && \\\n")
