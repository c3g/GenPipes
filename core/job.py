#!/usr/bin/env python

# Python Standard Modules
import collections
import logging
import os

# MUGQIC Modules
from config import *

log = logging.getLogger(__name__)

class Job:

    def __init__(self, input_files=[], output_files=[], module_entries = [], name="", command=""):
        # Remove undefined input/output files if any
        self._input_files = filter(None, input_files)
        self._output_files = filter(None, output_files)

        # Retrieve modules from config, removing duplicates but keeping the order
        self._modules = list(collections.OrderedDict.fromkeys([config.param(section, option) for section, option in module_entries]))

        self._name = name
        self._command = command

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
    def output_dir(self):
        return self._output_dir

    @property
    def input_files(self):
        return self._input_files

    @property
    def output_files(self):
        return self._output_files

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

    def abspath(self, file):
        tmp_file = os.path.expandvars(file)
        if not os.path.isabs(tmp_file):
            # File path is relative to the job output directory
            tmp_file = os.path.normpath(os.path.join(self.output_dir, tmp_file))
        return tmp_file

    def is_up2date(self):
        # If job has dependencies, job is not up to date
        if self.dependency_jobs:
            return False

        # Retrieve absolute paths for .done, input and output files to avoid redundant OS function calls
        abspath_done = self.abspath(self.done)
        abspath_input_files = [self.abspath(input_file) for input_file in self.input_files]
        abspath_output_files = [self.abspath(output_file) for output_file in self.output_files]

        # If any .done, input or output file is missing, job is not up to date
        for file in [abspath_done] + abspath_input_files + abspath_output_files:
            if not os.path.isfile(file):
                return False

        # Retrieve latest input file modification time i.e. maximum stat mtime
        latest_input_time = max([os.stat(input_file).st_mtime for input_file in abspath_input_files])

        # Same with earliest output file modification time
        earliest_output_time = min([os.stat(output_file).st_mtime for output_file in abspath_output_files])

        # If any input file is strictly more recent than all output files, job is not up to date
        # Use strictly '>' otherwise jobs like "cmd1 in1 > out1 && cmd2 out1 > out2" would always be out of date
        # (out1 is both output and input of the command)
        if latest_input_time > earliest_output_time:
            return False

        # If all previous tests passed, job is up to date
        return True


# Create a new job by concatenating a list of jobs together
def concat_jobs(jobs):

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

    job = Job(input_files, output_files)
    job.modules = modules

    # Merge commands
    job.command = " && \\\n".join([job_item.command for job_item in jobs])

    return job

# Create a new job by piping a list of jobs together
def pipe_jobs(jobs):

    job = Job(jobs[0].input_files, jobs[-1].output_files)

   # Merge all modules
    modules = []
    for job_item in jobs:
        modules.extend(job_item.modules)

    # Remove duplicates if any, keeping the order
    modules = list(collections.OrderedDict.fromkeys([module for module in modules]))
    job.modules = modules

    # Merge commands
    job.command = " | \\\n".join([job_item.command for job_item in jobs])

    return job
