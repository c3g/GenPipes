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

# Python Standard Modules
from collections import OrderedDict
import datetime
import logging
import os
import sys

# MUGQIC Modules
from .config import global_config_parser

log = logging.getLogger(__name__)

class Job(object):


    def __init__(self, input_files=[], output_files=[], module_entries=[],
                 name="", command="", report_files=[], multiqc_files=[],
                 removable_files=[], samples=[]):
        # Remove undefined input/output/removable files if any
        self._input_files = [_f for _f in input_files if _f]
        self._output_files = [_f for _f in output_files if _f]
        self._report_files = [_f for _f in report_files if _f]
        self._multiqc_files = [_f for _f in multiqc_files if _f]
        self._removable_files = [_f for _f in removable_files if _f]

        # Retrieve modules from config, removing duplicates but keeping the order
        self._modules = list(OrderedDict.fromkeys([global_config_parser.param(section, option)
                                                               for section, option in module_entries]))

        self._name = name
        self._command = command
        self._samples = samples

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def output_dir(self):
        return self._output_dir

    @output_dir.setter
    def output_dir(self, value):
        self._output_dir = value

    @property
    def input_files(self):
        return self._input_files

    @input_files.setter
    def input_files(self, value):
        self._input_files = value

    @property
    def output_files(self):
        return self._output_files

    @output_files.setter
    def output_files(self, value):
        self._output_files = value

    @property
    def report_files(self):
        return self._report_files

    @report_files.setter
    def report_files(self, value):
        self._report_files = value

    @property
    def multiqc_files(self):
        return self._multiqc_files

    @multiqc_files.setter
    def multiqc_files(self, value):
        self._multiqc_files = value

    @property
    def removable_files(self):
        return self._removable_files

    @removable_files.setter
    def removable_files(self, value):
        self._removable_files = value

    @property
    def done(self):
        return self._done

    @done.setter
    def done(self, value):
        self._done = value

    @property
    def dependency_jobs(self):
        return self._dependency_jobs

    @dependency_jobs.setter
    def dependency_jobs(self, value):
        self._dependency_jobs = value

    @property
    def modules(self):
        return self._modules

    @modules.setter
    def modules(self, value):
        self._modules = value

    @property
    def command(self):
        return self._command

    @command.setter
    def command(self, value):
        self._command = value

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, value):
        self._samples = value

    @property
    def command_with_modules(self):
        command = self.command
        if self.modules:
            command = "module purge && \\\n" + "module load " + " ".join(self.modules) + " && \\\n" + command
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
            log.debug("Job " + self.name + " NOT up to date")
            log.debug("Dependency jobs:\n  " + "\n  ".join([job.name for job in self.dependency_jobs]) + "\n")
            return False

        # Retrieve absolute paths for .done, input and output files to avoid redundant OS function calls
        abspath_done = self.abspath(self.done)
        abspath_input_files = [self.abspath(input_file) for input_file in self.input_files]
        abspath_output_files = [self.abspath(output_file) for output_file in self.output_files]

        # If any .done, input or output file is missing, job is not up to date
        for file in [abspath_done] + abspath_input_files + abspath_output_files:
            # Use 'exists' instead of 'isfile' since input/output files can be directories
            if not os.path.exists(file):
                log.debug("Job " + self.name + " NOT up to date")
                log.debug("Input, output or .done file missing: " + file)
                return False

        # Retrieve latest input file by modification time i.e. maximum stat mtime
        # Use lstat to avoid following symbolic links
        latest_input_file = max(abspath_input_files, key=lambda input_file: os.lstat(input_file).st_mtime)
        latest_input_time = os.lstat(latest_input_file).st_mtime

        # Same with earliest output file by modification time
        earliest_output_file = min(abspath_output_files, key=lambda output_file:os.lstat(output_file).st_mtime)
        earliest_output_time = os.lstat(earliest_output_file).st_mtime

        # If any input file is strictly more recent than all output files, job is not up to date
        if latest_input_time > earliest_output_time:
            log.debug("Job " + self.name + " NOT up to date")
            log.debug("Latest input file modification time: " + latest_input_file + " " + datetime.datetime.fromtimestamp(latest_input_time).isoformat() + " > earliest output file modification time: " + earliest_output_file + " " + datetime.datetime.fromtimestamp(earliest_output_time).isoformat() + "\n")
            return False

        # If all previous tests passed, job is up to date
        return True

# Create a new job by concatenating a list of jobs together
def concat_jobs(jobs, name="", input_dependency=[], output_dependency=[], samples=[], removable_files=[], report_files=[]):

    # Merge all input/output/report/removable files and modules
    input_files = []
    output_files = []
    multiqc_files = []
    modules = []
    for job_item in [job_item for job_item in jobs if job_item]:
        input_files.extend([input_file for input_file in job_item.input_files if input_file not in input_files and input_file not in output_files])
        output_files.extend([output_file for output_file in job_item.output_files if output_file not in output_files])
        report_files.extend([report_file for report_file in job_item.report_files if report_file not in report_files])
        multiqc_files.extend([multiqc_file for multiqc_file in job_item.multiqc_files if multiqc_file not in multiqc_files])
        removable_files.extend([removable_file for removable_file in job_item.removable_files if removable_file not in removable_files])
        modules.extend([module for module in job_item.modules if module not in modules])
        samples.extend([sample for sample in job_item.samples if sample not in samples])

    if input_dependency:
        input_files=input_dependency
    if output_dependency:
        output_files=output_dependency

    job = Job(input_files, output_files, name=name, report_files=report_files, multiqc_files=multiqc_files, removable_files=removable_files, samples=samples)
    job.modules = modules

    # Merge commands
    job.command = " && \\\n".join([job_item.command for job_item in jobs if job_item and job_item.command])

    return job

# Create a new job by piping a list of jobs together
def pipe_jobs(jobs, name="", input_dependency=[], output_dependency=[], samples=[], removable_files=[], report_files=[]):

    # Remove 'None' jobs from the job list, may happened in the case of conditional jobs
    jobs = [job for job in jobs if job]
    job = Job(jobs[0].input_files, jobs[-1].output_files, name=name)

    # Merge all report/removable files and modules
    multiqc_files = []
    modules = []
    for job_item in [job_item for job_item in jobs if job_item]:
        report_files.extend(job_item.report_files)
        multiqc_files.extend(job_item.multiqc_files)
        removable_files.extend(job_item.removable_files)
        modules.extend(job_item.modules)
        samples.extend([sample for sample in job_item.samples if sample not in samples])

    # Remove duplicates if any, keeping the order
    report_files = list(OrderedDict.fromkeys([report_file for report_file in report_files]))
    job.report_files = report_files
    multiqc_files = list(OrderedDict.fromkeys([multiqc_file for multiqc_file in multiqc_files]))
    job.multiqc_files = multiqc_files
    removable_files = list(OrderedDict.fromkeys([removable_file for removable_file in removable_files]))
    job.removable_files = removable_files
    modules = list(OrderedDict.fromkeys([module for module in modules]))
    job.modules = modules
    samples = list(OrderedDict.fromkeys([sample for sample in samples]))
    job.samples = samples

    if input_dependency:
        job.input_files=input_dependency
    if output_dependency:
        job.output_files=output_dependency

    # Merge commands
    job.command = " | \\\n".join([job_item.command for job_item in jobs])

    return job
