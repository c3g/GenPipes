#!/usr/bin/env python

# Python Standard Modules
import collections
import os

# MUGQIC Modules
from config import *

class Job:

    def __init__(self, config, input_files, output_files, module_entries = []):
        self._config = config
        # Remove undefined input/output files if any
        self._input_files = filter(None, input_files)
        self._output_files = filter(None, output_files)

        # Retrieve modules from config, removing duplicates but keeping the order
        self._modules = list(collections.OrderedDict.fromkeys([self.config.param(section, option) for section, option in module_entries]))

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
    def config(self):
        return self._config

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

#config = Config("/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/dnaSeq.abacus.ini")
#job = Job(
#    config,
#    ["/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/dnaSeq.abacus.ini", "/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/project.nanuq.csv"],
#    ["/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/output_file"],
#    [
#        ["trim", "moduleVersion.trimmomatic"],
#        ["aln", "moduleVersion.bwa"],
#        ["recalibration", "moduleVersion.picard"]
#    ])
#
#job.command = "ls -l"
#print(job.command_with_modules)
#print(job.is_up2date)
#job.show()
