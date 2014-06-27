#!/usr/bin/env python

# Python Standard Modules
import argparse
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
        self._step_map = []
        for step_dict in self.step_dict_map:
            step_name = step_dict['name'].__name__
            if self.steps_by_name(step_name):
                raise Exception("Error: step name \"" + step_name + "\" is not unique for this pipeline!")
            else:
                self.step_map.append(Step(step_dict['name'], step_dict['loop'] if 'loop' in step_dict else None))

        if re.search("^\d+([,-]\d+)*$", args.steps):
            self._steps = self.steps_by_range(args.steps)
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
    def step_dict_map(self):
        # Needs to be defined in pipeline child class
        raise NotImplementedError

    @property
    def step_map(self):
        return self._step_map

    @property
    def steps(self):
        return self._steps

    @property
    def jobs(self):
        jobs = []
        for step in self.steps:
            jobs.extend(step.jobs)
        return jobs

    def steps_by_name(self, name):
        return [step for step in self.step_map if step.name == name]

    def steps_by_range(self, range):
        return [self.step_map[i - 1] for i in parse_range(range)]

    def dependency_jobs(self, current_job):
        dependency_jobs = []
        for step in self.steps:
            for step_job in step.jobs:
                # If current job input files intersect with step job output files, step job is a dependency
                if set(current_job.input_files).intersection(set(step_job.output_files)):
                    dependency_jobs.append(step_job)
        return dependency_jobs

    def create_jobs(self):
        for step in self.steps:
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

    def submit_jobs(self):
        self.scheduler.submit(self)

    def show(self):
        print('Pipeline: ' + self.__class__.__name__)
        print('output_dir: ' + self.output_dir)
        print('scheduler: ' + self.scheduler.__class__.__name__)
        print('Steps:')
        for step in self.steps:
            print("Step: " + step.name)
            for job in step.jobs:
                print(job.id)
                for dependency_job in job.dependency_jobs:
                    print("  Dependency " + dependency_job.id + ", command: " + dependency_job.command)
                print("  Input files:" + ",".join(job.input_files))
                print("  Command: " + job.command)
                print("  Output files:" + ",".join(job.output_files))
                print()

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

    def __init__(self, step_dict_map):
        # Create ArgumentParser with numbered step list as epilog
        argparse.ArgumentParser.__init__(self, formatter_class=argparse.RawDescriptionHelpFormatter, epilog="Steps:\n" + "\n".join([str(idx + 1) + "- " + step_dict['name'].__name__ for idx, step_dict in enumerate(step_dict_map)]))

        self.add_argument("-c", "--config", help="config INI-style file", type=file, required=True)
        self.add_argument("-s", "--steps", help="step range e.g. '1-5', '3,6,7', '2,4-8'", required=True)
        self.add_argument("-o", "--output-dir", help="output directory (default: current)", default=os.getcwd())
        self.add_argument("-j", "--job-scheduler", help="job scheduler type (default: torque)", choices=["torque", "batch", "daemon"], default="torque")
        self.add_argument("-f", "--force", help="force creation of jobs even if up to date (default: false)", action="store_true")
        self.add_argument("-l", "--log", help="log level (default: info)", choices=["debug", "info", "warning", "error", "critical"], default="info")

    def parse_args(self):
        args = argparse.ArgumentParser.parse_args(self)
        logging.basicConfig(level=getattr(logging, args.log.upper()))
        config.parse_file(args.config)
        return args
