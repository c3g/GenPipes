#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules

def create_scheduler(type):
    if type == "torque":
        return TorqueScheduler()
    elif type == "batch":
        return BatchScheduler()
    elif type == "daemon":
        return DaemonScheduler()
    else:
        raise Exception("Error: scheduler type \"" + type + "\" is invalid!")

class Scheduler:
    def submit(self, pipeline):
        # Needs to be defined in scheduler child class
        raise NotImplementedError

class TorqueScheduler(Scheduler):
    def submit(self, pipeline):
        self.print_header(pipeline)
        for job in pipeline.jobs:
            print("echo(\"\n" + job.command_with_modules + " \\\n | qsub -N " + job.name + " depend=afterok:" +
                ":".join(["$" + dependency_job.id for dependency_job in job.dependency_jobs])) + "\")"

    def print_header(self, pipeline):
        print(
"""#!/bin/bash

OUTPUT_DIR=%s
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%%FT%%H.%%M.%%S`
JOB_LIST=$JOB_OUTPUT_DIR/%s_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
""" % (pipeline.output_dir, pipeline.__class__.__name__))
        

class BatchScheduler(Scheduler):
    def submit(self, pipeline):
        for job in pipeline.jobs:
            print(job.command_with_modules)

class DaemonScheduler(Scheduler):
    def submit(self, pipeline):
        for job in pipeline.jobs:
            print("daemon(" + job.command_with_modules + ")")
