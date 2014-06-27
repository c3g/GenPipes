#!/usr/bin/env python

# Python Standard Modules
import datetime
import os

# MUGQIC Modules
from config import *

# Output comment separator line
separator_line = "#" + "-" * 79

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
    def print_header(self, pipeline):
        print(
"""#!/bin/bash

{separator_line}
# {pipeline.__class__.__name__} Torque Job Submission Bash script
# Created on: {datetime}
# Steps:
{steps}
{separator_line}

OUTPUT_DIR={pipeline.output_dir}
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/{pipeline.__class__.__name__}_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR"""
            .format(
                separator_line=separator_line,
                pipeline=pipeline,
                steps="\n".join(["#   " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "" if step.jobs else "... skipping") for step in pipeline.step_range]),
                datetime=datetime.datetime.now()
            )
        )

    def submit(self, pipeline):
        self.print_header(pipeline)
        for step in pipeline.step_range:
            if step.jobs:
                print("""
{separator_line}
# STEP: {step.name}
{separator_line}
STEP={step.name}
mkdir -p $JOB_OUTPUT_DIR/$STEP""".format(separator_line=separator_line, step=step)
                )
                for job in step.jobs:
                    dependency_jobs = ":".join(["$" + dependency_job.id for dependency_job in job.dependency_jobs])
                    print("""
{separator_line}
# {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
JOB_DEPENDENCIES={dependency_jobs}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH""".format(
                            job=job,
                            dependency_jobs=dependency_jobs,
                            separator_line=separator_line
                        )
                    )

                    cmd = \
"""echo "rm -f $JOB_DONE && \\
{job.command_with_modules} && \\
MUGQIC_STATE=\$PIPESTATUS && echo MUGQICexitStatus:\$MUGQIC_STATE && \\
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi && exit \$MUGQIC_STATE" | \\
""".format(job=job)

                    cmd += \
                        config.param(step.name, 'clusterSubmitCmd') + " " + \
                        config.param(step.name, 'clusterOtherArg') + " " + \
                        config.param(step.name, 'clusterWorkDirArg') + " $OUTPUT_DIR " + \
                        config.param(step.name, 'clusterOutputDirArg') + " $JOB_OUTPUT " + \
                        config.param(step.name, 'clusterJobNameArg') + " $JOB_NAME " + \
                        config.param(step.name, 'clusterWalltime') + " " + \
                        config.param(step.name, 'clusterQueue') + " " + \
                        config.param(step.name, 'clusterCPU')
                    if dependency_jobs:
                        cmd += " " + config.param(step.name, 'clusterDependencyArg') + "$JOB_DEPENDENCIES"
                    cmd += " " + config.param(step.name, 'clusterSubmitCmdSuffix')

                    if config.param(step.name, 'clusterCmdProducesJobId'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST"

                    print cmd

class BatchScheduler(Scheduler):
    def submit(self, pipeline):
        for job in pipeline.jobs:
            print(job.command_with_modules)

class DaemonScheduler(Scheduler):
    def submit(self, pipeline):
        for job in pipeline.jobs:
            print("daemon(" + job.command_with_modules + ")")
