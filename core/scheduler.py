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
    else:
        raise Exception("Error: scheduler type \"" + type + "\" is invalid!")

class Scheduler:
    def submit(self, pipeline):
        # Needs to be defined in scheduler child class
        raise NotImplementedError

    def print_header(self, pipeline):
        print(
"""#!/bin/bash

{separator_line}
# {pipeline.__class__.__name__} {scheduler.__class__.__name__} Job Submission Bash script
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
                scheduler=self,
                steps="\n".join(["#   " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "" if step.jobs else "... skipping") for step in pipeline.step_range]) + \
                "\n#   TOTAL: " + str(len(pipeline.jobs)) + " job" + ("s" if len(pipeline.jobs) > 1 else "" if pipeline.jobs else "... skipping"),
                datetime=datetime.datetime.now()
            )
        )

    def print_step(self, step):
        print("""
{separator_line}
# STEP: {step.name}
{separator_line}
STEP={step.name}
mkdir -p $JOB_OUTPUT_DIR/$STEP""".format(separator_line=separator_line, step=step)
        )

class TorqueScheduler(Scheduler):
    def submit(self, pipeline):
        self.print_header(pipeline)
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    if job.dependency_jobs:
                        # Chunk JOB_DEPENDENCIES on multiple lines to avoid lines too long
                        max_dependencies_per_line = 50
                        dependency_chunks = [job.dependency_jobs[i:i + max_dependencies_per_line] for i in range(0, len(job.dependency_jobs), max_dependencies_per_line)]
                        job_dependencies = "JOB_DEPENDENCIES=" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunks[0]])
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += "\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES:" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunk])
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    print("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line
                        )
                    )

                    cmd = \
"""echo "rm -f $JOB_DONE && \\
{job.command_with_modules}
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \\
""".format(job=job)

                    cmd += \
                        config.param(step.name, 'cluster_submit_cmd') + " " + \
                        config.param(step.name, 'cluster_other_arg') + " " + \
                        config.param(step.name, 'cluster_work_dir_arg') + " $OUTPUT_DIR " + \
                        config.param(step.name, 'cluster_output_dir_arg') + " $JOB_OUTPUT " + \
                        config.param(step.name, 'cluster_job_name_arg') + " $JOB_NAME " + \
                        config.param(step.name, 'cluster_walltime') + " " + \
                        config.param(step.name, 'cluster_queue') + " " + \
                        config.param(step.name, 'cluster_cpu')
                    if job.dependency_jobs:
                        cmd += " " + config.param(step.name, 'cluster_dependency_arg') + "$JOB_DEPENDENCIES"
                    cmd += " " + config.param(step.name, 'cluster_submit_cmd_suffix')

                    if config.param(step.name, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST"

                    print cmd

class BatchScheduler(Scheduler):
    def submit(self, pipeline):
        self.print_header(pipeline)
        print("SEPARATOR_LINE=`seq -s - 80 | sed 's/[0-9]//g'`")
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    print("""
{separator_line}
# JOB: {job.name}
{separator_line}
JOB_NAME={job.name}
JOB_DONE={job.done}
printf "\\n$SEPARATOR_LINE\\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H.%M.%S`" && \\
rm -f $JOB_DONE && \\
{command_with_modules}
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H.%M.%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi""".format(
                            job=job,
                            separator_line=separator_line,
                            command_with_modules=re.sub(r"\\(.)", r"\1", job.command_with_modules)
                        )
                    )
