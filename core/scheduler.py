#!/usr/bin/env python

# Python Standard Modules
import json
import os

# MUGQIC Modules
from config import *

# Output comment separator line
separator_line = "#" + "-" * 79

def create_scheduler(type):
    if type == "pbs":
        return PBSScheduler()
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

    def print_header(self, pipeline):
        print("""\
#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

{separator_line}
# {pipeline.__class__.__name__} {scheduler.__class__.__name__} Job Submission Bash script
# Version: {pipeline.version}
# Created on: {pipeline.timestamp}
# Steps:
{steps}
{separator_line}"""
            .format(
                separator_line=separator_line,
                pipeline=pipeline,
                scheduler=self,
                steps="\n".join(["#   " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "" if step.jobs else "... skipping") for step in pipeline.step_range]) + \
                "\n#   TOTAL: " + str(len(pipeline.jobs)) + " job" + ("s" if len(pipeline.jobs) > 1 else "" if pipeline.jobs else "... skipping")
            )
        )

        if pipeline.jobs:
            print(
"""
OUTPUT_DIR={pipeline.output_dir}
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/{pipeline.__class__.__name__}_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
"""
                .format(
                    pipeline=pipeline,
                )
            )

    def print_step(self, step):
        print("""
{separator_line}
# STEP: {step.name}
{separator_line}
STEP={step.name}
mkdir -p $JOB_OUTPUT_DIR/$STEP
""".format(separator_line=separator_line, step=step)
        )

class PBSScheduler(Scheduler):
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
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << '{limit_string}'
{job.command_with_modules}
{limit_string}
)""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line,
                            limit_string=os.path.basename(job.done)
                        )
                    )

                    cmd = """\
echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \\
""".format(job=job)

                    # Cluster settings section must match job name prefix before first "."
                    # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                    job_name_prefix = job.name.split(".")[0]
                    cmd += \
                        config.param(job_name_prefix, 'cluster_submit_cmd') + " " + \
                        config.param(job_name_prefix, 'cluster_other_arg') + " " + \
                        config.param(job_name_prefix, 'cluster_work_dir_arg') + " $OUTPUT_DIR " + \
                        config.param(job_name_prefix, 'cluster_output_dir_arg') + " $JOB_OUTPUT " + \
                        config.param(job_name_prefix, 'cluster_job_name_arg') + " $JOB_NAME " + \
                        config.param(job_name_prefix, 'cluster_walltime') + " " + \
                        config.param(job_name_prefix, 'cluster_queue') + " " + \
                        config.param(job_name_prefix, 'cluster_cpu')
                    if job.dependency_jobs:
                        cmd += " " + config.param(job_name_prefix, 'cluster_dependency_arg') + "$JOB_DEPENDENCIES"
                    cmd += " " + config.param(job_name_prefix, 'cluster_submit_cmd_suffix')

                    if config.param(job_name_prefix, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST\n"

                    print cmd

        # Check cluster maximum job submission
        cluster_max_jobs = config.param('DEFAULT', 'cluster_max_jobs', type='posint', required=False)
        if cluster_max_jobs and len(pipeline.jobs) > cluster_max_jobs:
            log.warning("Number of jobs: " + str(len(pipeline.jobs)) + " > Cluster maximum number of jobs: " + str(cluster_max_jobs) + "!")

class BatchScheduler(Scheduler):
    def submit(self, pipeline):
        self.print_header(pipeline)
        if pipeline.jobs:
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
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \\
rm -f $JOB_DONE && \\
{job.command_with_modules}
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi
""".format(
                            job=job,
                            separator_line=separator_line
                        )
                    )

class DaemonScheduler(Scheduler):
    def submit(self, pipeline):
        print self.json(pipeline)

    def json(self, pipeline):
        return json.dumps(
            {'pipeline': [
                {step.name: [
                    {
                        "job_name": job.name,
                        "job_id": job.id,
                        "job_command": job.command_with_modules,
                        "job_dependencies": [dependency_job.id for dependency_job in job.dependency_jobs],
                        "job_cluster_options": 
                            # Cluster settings section must match job name prefix before first "."
                            # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                            config.param(job.name.split(".")[0], 'cluster_submit_cmd') + " " + \
                            config.param(job.name.split(".")[0], 'cluster_other_arg') + " " + \
                            config.param(job.name.split(".")[0], 'cluster_work_dir_arg') + " " + pipeline.output_dir + " " + \
                            config.param(job.name.split(".")[0], 'cluster_output_dir_arg') + " " + os.path.join(pipeline.output_dir, "job_output", step.name, job.name + ".o") + " " + \
                            config.param(job.name.split(".")[0], 'cluster_job_name_arg') + " " + job.name + " " + \
                            config.param(job.name.split(".")[0], 'cluster_walltime') + " " + \
                            config.param(job.name.split(".")[0], 'cluster_queue') + " " + \
                            config.param(job.name.split(".")[0], 'cluster_cpu')
                    } for job in step.jobs
                ]} for step in pipeline.step_range
            ]}, indent=4)
