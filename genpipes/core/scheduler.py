################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################
import json
import logging
import os
import re
import stat
import sys
import tempfile
import textwrap
from uuid import uuid4

from ..utils import utils

from .config import global_conf

# Output comment separator line
separator_line = "#" + "-" * 79

logger = logging.getLogger(__name__)


def create_scheduler(s_type, config_files, container=None, genpipes_file=None):
    if s_type == "pbs":
        return PBSScheduler(config_files, container=container, genpipes_file=genpipes_file)
    elif s_type == "batch":
        return BatchScheduler(config_files, container=container, genpipes_file=genpipes_file)
    elif s_type == "daemon":
        return DaemonScheduler(config_files, genpipes_file=genpipes_file)
    elif s_type == "slurm":
        return SlurmScheduler(config_files, container=container, genpipes_file=genpipes_file)
    else:
        raise Exception("Error: scheduler type \"" + s_type + "\" is invalid!")


class Scheduler:
    """
    Common Scheduler
    """
    def __init__(self, config_files, container=None, genpipes_file=None, **kwargs):

        self.name = 'generic'
        self._config_files = config_files
        self._container = container
        self._host_cvmfs_cache = None
        self._cvmfs_cache = None
        self._bind = None
        self._submit_cmd = None
        if genpipes_file is None:
            self.genpipes_file = sys.stdout
        else:
            self.genpipes_file = genpipes_file
            if genpipes_file is not sys.stdout:
                # make sure it is user and group executable
                st = os.stat(genpipes_file.name)
                os.chmod(genpipes_file.name, st.st_mode | stat.S_IEXEC | stat.S_IXGRP)
        self.write = self.genpipes_file.write
        self.flush = self.genpipes_file.flush

    @property
    def submit_cmd(self):
        if self._submit_cmd is None:
            raise NotImplementedError('_submit_cmd needs to be implemented for {} class'.format(self.__class__))
        return self._submit_cmd

    def walltime(self, job_name_prefix):
        raise NotImplementedError

    def memory(self, job_name_prefix):
        raise NotImplementedError

    def gpu_type(self, job_name_prefix):
        gpu_type = global_conf.global_get(job_name_prefix, 'cluster_gpu_type', required=False)
        return ''.join(gpu_type.split())

    def gpu(self, job_name_prefix):
        return global_conf.global_get(job_name_prefix, 'cluster_gpu', required=False)

    def dependency_arg(self, job_name_prefix):
        # be careful "after" is a subset of the other stings and must be at the end of the list.
        supported = ['afterany', 'afternotok', 'afterok', 'after']
        dep_str = global_conf.global_get(job_name_prefix, 'cluster_dependency_arg')
        for condition in supported:
            if condition in dep_str:
                return condition
        raise ValueError(f'{dep_str} not part of cluster_dependency_arg supported value {supported}')

    def cpu(self, job_name_prefix):
        cpu_str = global_conf.global_get(job_name_prefix, 'cluster_cpu', required=True)
        try:
            if "ppn" in cpu_str or '-c' in cpu_str:
                # to be back compatible
                cpu = re.search(r"(ppn=|-c\s)([0-9]+)", cpu_str).groups()[1]
            else:
                cpu = re.search(r"[0-9]+", cpu_str).group()
        except AttributeError:
            raise ValueError(f'"{cpu_str}" is not a valid entry for "cluster_cpu" ({job_name_prefix=})')
        return cpu

    def node(self, job_name_prefix):
        # code run on 1 node by default
        node = 1
        node_str = global_conf.global_get(job_name_prefix, 'cluster_node', required=False)

        cpu_str = None
        if not node_str:
            cpu_str = global_conf.global_get(job_name_prefix, 'cluster_cpu', required=False)
        if cpu_str:
            try:
                if "nodes" in cpu_str or '-N' in cpu_str:
                    # to be back compatible
                    return re.search(r"(nodes=|-N\s*)([0-9]+)",cpu_str).groups()[1]
            except AttributeError as exc:
                raise ValueError(f'"{cpu_str}" is not a valid entry for "cluster_cpu"') from exc
        try:
            return re.search("[0-9]+", node_str).group()
        except AttributeError:
            return node

    def submit(self, pipeline):
        # Needs to be defined in scheduler child class
        raise NotImplementedError

    @property
    def container(self):
        return self._container

    @property
    def host_cvmfs_cache(self):

        if self._host_cvmfs_cache is None:

            self._host_cvmfs_cache = global_conf.global_get("container", 'host_cvmfs_cache',
                                                     required=False, param_type="string")
            if not self._host_cvmfs_cache:

                tmp_dir = global_conf.global_get("DEFAULT", 'tmp_dir', required=True)
                tmp_dir = utils.expandvars(tmp_dir)

                if not tmp_dir:
                    tmp_dir = None

                self._host_cvmfs_cache = tempfile.mkdtemp(prefix="genpipes_cvmfs_", dir=tmp_dir)

        return self._host_cvmfs_cache

    @property
    def cvmfs_cache(self):

        if self._cvmfs_cache is None:
            self._cvmfs_cache = global_conf.global_get("container", 'cvmfs_cache', required=False, param_type="string")
            if not self._cvmfs_cache:
                self._cvmfs_cache = "/cvmfs-cache"

        return self._cvmfs_cache

    @property
    def bind(self):
        if self._bind is None:
            self._bind = global_conf.global_get("container", 'bind_list', required=False, param_type='list')

            if not self._bind:
                self._bind = ['/tmp', '/home']
        return self._bind

    @property
    def disable_modulercfile(self):
        if self.container:
            return 'unset MODULERCFILE'
        return ""

    @property
    def container_line(self):
        if self.container:
            if self.container.type == 'docker':
                v_opt = ' -v {}:{}'.format(self.host_cvmfs_cache, self.cvmfs_cache)
                for b in self.bind:
                    v_opt += ' -v {0}:{0}'.format(b)
                network = " --network host"
                user = " --user $UID:$GROUPS"

                return (f'docker run --env-file <( env| cut -f1 -d= ) --privileged {network} {user} {v_opt} {self.container.name} ')

            elif self.container.type == 'singularity':
                b_opt = f' -B {self.host_cvmfs_cache}:{self.cvmfs_cache}'
                for b in self.bind:
                    b_opt += ' -B {0}:{0}'.format(b)

                return (f"singularity run {b_opt} {self.container.name}   ")

            elif self.container.type == 'wrapper':

                return f"{self.container.name} "
        else:

            return ""


    def fail_on_pattern(self, job_name_prefix):
        pattern = global_conf.global_get(job_name_prefix, 'fail_on_pattern',required=False)

        if not pattern:
            return ("", "")
        else:
            tmp_dir = global_conf.global_get("DEFAULT", 'tmp_dir', required=True)

            append_command = f" | tee {tmp_dir}/${{JOB_NAME}}_${{TIMESTAMP}}.o "
            test_condition = """
grep {pattern} {tmp_dir}/${{JOB_NAME}}_${{TIMESTAMP}}.o
NO_PROBLEM_IN_LOG=\\$?

  if [[  \\$NO_PROBLEM_IN_LOG == 0 ]] ; then
   echo {pattern} found in job log, forcing error 
   GenPipes_STATE=74
fi
""".format(pattern=pattern, tmp_dir=tmp_dir)

            return (append_command, test_condition)



    def print_header(self, pipeline,shebang='/bin/bash'):

        self.genpipes_file.write("""#!{shebang}
# Exit immediately on error
{scheduler.disable_modulercfile}
set -eu -o pipefail

{separator_line}
# {pipeline.__class__.__name__} {scheduler.name} Job Submission Bash script
# Version: {pipeline.version}
# Created on: {pipeline.timestamp}
# Steps:
{steps}
{separator_line}
"""
            .format(
                shebang=shebang,
                separator_line=separator_line,
                pipeline=pipeline,
                scheduler=self,
                steps="\n".join(["#   " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "" if step.jobs else "... skipping") for step in pipeline.step_to_execute]) + \
                "\n#   TOTAL: " + str(len(pipeline.jobs)) + " job" + ("s" if len(pipeline.jobs) > 1 else "" if pipeline.jobs else "... skipping")
            )
        )

        if pipeline.jobs:
            self.genpipes_file.write("""
OUTPUT_DIR={pipeline.output_dir}
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP={pipeline.timestamp}
JOB_LIST=$JOB_OUTPUT_DIR/{pipeline.__class__.__name__}{protocol}.job_list.$TIMESTAMP
export CONFIG_FILES="{config_files}"
mkdir -p $OUTPUT_DIR
"""
                .format(
                    pipeline=pipeline,
                    protocol="." + pipeline.protocol if pipeline.protocol else "",
                    config_files=",".join([c.name for c in self._config_files])
                )
            )
        if pipeline.json:
            json_files = []
            for step in pipeline.step_to_execute:
                for job in step.jobs:
                    for sample in job.samples:
                        json_files.append(os.path.join(pipeline.output_dir, "json", sample.json_file))
            json_files = list(set(json_files))
            for j_file in json_files:
                self.genpipes_file.write("""
                    sed -i "s/\\"submission_date\\": \\"\\",/\\"submission_date\\": \\"$TIMESTAMP\\",/" {file}
""".format(file=j_file))

        ## Print a copy of sample JSONs for the genpipes dashboard
        if pipeline.json and pipeline.portal_output_dir != "":
            copy_commands = []
            test_copy_commands = []
            for i, sample in enumerate(pipeline.sample_list):
                unique_uuid = str(uuid4())
                input_file = pipeline.sample_paths[i]
                output_file = os.path.join(pipeline.portal_output_dir, '$USER.' + sample.name + '.' + unique_uuid + '.json')
                #test_output_file = os.path.join("/lb/project/mugqic/analyste_dev/portal_test_dir/", '$USER.' + sample.name + '.' + unique_uuid + '.json')
                copy_commands.append("cp \"{input_file}\" \"{output_file}\"".format(
                    input_file=input_file, output_file=output_file))
                #test_copy_commands.append("cp \"{input_file}\" \"{output_file}\"".format(
                    #input_file=input_file, output_file=test_output_file))
            self.genpipes_file.write(textwrap.dedent("""
                #------------------------------------------------------------------------------
                # Print a copy of sample JSONs for the genpipes dashboard
                #------------------------------------------------------------------------------
                {copy_commands}
            """).format(copy_commands='\n'.join(copy_commands)))

        if pipeline.jobs:
            self.genpipes_file.write("cd $OUTPUT_DIR")


    def print_step(self, step):
        self.genpipes_file.write("""
{separator_line}
# STEP: {step.name}
{separator_line}
STEP={step.name}
mkdir -p $JOB_OUTPUT_DIR/$STEP
""".format(separator_line=separator_line, step=step)
                                 )

    def job2json(self, pipeline, step, job, job_status):
        if not pipeline.json:
            return ""

        json_file_list = ",".join([os.path.join(pipeline.output_dir, "json", sample.json_file) for sample in job.samples])
        return """\
{job2json_script} \\
  -u \\"$USER\\" \\
  -c \\"{config_files}\\" \\
  -s \\"{step.name}\\" \\
  -j \\"$JOB_NAME\\" \\
  -d \\"$JOB_DONE\\" \\
  -l \\"$JOB_OUTPUT\\" \\
  -o \\"{jsonfiles}\\" \\
  -f {status} {command_separator}
""".format(
            job2json_script="job2json.py",
            module_python=global_conf.global_get('DEFAULT', 'module_python'),
            step=step,
            jsonfiles=json_file_list,
            config_files=",".join([ os.path.abspath(c.name) for c in self._config_files ]),
            status=job_status,
            command_separator="&&" if (job_status=='\\"running\\"') else ""
        ) if json_file_list else ""

    def job2json_project_tracking(self, pipeline, job, job_status):
        if not pipeline.project_tracking_json:
            return ""

        pipeline_output_dir = pipeline.output_dir
        json_folder = os.path.join(pipeline_output_dir, "json")
        timestamp = pipeline.timestamp
        try:
            json_outfile = os.path.join(json_folder, f"{pipeline.__class__.__name__}.{pipeline.protocol}_{timestamp}.json")
        except AttributeError:
            json_outfile = os.path.join(json_folder, f"{pipeline.__class__.__name__}_{timestamp}.json")

        # The project tracking json file is exported here as a variable for use by job2json_project_tracking.py.
        # Avoids restarts in steps that used to reference the json file name with timestamp in the name.
        return """\
{job2json_project_tracking_script} \\
  -s \\"{samples}\\" \\
  -r \\"{readsets}\\" \\
  -j \\"{job_name}\\" \\{metrics}
  -o \\"{json_outfile}\\" \\
  -f {status}
export PT_JSON_OUTFILE=\\"{json_outfile}\\" {command_separator}
""".format(
            module_python=global_conf.global_get('DEFAULT', 'module_python'),
            job2json_project_tracking_script="job2json_project_tracking.py",
            samples=",".join([sample.name for sample in job.samples]),
            readsets=",".join([readset.name for readset in job.readsets]),
            job_name=job.name,
            metrics=('\n  -m \\"' + ','.join(job.metrics) + '\\" \\') if job.metrics else '',
            json_outfile=json_outfile,
            status=job_status,
            command_separator="&&" if (job_status=='\\"RUNNING\\"') else ""
        ) if json_outfile else ""


class PBSScheduler(Scheduler):
    """
    PBS Scheduler
    """
    def __init__(self, *args, **kwargs):
        super(PBSScheduler, self).__init__(*args, **kwargs)
        # should be fed in the arguments but hey lets do that first.
        self.config = global_conf
        self._submit_cmd = 'qsub'
        self.name = 'PBS'

    def walltime(self, job_name_prefix):
        walltime = global_conf.global_get(job_name_prefix, 'cluster_walltime')
        # force the DD-HH:MM[:00] format to HH:MM[:00]
        time = utils.time_to_datetime(walltime)
        sec = int(time.seconds % 60)
        minutes = int(((time.seconds - sec) / 60) % 60)
        hours = int((time.seconds - sec - 60 * minutes) / 3600 + time.days * 24)
        return '-l walltime={:02d}:{:02d}:{:02d}'.format(hours, minutes, sec)

    def dependency_arg(self, job_name_prefix):
        condition = super().dependency_arg(job_name_prefix)
        return f'-W depend={condition}:'

    def memory(self, job_name_prefix, adapt=None, info=False):
        mem_str = global_conf.global_get(job_name_prefix, 'cluster_mem', required=False)
        try:
            mem = re.search("[0-9]+[a-zA-Z]*", mem_str).group()
        except AttributeError:
            return " "

        if adapt is not None:
            return ''

        if 'per' in mem_str.lower() and 'cpu' in mem_str.lower():
            option = '-l pmem='
            per_core = True
        else:
            option = '-l mem='
            per_core = False

        if info:
            return per_core, mem

        return f"{option}{mem}"

    def cpu(self, job_name_prefix, adapt=None):
        cpu = int(super().cpu(job_name_prefix))

        mem_info = self.memory(job_name_prefix, info=True)
        if adapt and mem_info:
            adapt_mem = int(re.search("[0-9]+", adapt).group())
            if 'G' in adapt:
                adapt_mem = adapt_mem * 1024

            [per_cpu, mem_str] = mem_info
            mem = int(re.search("[0-9]+", mem_str).group())
            if 'G' in mem_str:
                mem = mem * 1024
            if per_cpu:
                mem = mem * cpu

            import math
            cpu_ = math.ceil(mem/adapt_mem)
            cpu = max(cpu, cpu_)

        node = self.node(job_name_prefix)
        gpu = self.gpu(job_name_prefix)
        if gpu:
            return f"-l nodes={node}:ppn={cpu}:gpu{gpu}"
        else:
            return f"-l nodes={node}:ppn={cpu}"

    def submit(self, pipeline):
        self.print_header(pipeline)
        self.genpipes_file.flush()
        for step in pipeline.step_to_execute:
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

                    job_name_prefix = job.name.split(".")[0]
                    config_step_wrapper = global_conf.global_get(job_name_prefix, 'step_wrapper', required=False)

                    #sleepTime = random.randint(10, 100)
                    self.genpipes_file.write("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.sh
cat << '{limit_string}' > $COMMAND
{job.command_with_modules}
{limit_string}
chmod 755 $COMMAND
""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line,
                            limit_string=os.path.basename(job.done)
                        )
                    )

                    self.genpipes_file.flush()

                    cmd = """\
echo "rm -f $JOB_DONE && {job2json_project_tracking_start} {job2json_start} {step_wrapper} {container_line} $COMMAND {fail_on_pattern0}
GenPipes_STATE=\\$PIPESTATUS
echo GenPipesExitStatus:\\$GenPipes_STATE
{job2json_end}
{job2json_project_tracking_end}
{fail_on_pattern1}
if [ \\$GenPipes_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \\$GenPipes_STATE" | \\
""".format(
                        container_line=self.container_line,
                        job2json_project_tracking_start=self.job2json_project_tracking(pipeline, job, '\\"RUNNING\\"'),
                        job2json_project_tracking_end=self.job2json_project_tracking(pipeline, job, '\\$GenPipes_STATE'),
                        job2json_start=self.job2json(pipeline, step, job, '\\"running\\"'),
                        job2json_end=self.job2json(pipeline, step, job, '\\$GenPipes_STATE'),
                        step_wrapper=config_step_wrapper,
                        fail_on_pattern0=self.fail_on_pattern(job_name_prefix)[0],
                        fail_on_pattern1=self.fail_on_pattern(job_name_prefix)[1]
                    )
                        #sleep_time=sleepTime

                    # Cluster settings section must match job name prefix before first "."
                    # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                    job_name_prefix = job.name.split(".")[0]
                    cmd += \
                        self.submit_cmd + " " + \
                        global_conf.global_get(job_name_prefix, 'cluster_other_arg') + " " + \
                        global_conf.global_get(job_name_prefix, 'cluster_work_dir_arg') + " $OUTPUT_DIR " + \
                        global_conf.global_get(job_name_prefix, 'cluster_output_dir_arg') + " $JOB_OUTPUT " + \
                        global_conf.global_get(job_name_prefix, 'cluster_job_name_arg') + " $JOB_NAME " + \
                        self.walltime(job_name_prefix) + " " + \
                        self.memory(job_name_prefix, adapt=pipeline.force_mem_per_cpu) + " " + \
                        self.cpu(job_name_prefix, adapt=pipeline.force_mem_per_cpu) + " " + \
                        global_conf.global_get(job_name_prefix, 'cluster_queue') + " "

                    if job.dependency_jobs:
                        cmd += " " + self.dependency_arg(job_name_prefix) + "$JOB_DEPENDENCIES"
                    cmd += " " + global_conf.global_get(job_name_prefix, 'cluster_submit_cmd_suffix')

                    if global_conf.global_get(job_name_prefix, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST\n"

                    self.genpipes_file.write(cmd)
                    self.genpipes_file.flush()

        # Check cluster maximum job submission
        cluster_max_jobs = global_conf.global_get('DEFAULT', 'cluster_max_jobs', param_type='posint', required=False)
        if cluster_max_jobs and len(pipeline.jobs) > cluster_max_jobs:
            logging.warning(f"Number of jobs: {str(len(pipeline.jobs))} > Cluster maximum number of jobs: {str(cluster_max_jobs)}!")


class BatchScheduler(Scheduler):
    """
    No Scheduler: batch mode
    """
    def __init__(self, *args, **kwargs):
        super(BatchScheduler, self).__init__(*args, **kwargs)
        self.name = 'Batch'

    def job2json(self, pipeline, step, job, job_status):
        if not pipeline.json:
            return ""

        json_file_list = ",".join([os.path.join(pipeline.output_dir, "json", sample.json_file) for sample in job.samples])
        return """\
module load {module_python}
{job2json_script} \\
  -u \"$USER\" \\
  -c \"{config_files}\" \\
  -s \"{step.name}\" \\
  -j \"$JOB_NAME\" \\
  -d \"$JOB_DONE\" \\
  -l \"$JOB_OUTPUT\" \\
  -o \"{jsonfiles}\" \\
  -f {status}
module unload {module_python} {command_separator}
""".format(
            job2json_script="job2json.py",
            module_python=global_conf.global_get('DEFAULT', 'module_python'),
            step=step,
            jsonfiles=json_file_list,
            config_files=",".join([ os.path.abspath(c.name) for c in self._config_files ]),
            status=job_status,
            command_separator="&&" if (job_status=='"running"') else ""
        ) if json_file_list else ""

    def submit(self, pipeline):
        logger.info('\n\t To run the script use: \n\t"{}  ./<command>.sh"'.format(
            self.container_line))
        self.print_header(pipeline)
        if pipeline.jobs:
            self.genpipes_file.write("SEPARATOR_LINE=`seq -s - 80 | sed 's/[0-9]//g'`\n")
        for step in pipeline.step_to_execute:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    job_name_prefix = job.name.split(".")[0]
                    config_step_wrapper = global_conf.global_get(job_name_prefix, 'step_wrapper', required=False)

                    self.genpipes_file.write("""
{separator_line}
# JOB: {job.name}
{separator_line}
JOB_NAME={job.name}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.sh
cat << '{limit_string}' > $COMMAND
#!/bin/bash
set -eu -o pipefail
{job.command_with_modules}
{limit_string}
chmod 755 $COMMAND
printf "\\n$SEPARATOR_LINE\\n"
echo "Begin GenPipes Job $JOB_NAME at `date +%FT%H:%M:%S`" && \\
rm -f $JOB_DONE && {job2json_project_tracking_start} {job2json_start} {step_wrapper} $COMMAND &> $JOB_OUTPUT
GenPipes_STATE=$?
echo "End GenPipes Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo GenPipesExitStatus:$GenPipes_STATE
{job2json_end}
{job2json_project_tracking_end}
if [ $GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $GenPipes_STATE ; fi
""".format(
                            job=job,
                            limit_string=os.path.basename(job.done),
                            separator_line=separator_line,
                            job2json_project_tracking_start=self.job2json_project_tracking(pipeline, job, '\\"RUNNING\\"'),
                            job2json_project_tracking_end=self.job2json_project_tracking(pipeline, job, '\\$GenPipes_STATE'),
                            job2json_start=self.job2json(pipeline, step, job, '\\"running\\"'),
                            job2json_end=self.job2json(pipeline, step, job, '\\$GenPipes_STATE'),
                            step_wrapper=config_step_wrapper
                        )
                    )


class SlurmScheduler(Scheduler):
    """
    Slurm Scheduler
    """
    def __init__(self, *args, **kwargs):
        super(SlurmScheduler, self).__init__(*args, **kwargs)
        self.name = 'SLURM'
        self.config = global_conf
        self._submit_cmd = 'sbatch'

    def walltime(self, job_name_prefix):
        walltime = global_conf.global_get(job_name_prefix, 'cluster_walltime')
        # force the DD-HH:MM[:00] format to HH:MM[:00]
        time = utils.time_to_datetime(walltime)
        sec = int(time.seconds % 60)
        minutes = int(((time.seconds - sec) / 60) % 60)
        hours = int((time.seconds - sec - 60 * minutes) / 3600 + time.days * 24)
        return '--time={:02d}:{:02d}:{:02d}'.format(hours, minutes, sec)

    def gpu(self, job_name_prefix):
        n_gpu = super().gpu(job_name_prefix)
        gpu_type = self.gpu_type(job_name_prefix)
        if gpu_type and n_gpu:
            return f'--gres=gpu:{gpu_type}:{n_gpu}'
        elif n_gpu:
            return f'--gres=gpu:{n_gpu}'
        else:
            return ''

    def dependency_arg(self, job_name_prefix):
        condition = super().dependency_arg(job_name_prefix)
        return f'--depend={condition}:'

    def memory(self, job_name_prefix):
        config_str = 'cluster_mem'
        mem_str = global_conf.global_get(job_name_prefix, config_str, required=False)
        try:
            mem = re.search("[0-9]+[a-zA-Z]*", mem_str).group()
        except AttributeError:
            return " "

        if 'per' in mem_str.lower() and 'cpu' in mem_str.lower():
            option = '--mem-per-cpu'
        else:
            option = '--mem'
        return f"{option} {mem}"

    def cpu(self, job_name_prefix):
        cpu = super().cpu(job_name_prefix)
        return f'-c {cpu}'

    def node(self, job_name_prefix):
        node = super().node(job_name_prefix)
        return f'-N {node}'

    def submit(self, pipeline):
        self.print_header(pipeline)
        self.genpipes_file.flush()

        for step in pipeline.step_to_execute:
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

                    job_name_prefix = job.name.split(".")[0]
                    config_step_wrapper = global_conf.global_get(job_name_prefix, 'step_wrapper', required=False)

                    self.genpipes_file.write("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.sh
cat << '{limit_string}' > $COMMAND
{job.command_with_modules}
{limit_string}
chmod 755 $COMMAND
""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line,
                            limit_string=os.path.basename(job.done)
                    )
                        )
                    self.genpipes_file.flush()

                    cmd = """\
echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \\$SLURM_JOBID
sstat -j \\$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE && {job2json_project_tracking_start} {job2json_start} {step_wrapper} {container_line}  $COMMAND {fail_on_pattern0}
GenPipes_STATE=\\$PIPESTATUS
echo GenPipesExitStatus:\\$GenPipes_STATE
{job2json_end}
{job2json_project_tracking_end}
{fail_on_pattern1}
if [ \\$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \\$SLURM_JOBID
sstat -j \\$SLURM_JOBID.batch
echo '#######################################'
exit \\$GenPipes_STATE" | \\
""".format(
                        job=job,
                        job2json_project_tracking_start=self.job2json_project_tracking(pipeline, job, '\\"RUNNING\\"'),
                        job2json_project_tracking_end=self.job2json_project_tracking(pipeline, job, '\\$GenPipes_STATE'),
                        job2json_start=self.job2json(pipeline, step, job, '\\"running\\"'),
                        job2json_end=self.job2json(pipeline, step, job, '\\$GenPipes_STATE') ,
                        container_line=self.container_line,
                        step_wrapper=config_step_wrapper,
                        fail_on_pattern0=self.fail_on_pattern(job_name_prefix)[0],
                        fail_on_pattern1=self.fail_on_pattern(job_name_prefix)[1]
)

                    # Cluster settings section must match job name prefix before first "."
                    # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                    job_name_prefix = job.name.split(".")[0]
                    cmd += \
                        self.submit_cmd + " " + \
                        global_conf.global_get(job_name_prefix, 'cluster_other_arg') + " " + \
                        global_conf.global_get(job_name_prefix, 'cluster_work_dir_arg') + " $OUTPUT_DIR " + \
                        global_conf.global_get(job_name_prefix, 'cluster_output_dir_arg') + " $JOB_OUTPUT " + \
                        global_conf.global_get(job_name_prefix, 'cluster_job_name_arg') + " $JOB_NAME " + \
                        self.walltime(job_name_prefix) + " " + \
                        self.memory(job_name_prefix) + " " + \
                        self.cpu(job_name_prefix) + " " + \
                        self.node(job_name_prefix) + " " + \
                        self.gpu(job_name_prefix) + " " + \
                        global_conf.global_get(job_name_prefix, 'cluster_queue') + " "

                    if job.dependency_jobs:
                        cmd += " " + self.dependency_arg(job_name_prefix) + "$JOB_DEPENDENCIES"
                    cmd += " " + global_conf.global_get(job_name_prefix, 'cluster_submit_cmd_suffix')

                    if global_conf.global_get(job_name_prefix, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST\n"

                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME submitted\""
                    #add 0.2s sleep to let slurm submiting the job correctly
                    cmd += "\nsleep 0.1\n"

                    self.genpipes_file.write(cmd)
                    self.genpipes_file.flush()

        logger.info("\nGenpipes file generated\"")
        # Check cluster maximum job submission
        cluster_max_jobs = global_conf.global_get('DEFAULT', 'cluster_max_jobs', param_type='posint', required=False)
        if cluster_max_jobs and len(pipeline.jobs) > cluster_max_jobs:
            logger.warning(f"Number of jobs: {str(len(pipeline.jobs))} > Cluster maximum number of jobs: {str(cluster_max_jobs)} !")


class DaemonScheduler(Scheduler):

    def __init__(self, *args, **kwargs):
        super(DaemonScheduler, self).__init__(*args, **kwargs)
        self.name = 'DAEMON'

    def submit(self, pipeline):
        self.genpipes_file.write(self.json(pipeline))

    def json(self, pipeline):
        #with open('sample.json', 'w') as json_file:
        #json.dump(the_dump, json_file)
        return json.dumps(
            {'pipeline': {
                'output_dir': pipeline.output_dir,
                'samples': [{
                    'name': sample.name,
                    'readsets': [{
                        "name": readset.name,
                        "library": readset.library,
                        "runType": readset.run_type,
                        "run": readset.run,
                        "lane": readset.lane,
                        "adapter1": readset.adapter1,
                        "adapter2": readset.adapter2,
                        "qualityoffset": readset.quality_offset,
                        "bed": [bed for bed in readset.beds],
                        "fastq1": readset.fastq1,
                        "fastq2": readset.fastq2,
                        "bam": readset.bam,
                     } for readset in pipeline.readsets if readset.sample.name == sample.name]
                } for sample in pipeline.samples],
                'steps': [{
                    'name': step.name,
                    'jobs': [{
                        "job_name": job.name,
                        "job_id": job.id,
                        "job_command": job.command_with_modules,
                        "job_input_files": job.input_files,
                        "job_output_files": job.output_files,
                        "job_dependencies": [dependency_job.id for dependency_job in job.dependency_jobs],
                        "job_cluster_options": {
                            # Cluster settings section must match job name prefix before first "."
                            # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                            'cluster_submit_cmd': global_conf.global_get(job.name.split(".")[0], 'cluster_submit_cmd'),
                            'cluster_other_arg': global_conf.global_get(job.name.split(".")[0], 'cluster_other_arg'),
                            'cluster_work_dir_arg': global_conf.global_get(job.name.split(".")[0], 'cluster_work_dir_arg') + " " + pipeline.output_dir,
                            'cluster_output_dir_arg': global_conf.global_get(job.name.split(".")[0], 'cluster_output_dir_arg') + " " + os.path.join(pipeline.output_dir, "job_output", step.name, job.name + ".o"),
                            'cluster_job_name_arg': global_conf.global_get(job.name.split(".")[0], 'cluster_job_name_arg') + " " + job.name,
                            'cluster_walltime': global_conf.global_get(job.name.split(".")[0], 'cluster_walltime'),
                            'cluster_queue': global_conf.global_get(job.name.split(".")[0], 'cluster_queue'),
                            'cluster_cpu': global_conf.global_get(job.name.split(".")[0], 'cluster_cpu')
                        },
                        "job_done": job.done
                    } for job in step.jobs]
                } for step in pipeline.step_to_execute]
            }}, indent=4)
