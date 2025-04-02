################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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

import logging
import os
import re
import stat
import sys
import tempfile

from .utils import expandvars, time_to_datetime

from .config import global_conf

# Output comment separator line
SEPARATOR_LINE = "#" + "-" * 79

log = logging.getLogger(__name__)

def create_scheduler(s_type, config_files, container=None, genpipes_file=None):
    """
    Create a scheduler object
    Args:
        s_type (str): Scheduler type
        config_files (list): List of configuration files
        container (Container): Container object
        genpipes_file (file): GenPipes file
    Returns:
        Scheduler: Scheduler object
    """
    if s_type == "pbs":
        return PBSScheduler(config_files, container=container, genpipes_file=genpipes_file)
    if s_type == "batch":
        return BatchScheduler(config_files, container=container, genpipes_file=genpipes_file)
    if s_type == "slurm":
        return SlurmScheduler(config_files, container=container, genpipes_file=genpipes_file)
    raise Exception(f"""Error: scheduler type "{s_type}" is invalid!""")


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
            raise NotImplementedError(f'submit_cmd needs to be implemented for {self.__class__} class')
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
            self._host_cvmfs_cache = global_conf.global_get("container", 'host_cvmfs_cache', required=False, param_type="string")
            if not self._host_cvmfs_cache:
                tmp_dir = global_conf.global_get("DEFAULT", 'tmp_dir', required=True)
                tmp_dir = expandvars(tmp_dir)
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
        """
        Return the container line
        Returns:
            str: Container line
        """
        if self.container:
            if self.container.type == 'docker':
                v_opt = f' -v {self.host_cvmfs_cache}:{self.cvmfs_cache}'
                for b in self.bind:
                    v_opt += f' -v {b}:{b}'
                network = " --network host"
                user = " --user $UID:$GROUPS"
                return f'docker run --env-file <( env| cut -f1 -d= ) --privileged {network} {user} {v_opt} {self.container.name} '

            if self.container.type == 'singularity':
                b_opt = f' -B {self.host_cvmfs_cache}:{self.cvmfs_cache}'
                for b in self.bind:
                    b_opt += f' -B {b}:{b}'
                return f"singularity run {b_opt} {self.container.name} "

            if self.container.type == 'wrapper':
                return f"{self.container.name} "
        else:
            return ""


    def fail_on_pattern(self, job_name_prefix):
        """
        Return the fail_on_pattern command and test_condition
        Args:
            job_name_prefix (str): Job name prefix
        Returns:
            tuple: fail_on_pattern command and test_condition
        """
        pattern = global_conf.global_get(job_name_prefix, 'fail_on_pattern',required=False)

        if not pattern:
            return ("", "")

        tmp_dir = global_conf.global_get("DEFAULT", 'tmp_dir', required=True)
        append_command = f" | tee {tmp_dir}/${{JOB_NAME}}_${{TIMESTAMP}}.o "
        test_condition = f"""
grep {pattern} {tmp_dir}/${{JOB_NAME}}_${{TIMESTAMP}}.o
NO_PROBLEM_IN_LOG=\\$?

  if [[  \\$NO_PROBLEM_IN_LOG == 0 ]] ; then
   echo {pattern} found in job log, forcing error
   GenPipes_STATE=74
fi
"""
        return (append_command, test_condition)



    def print_header(self, pipeline, shebang='/bin/bash'):
        """
        Print the header of the GenPipes file
        Args:
            pipeline (Pipeline): Pipeline object
            shebang (str): Shebang to use in the script
        """

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
                separator_line=SEPARATOR_LINE,
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
        if pipeline.jobs:
            self.genpipes_file.write("cd $OUTPUT_DIR")


    def print_step(self, step):
        """
        Print the step in the GenPipes file
        Args:
            step (Step): Step object
        """
        self.genpipes_file.write(f"""
{SEPARATOR_LINE}
# STEP: {step.name}
{SEPARATOR_LINE}
STEP={step.name}
mkdir -p $JOB_OUTPUT_DIR/$STEP
"""
                                 )

    def job2json_project_tracking(self, pipeline, job, job_status):
        """
        Generate the job2json_project_tracking command for the job
        Args:
            pipeline (Pipeline): Pipeline
            job (Job): Job
            job_status (str): Job status
        Returns:
            str: job2json_project_tracking command
        """
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
            job2json_project_tracking_script="genpipes tools job2json_project_tracking",
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
        time = time_to_datetime(walltime)
        sec = int(time.seconds % 60)
        minutes = int(((time.seconds - sec) / 60) % 60)
        hours = int((time.seconds - sec - 60 * minutes) / 3600 + time.days * 24)
        return f'-l walltime={hours:02d}:{minutes:02d}:{sec:02d}'

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

        # Checking if this is CIT and if the cpu is 1, then we will increase it to 2 to avoid a problem with cpulimit exceeding cpu usage
        continuous_integration_testing = 'GENPIPES_CIT' in os.environ
        if continuous_integration_testing:
            cpu += 1

        node = self.node(job_name_prefix)
        gpu = self.gpu(job_name_prefix)
        if gpu:
            return f"-l nodes={node}:ppn={cpu}:gpus={gpu}"
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
                        job_dependencies = f"""JOB_DEPENDENCIES={":".join([f"${dependency_job.id}" for dependency_job in dependency_chunks[0]])}"""
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += f"""\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES:{":".join([f"${dependency_job.id}" for dependency_job in dependency_chunk])}"""
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    job_name_prefix = job.name.split(".")[0]
                    log.debug(f"For {job.name} the ini section considered is [{job_name_prefix}]")
                    config_step_wrapper = global_conf.global_get(job_name_prefix, 'step_wrapper', required=False)

                    self.genpipes_file.write(f"""\
{SEPARATOR_LINE}
# JOB: {job.id}: {job.name}
{SEPARATOR_LINE}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
SCIENTIFIC_FILE=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.sh
SUBMISSION_FILE=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.submit
# Create the scientific file {os.path.basename(job.done)} is used as EOF marker to avoid issues with special characters like $ in the command
cat << '{os.path.basename(job.done)}' > $SCIENTIFIC_FILE
{job.command_with_modules}
{os.path.basename(job.done)}
chmod 755 $SCIENTIFIC_FILE
"""
                    )

                    self.genpipes_file.flush()

                    # job_name_prefix = job.name.split(".")[0]
                    if job.dependency_jobs:
                        dependencies = f"#PBS {self.dependency_arg(job_name_prefix)}$JOB_DEPENDENCIES"
                    else:
                        dependencies = ""
                    memory = self.memory(job_name_prefix, adapt=pipeline.force_mem_per_cpu)
                    if memory:
                        memory = f"#PBS {memory}"
                    else:
                        memory = ""
                    # Epilogue script have to be updated by IT for Abacus because having it anywhere else than HOME fails
                    # /!\ Make sure anytime epilogue script is updated, IT is informed so they update the script on Abacus/!\
                    cmd = f"""\
# Create the submission file
echo "#!/bin/bash
#PBS -T GenPipes
#PBS {global_conf.global_get(job_name_prefix, 'cluster_other_arg')} {global_conf.global_get(job_name_prefix, 'cluster_queue')}
#PBS -d $OUTPUT_DIR
#PBS -j oe
#PBS -o $JOB_OUTPUT
#PBS -N $JOB_NAME
#PBS {self.walltime(job_name_prefix)}
#PBS {self.cpu(job_name_prefix, adapt=pipeline.force_mem_per_cpu)}
{memory}
{dependencies}
{os.path.dirname(os.path.abspath(__file__))}/prologue.py
{self.job2json_project_tracking(pipeline, job, "RUNNING")}
{config_step_wrapper} {self.container_line} bash $SCIENTIFIC_FILE
GenPipes_STATE=\\$PIPESTATUS
echo GenPipesExitStatus:\\$GenPipes_STATE
{self.job2json_project_tracking(pipeline, job, '\\$GenPipes_STATE')}
if [ \\$GenPipes_STATE -eq 0 ]; then
    touch $JOB_DONE
fi
exit \\$GenPipes_STATE" > $SUBMISSION_FILE
# Submit the job and get the job id
{job.id}=$({self.submit_cmd} $SUBMISSION_FILE | awk '{{print $1}}')
# Write job parameters in job list file
echo "${job.id}\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST
echo "Submitted job with ID: ${job.id}"
# sleep to let the scheduler submiting the job correctly
sleep 0.1
"""
                    # Cluster settings section must match job name prefix before first "."
                    # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"

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

    def submit(self, pipeline):
        log.info(f'\n\t To run the script use: \n\t"{self.container_line}  ./<command>.sh"')
        self.print_header(pipeline)
        if pipeline.jobs:
            self.genpipes_file.write("\nSEPARATOR_LINE=`seq -s - 80 | sed 's/[0-9]//g'`\n")
        for step in pipeline.step_to_execute:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    job_name_prefix = job.name.split(".")[0]
                    log.debug(f"For {job.name} the ini section considered is [{job_name_prefix}]")
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
rm -f $JOB_DONE && {job2json_project_tracking_start} {step_wrapper} $COMMAND &> $JOB_OUTPUT
GenPipes_STATE=$?
echo "End GenPipes Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo GenPipesExitStatus:$GenPipes_STATE
{job2json_project_tracking_end}
if [ $GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $GenPipes_STATE ; fi
""".format(
                            job=job,
                            limit_string=os.path.basename(job.done),
                            separator_line=SEPARATOR_LINE,
                            job2json_project_tracking_start=self.job2json_project_tracking(pipeline, job, '\\"RUNNING\\"'),
                            job2json_project_tracking_end=self.job2json_project_tracking(pipeline, job, '\\$GenPipes_STATE'),
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
        time = time_to_datetime(walltime)
        sec = int(time.seconds % 60)
        minutes = int(((time.seconds - sec) / 60) % 60)
        hours = int((time.seconds - sec - 60 * minutes) / 3600 + time.days * 24)
        return f'--time={hours:02d}:{minutes:02d}:{sec:02d}'

    def gpu(self, job_name_prefix):
        n_gpu = super().gpu(job_name_prefix)
        gpu_type = self.gpu_type(job_name_prefix)
        if gpu_type and n_gpu:
            return f'--gres=gpu:{gpu_type}:{n_gpu}'
        if n_gpu:
            return f'--gres=gpu:{n_gpu}'
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
                        job_dependencies = f"""JOB_DEPENDENCIES={":".join([f"${dependency_job.id}" for dependency_job in dependency_chunks[0]])}"""
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += f"""\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES:{":".join([f"${dependency_job.id}" for dependency_job in dependency_chunk])}"""
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    job_name_prefix = job.name.split(".")[0]
                    log.debug(f"For {job.name} the ini section considered is [{job_name_prefix}]")
                    config_step_wrapper = global_conf.global_get(job_name_prefix, 'step_wrapper', required=False)

                    self.genpipes_file.write(f"""\
{SEPARATOR_LINE}
# JOB: {job.id}: {job.name}
{SEPARATOR_LINE}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
SCIENTIFIC_FILE=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.sh
SUBMISSION_FILE=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.submit
# Create the scientific file {os.path.basename(job.done)} is used as EOF marker to avoid issues with special characters like $ in the command
cat << '{os.path.basename(job.done)}' > $SCIENTIFIC_FILE
{job.command_with_modules}
{os.path.basename(job.done)}
chmod 755 $SCIENTIFIC_FILE
"""
                    )

                    self.genpipes_file.flush()

                    # job_name_prefix = job.name.split(".")[0]
                    if job.dependency_jobs:
                        dependencies = f"#SBATCH {self.dependency_arg(job_name_prefix)}$JOB_DEPENDENCIES"
                    else:
                        dependencies = ""

                    cmd = f"""\
# Create the submission file
echo "#!/bin/bash
#SBATCH {global_conf.global_get(job_name_prefix, 'cluster_other_arg')} {global_conf.global_get(job_name_prefix, 'cluster_queue')}
#SBATCH -D $OUTPUT_DIR
#SBATCH -o $JOB_OUTPUT
#SBATCH -J $JOB_NAME
#SBATCH {self.walltime(job_name_prefix)}
#SBATCH {self.memory(job_name_prefix)}
#SBATCH {self.cpu(job_name_prefix)} {self.gpu(job_name_prefix)}
{dependencies}
EPILOGUE_SCRIPT="{os.path.dirname(os.path.abspath(__file__))}/epilogue.py"
trap "\\$EPILOGUE_SCRIPT" EXIT
{os.path.dirname(os.path.abspath(__file__))}/prologue.py
{self.job2json_project_tracking(pipeline, job, "RUNNING")}
srun --wait=0 {config_step_wrapper} {self.container_line} bash $SCIENTIFIC_FILE
GenPipes_STATE=\\$PIPESTATUS
echo GenPipesExitStatus:\\$GenPipes_STATE
{self.job2json_project_tracking(pipeline, job, '\\$GenPipes_STATE')}
if [ \\$GenPipes_STATE -eq 0 ]; then
    touch $JOB_DONE
fi
exit \\$GenPipes_STATE" > $SUBMISSION_FILE
# Submit the job and get the job id
{job.id}=$({self.submit_cmd} $SUBMISSION_FILE | awk '{{print $4}}')
# Write job parameters in job list file
echo "${job.id}\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST
echo "Submitted job with ID: ${job.id}"
# sleep to let the scheduler submiting the job correctly
sleep 0.1
"""
                    # Cluster settings section must match job name prefix before first "."
                    # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"

                    self.genpipes_file.write(cmd)
                    self.genpipes_file.flush()

        log.info("\nGenpipes file generated\"")
        # Check cluster maximum job submission
        cluster_max_jobs = global_conf.global_get('DEFAULT', 'cluster_max_jobs', param_type='posint', required=False)
        if cluster_max_jobs and len(pipeline.jobs) > cluster_max_jobs:
            log.warning(f"Number of jobs: {str(len(pipeline.jobs))} > Cluster maximum number of jobs: {str(cluster_max_jobs)} !")
