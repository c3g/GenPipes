#!/usr/bin/env python3

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

import argparse
import collections
import csv
import datetime
import logging
import os
import sys
import re
import subprocess

logger = logging.getLogger(__name__)


class JobStat(object):
    JOBID = 'JobId'
    JOBSTATE = 'JobState'
    SUBMITTIME = 'SubmitTime'
    STARTTIME = 'StartTime'
    ENDTIME = 'EndTime'
    NUMCPUS = 'NumCPUs'
    MEM = 'mem'
    REQUEUE = 'Requeue'
    RESTARTS = 'Restarts'
    RUNTIME = 'RunTime'
    REMOTE_LIST = {'narval': 'narval.calculcanada.ca',
              'beluga': 'beluga.calculcanada.ca',
              'cedar': 'cedar.calculcanada.ca',
              'graham': 'graham.calculcanada.ca'}

    # Fields extracted from sacct
    SJOBID = 'JobID'
    SSTATE = 'State'
    SELAPSEDRAW = 'ElapsedRaw'

    COMPLETED = 'COMPLETED'

    registery = {}
    _REMOTE = None

    def __init__(self, step_output_file, jobid=None, dependency=None, name=None,
                 parsed_after_header=50, remote_hpc=None):
        self.completed = None
        self._path = step_output_file
        if remote_hpc == "abacus":
            self.jobid = int(jobid.split(".")[0])
        else:
            self.jobid = int(jobid)
        self._dependency = dependency
        self._name = name
        self._prologue = {}
        self._epilogue = {}
        self._delta = None
        self._delta_with_overhead = None
        # Real slurm job number, if not none. There is a problem with
        # the pipeline output path
        self.output_log_id = []
        self._start = None
        self._stop = None
        self._sacct_start = None
        self._sacct_stop = None
        self._parsed_after_header = parsed_after_header
        self._mugqic_exit_status = None
        self.set_remote(remote_hpc)
        self._slurm_state = None
        self._sub_id = None
        self.log_missing = None
        self._elapsed_sec = None
        self._sacct_val = None
        if os.path.isfile(self._path):
            self.log_missing = False
            if remote_hpc == "abacus":
                self.fill_from_file_pbs(self._path)
            else:
                self.fill_from_file_slurm(self._path)
        else:
            self.log_missing = True
            logger.warning(f'{self._path} output path is missing')

        self.registery[self.jobid] = self
        # self.get_job_status()


    def __str__(self):
        return f'{self.jobid} {self.output_log_id} {self.name} {self.n_cpu} {self.mem} {self.runtime_job} {self.log_file_exit_status}'

    @classmethod
    def get_remote(cls):
        return cls._REMOTE

    @classmethod
    def set_remote(cls, remote_hpc):
        cls._REMOTE = remote_hpc

    @property
    def path(self):
        return os.path.realpath(self._path)

    @classmethod
    def sacct(cls, jobid):
        cmd = ["sacct", "--format=ALL", "-P", "--delimiter", "^", "-j", str(jobid)]
        if cls.get_remote() is not None:
            cmd = ['ssh', '-o', "StrictHostKeyChecking no", cls.REMOTE_LIST[cls.get_remote()]] + cmd
        return subprocess.check_output(cmd).decode("utf-8")

    @classmethod
    def set_all_status_slurm(cls):
        ids = ','.join([str(i) for i in cls.registery.keys()])
        raw_output = cls.sacct(ids)
        lines = raw_output.rstrip().splitlines()
        keys = lines[0].strip().split("^")
        job_i = [i for i, k in enumerate(keys) if k == cls.SJOBID][0]
        bidon = collections.defaultdict(list)
        for line in lines[1:]:
            slurm_steps = line.strip().split("^")
            job_id = slurm_steps[job_i].split('.')[0]
            bidon[job_id].append(slurm_steps)
        for job_id, steps in bidon.items():
            acct = dict(zip(keys, zip(*steps)))
            cls.registery[int(job_id)].set_job_status_slurm(acct)

    def set_job_status_slurm(self, sacct_val):
        self._sacct_val = sacct_val
        self._slurm_state = sacct_val[self.SSTATE]
        self._sub_id = sacct_val[self.SJOBID]

        if all([s == self.COMPLETED for s in self._slurm_state]):
            self.completed = True
        else:
            self.completed = False
            logger.error(f'The job {self.jobid} has a {self._slurm_state} slurm state ')

        self._elapsed_sec = [int(t) if t is not None else 0 for t in sacct_val[self.SELAPSEDRAW]]

    @property
    def slurm_state(self):
        return dict(zip(self._sub_id, self._slurm_state))

    def fill_from_file_slurm(self, log_file_path):
        # parse stdout and get all the prologue and epilogues values
        to_parse = False
        n = 0
        lines_to_parse = []
        if os.path.isfile(log_file_path):  # Make sure the file exists
            for line in open(log_file_path, encoding='utf-8', errors='ignore'):
                if "SLURM FAKE PROLOGUE" in line:
                    n = 0
                    to_parse = True
                elif "SLURM FAKE EPILOGUE" in line:
                    n = 0
                    to_parse = True
                elif "MUGQICexitStatus" in line:
                    self._mugqic_exit_status = line.strip('\n')

                if n >= self._parsed_after_header:
                    to_parse = False

                if to_parse:
                    lines_to_parse.append(line)
                    n += 1

            to_parse = ' '.join(lines_to_parse)
            bidon = re.findall(r"(\w+)=(\S+)\s", to_parse)
            all_value = {}
            for k, v in bidon:
                all_value.setdefault(k, []).append(v)

            try:
                fake_pro_epi = [i for i, x in enumerate(all_value[self.JOBID])
                                if self.jobid == int(x)]
            except KeyError:
                logger.warning(f'{log_file_path} has no jobID log')
                fake_pro_epi = []

            if len(fake_pro_epi) == 2:
                self.output_log_id = list(set(all_value[self.JOBID]))
                pro, epi = fake_pro_epi
            elif len(fake_pro_epi) >= 2:

                logger.error(f'{self.jobid} log is not in {self._path}, there is a mismatch '
                             'between the job number and the output')

                # Get the closest number
                diff = [abs(int(v) - self.jobid) for v in all_value[self.JOBID]]
                self.output_log_id = list(set(all_value[self.JOBID]))

                pro, epi = [i for i, x in enumerate(all_value[self.JOBID])
                            if int(all_value[self.JOBID][diff.index(min(diff))]) == int(x)]

            elif len(fake_pro_epi) == 1:
                # no epilogue!
                # Just put prologue in it
                pro = fake_pro_epi[0]
                epi = pro
                logger.error(f'No epilogue in log file {log_file_path}')
            else:
                pro = None
                epi = None

            for k, v in all_value.items():
                try:
                    self._prologue[k] = v[pro]
                    self._epilogue[k] = v[epi]
                except (IndexError, TypeError):
                    logger.warning(f'{k} = {v} is not a slurm prologue/epilogue value or is ambiguous')

            tres = re.findall(r"TRES=cpu=(\d+),mem=(\d+\.?\d*\w)", to_parse)
            try:
                self._prologue[self.NUMCPUS] = tres[pro][0]
                self._epilogue[self.NUMCPUS] = tres[epi][0]
                self._prologue[self.MEM] = tres[pro][1]
                self._epilogue[self.MEM] = tres[epi][1]
            except (IndexError, TypeError):
                logger.warning('epilogue and or prologue is not in the log')


    def fill_from_file_pbs(self, log_file_path):
        # parse stdout and get all the prologue and epilogues values
        pattern_start = re.compile(r"^Begin PBS Prologue")
        pattern_stop = re.compile(r"^End PBS Prologue")
        to_parse = False
        lines_to_parse = []
        if os.path.isfile(log_file_path):  # Make sure the file exists
            for line in open(log_file_path, encoding='utf-8', errors='ignore'):
                if pattern_start.match(line):
                    to_parse = True
                    pattern_start = re.compile(r"^Begin PBS Epilogue")
                if pattern_stop.search(line):
                    to_parse = False
                    pattern_stop = re.compile(r"^Begin PBS Epilogue")
                if to_parse:
                    lines_to_parse.append(line.rstrip('\n'))

            exit()
            to_parse = ' '.join(lines_to_parse)
            bidon = re.findall(r"(\w+)=(\S+)\s", to_parse)
            all_value = {}
            for k, v in bidon:
                all_value.setdefault(k, []).append(v)

            try:
                fake_pro_epi = [i for i, x in enumerate(all_value[self.JOBID])
                                if self.jobid == int(x)]
            except KeyError:
                logger.warning(f'{log_file_path} has no jobID log')
                fake_pro_epi = []

            if len(fake_pro_epi) == 2:
                self.output_log_id = list(set(all_value[self.JOBID]))
                pro, epi = fake_pro_epi
            elif len(fake_pro_epi) >= 2:

                logger.error(f'{self.jobid} log is not in {self._path}, there is a mismatch '
                             'between the job number and the output')

                # Get the closest number
                diff = [abs(int(v) - self.jobid) for v in all_value[self.JOBID]]
                self.output_log_id = list(set(all_value[self.JOBID]))

                pro, epi = [i for i, x in enumerate(all_value[self.JOBID])
                            if int(all_value[self.JOBID][diff.index(min(diff))]) == int(x)]

            elif len(fake_pro_epi) == 1:
                # no epilogue!
                # Just put prologue in it
                pro = fake_pro_epi[0]
                epi = pro
                logger.error(f'No epilogue in log file {log_file_path}')
            else:
                pro = None
                epi = None

            for k, v in all_value.items():
                try:
                    self._prologue[k] = v[pro]
                    self._epilogue[k] = v[epi]
                except (IndexError, TypeError):
                    logger.warning(f'{k} = {v} is not a slurm prologue/epilogue value or is ambiguous')

            tres = re.findall(r"TRES=cpu=(\d+),mem=(\d+\.?\d*\w)", to_parse)
            try:
                self._prologue[self.NUMCPUS] = tres[pro][0]
                self._epilogue[self.NUMCPUS] = tres[epi][0]
                self._prologue[self.MEM] = tres[pro][1]
                self._epilogue[self.MEM] = tres[epi][1]
            except (IndexError, TypeError):
                logger.warning('epilogue and or prologue is not in the log')

    @property
    def log_file_exit_status(self):
        try:
            return self._mugqic_exit_status.split(':')[1]
        except AttributeError:
            return None

    @property
    def runtime_job(self):
        try:
            if self._delta is None:
                self._delta = self.stop - self.start

            return self._delta
        except (KeyError, ValueError):
            return datetime.timedelta(0)

    @property
    def runtime_with_overhead(self):

        try:
            return datetime.timedelta(seconds=max(self._elapsed_sec))
        except ValueError:
            return datetime.timedelta(0)

    @property
    def start(self):
        if self._start is None:
            start = [int(s) for s in self._prologue[self.RUNTIME].split(':')]
            self._start = datetime.timedelta(hours=start[0],
                                             minutes=start[1], seconds=start[2])
        return self._start

    @property
    def stop(self):
        if self._stop is None:
            stop = [int(s) for s in self._epilogue[self.RUNTIME].split(':')]
            self._stop = datetime.timedelta(hours=stop[0], minutes=stop[1], seconds=stop[2])
        return self._stop

    @property
    def prologue_time(self):
        try:
            return self.start_time + self.start
        except (TypeError, KeyError, ValueError):
            return None

    @property
    def epilogue_time(self):
        try:
            return self.start_time + self.stop
        except (TypeError, KeyError, ValueError):
            return None

    @property
    def start_time(self):
        try:
            return datetime.datetime.strptime(self._prologue[self.STARTTIME], '%Y-%m-%dT%H:%M:%S')
        except KeyError:
            return None

    @property
    def end_time(self):
        try:
            return datetime.datetime.strptime(self._epilogue[self.ENDTIME], '%Y-%m-%dT%H:%M:%S')
        except KeyError:
            return None

    @property
    def n_cpu(self):
        try:
            return int(self._prologue[self.NUMCPUS])
        except KeyError:
            return 0

    @property
    def mem(self):
        try:
            return self._prologue[self.MEM]
        except KeyError:
            return None

    @property
    def name(self):
        return self._name


def get_report(job_list_tsv=None, remote_hpc=None):
    job_output_path = os.path.dirname(job_list_tsv)

    with open(job_list_tsv) as tsvin:
        jobs = csv.reader(tsvin, delimiter='\t')

        report = []
        for job in jobs:
            logger.info(f'loading {job}')
            report.append(
                JobStat(
                    step_output_file=os.path.join(job_output_path, job[3]),
                    name=job[1],
                    jobid=job[0],
                    dependency=job[2],
                    remote_hpc=remote_hpc
                    )
                )
        if remote_hpc == "abacus":
            pass
        else:
            JobStat.set_all_status_slurm()
    return report


def print_report(report, to_stdout=True, to_tsv=None):
    min_start = []
    max_stop = []
    total_machine = datetime.timedelta(0)
    total_machine_core = datetime.timedelta(0)
    ok_message = []
    bad_message = []

    header = (
        'id',
        'log_from_job',
        'same_id',
        'name',
        'slurm_prologue',
        'slurm_main',
        'slurm_epilogue',
        'custom_exit',
        'output_file_path',
        'outpout_file_exist',
        'runtime',
        'runtime_with_overhead',
        'start_time',
        'stop_time'
        )
    data_table = []

    for job_object in report:
        pro, batch, epi = None, None, None
        for k, v in job_object.slurm_state.items():
            if 'batch' in k:
                batch = v
            elif 'extern' in k:
                epi = v
            else:
                pro = v
        same = True
        if len(job_object.output_log_id) != 1:
            same = False
        elif job_object.jobid != int(job_object.output_log_id[0]):
            same = False

        data_table.append((
            job_object.jobid,
            job_object.output_log_id,
            same,
            job_object.name,
            pro,
            batch,
            epi,
            job_object.log_file_exit_status,
            job_object.path,
            not job_object.log_missing,
            job_object.runtime_job,
            job_object.runtime_with_overhead,
            job_object.start_time,
            job_object.stop_time,
            ))

        if to_stdout:
            if not job_object.completed:
                if job_object.log_missing:
                    bad_message.append(f"Job Name {job_object.name}, jobid {job_object.jobid} has no log file {job_object.path}. May not be started.")
                else:
                    bad_message.append(f"Job Name {job_object.name}, jobid {job_object.jobid} has not completed with status {job_object.slurm_state.keys()}, log file status {job_object.log_file_exit_status}")
            elif (job_object.output_log_id and job_object.jobid != int(job_object.output_log_id[0])) or len(job_object.output_log_id) > 1:
                bad_message.append(f"Job Name {job_object.name}, jobid {job_object.jobid} has log from job {job_object.output_log_id} in its log file {job_object.path}")
            else:
                ok_message.append(f"Job Name {job_object.name}, jobid {job_object.jobid} is all good!")

            if job_object.log_missing:
                total_machine = datetime.timedelta(0)
                total_machine_core = datetime.timedelta(0)
            else:
                if job_object.prologue_time is not None:
                    min_start.append(job_object.prologue_time)
                if job_object.epilogue_time is not None:
                    max_stop.append(job_object.epilogue_time)
                total_machine += job_object.runtime_job
                total_machine_core += job_object.runtime_job * job_object.n_cpu

    if to_stdout:
        try:
            min_start = min(min_start)
            max_stop = max(max_stop)
            # total on node
            total_human = max_stop - min_start
        except ValueError:
            logger.error('No time recorded on that job')
            total_human = None

        print(f"""
    Cumulative time spent on compute nodes: {total_machine}
    Cumulative core time: {total_machine_core}
    Human time from beginning of pipeline to its end: {total_human}"""
    )

    if to_tsv:
        if not isinstance(to_tsv, str):
            to_tsv = "./log_report.tsv"
        with open(to_tsv, 'w') as output_file:
            writer = csv.writer(output_file, delimiter='\t')
            writer.writerow(header)
            for row in data_table:
                writer.writerow(row)


def fine_grain(stats):
    total = datetime.timedelta(0)
    max_r = datetime.timedelta(0)
    max_h = datetime.timedelta(0)
    for job in stats:
        if "NA12878_PCRFRee_gatk4" in job.name:
            if 'realigner' in job.name:
                mar_r = max(max_r, job.runtime_job)
            if 'haplotype_caller' in job.name:
                mar_h = max(max_h, job.runtime_job)
            else:
                total = total + job.runtime_job

        elif job.jobid in [336, 338]:
            total = total + job.runtime_job

    total = total + max_h + max_r

    print(total.total_seconds() / 3600.)

    total = datetime.timedelta(0)
    max_r = datetime.timedelta(0)
    max_h = datetime.timedelta(0)
    for job in stats:
        if "NA12878_PCRFRee_20170509" in job.name:
            if 'realigner' in job.name:
                mar_r = max(max_r, job.runtime_job)
            if 'haplotype_caller' in job.name:
                mar_h = max(max_h, job.runtime_job)
            elif job.jobid in [336, 338]:
                pass
            else:
                total = total + job.runtime_job

    total = total + max_h + max_r

    print(total.total_seconds() / 3600.)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job_list_path', help="Path to a GenPipes job list")
    parser.add_argument('--remote', '-r', help="Remote HPC where the job was ran",
                        choices=['beluga', 'cedar', 'narval', 'abacus'],
                        default=None)
    parser.add_argument('--loglevel', help="Standard Python log level",
                        choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"],
                        default='ERROR')
    parser.add_argument('--tsv', help="output to tsv file")
    parser.add_argument('--quiet', '-q', help="No report printed to terminal",
                        action='store_true', default=False)

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)

    path = args.job_list_path
    remote = args.remote
    to_tsv = args.tsv
    to_stdout = not args.quiet

    stats = get_report(path, remote_hpc=remote)
    print_report(stats, to_tsv=to_tsv, to_stdout=to_stdout)
