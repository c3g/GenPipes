#!/usr/bin/env python3

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
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import argparse
import csv
import datetime
import logging
import os
import re
from collections import OrderedDict

# Configure logging at the module level
logger = logging.getLogger(__name__)

class JobStat():
    """
    Class to store the job information
    """
    def __init__(self, output_file, job_id=None, dependency=None, job_name=None):
        self.output_file = output_file
        self.job_id = job_id
        self.dependency = dependency
        self.job_name = job_name
        self.user = None
        self.node = None
        self.priority = None
        self.status = None
        self.submit_time = None
        self.eligible_time = None
        self.start_time = None
        self.queue_time = None
        self.end_time = None
        self.total_time = None
        self.time_limit = None
        self.time_efficiency = None
        self.n_cpu = None
        self.cpu_time = None
        self.cpu_efficiency = None
        self.mem = None
        self.ave_mem = None
        self.max_mem = None
        self.mem_efficiency = None
        self.ave_diskr = None
        self.max_diskr = None
        self.ave_diskw = None
        self.max_diskw = None
        if os.path.isfile(self.output_file):
            self.fill_from_file(self.output_file)
        else:
            logger.warning(f'Job Name {self.job_name} with ID {self.job_id} has no log file {self.output_file}. Skipping this job in report...')
            self.output_file = None
            self.status = 'PENDING'

    def __repr__(self):
        return f'<JobStat {self.job_id} {self.job_name}>'

    @property
    def path(self):
        """
        Return the full path to the output file
        """
        if self.output_file:
            return os.path.realpath(self.output_file)

    def fill_from_file(self, log_file_path):
        """
        Fill the object with information from the PROLOGUE and EPILOGUE of the .o file
        Args:
            log_file_path: path to the output file
        """
        with open(log_file_path) as log_file:
            for line in log_file:
                # Prologue section
                if line.startswith('PROLOGUE - '):
                    user = re.search(r'User:\s+(.+$)', line)
                    if user:
                        self.user = user.group(1)
                    node = re.search(r'Node\(s\):\s+(.+$)', line)
                    if node:
                        self.node = node.group(1)
                    priority = re.search(r'Priority:\s+(.+$)', line)
                    if priority:
                        self.priority = priority.group(1)
                    self.status = "RUNNING"
                    submit_time = re.search(r'Submit Time:\s+(.+$)', line)
                    if submit_time:
                        self.submit_time = submit_time.group(1)
                    time_limit = re.search(r'Time Limit:\s+(.+$)', line)
                    if time_limit:
                        self.time_limit = time_limit.group(1)
                    n_cpu = re.search(r'Number of CPU\(s\) Requested:\s+(.+$)', line)
                    if n_cpu:
                        self.n_cpu = int(n_cpu.group(1))
                    mem = re.search(r'Memory Requested:\s+(.+$)', line)
                    if mem:
                        self.mem = mem.group(1)
                # Epilogue section
                if line.startswith('EPILOGUE - '):
                    status = re.search(r'Status:\s+(.+$)', line)
                    if status:
                        self.status = status.group(1)
                    eligible_time = re.search(r'Eligible Time:\s+(.+$)', line)
                    if eligible_time:
                        self.eligible_time = eligible_time.group(1)
                    start_time = re.search(r'Start Time:\s+(.+$)', line)
                    if start_time:
                        self.start_time = start_time.group(1)
                    queue_time = re.search(r'Time Spent in Queue:\s+(.+$)', line)
                    if queue_time:
                        self.queue_time = queue_time.group(1)
                    end_time = re.search(r'End Time:\s+(.+$)', line)
                    if end_time:
                        self.end_time = end_time.group(1)
                    total_time = re.search(r'Total Wall-clock Time:\s+(.+$)', line)
                    if total_time:
                        self.total_time = total_time.group(1)
                    time_efficiency = re.search(r'Time Efficiency \(% of Wall-clock Time to Time Limit\):\s+(.+$)', line)
                    if time_efficiency:
                        self.time_efficiency = time_efficiency.group(1)
                    cpu_time = re.search(r'Total CPU Time:\s+(.+$)', line)
                    if cpu_time:
                        self.cpu_time = cpu_time.group(1)
                    cpu_efficiency = re.search(r'CPU Efficiency \(% of CPU Time to Wall-clock Time\):\s+(.+$)', line)
                    if cpu_efficiency:
                        self.cpu_efficiency = cpu_efficiency.group(1)
                    ave_mem = re.search(r'Average Memory Usage:\s+(.+$)', line)
                    if ave_mem:
                        self.ave_mem = ave_mem.group(1)
                    max_mem = re.search(r'Maximum Memory Usage:\s+(.+$)', line)
                    if max_mem:
                        self.max_mem = max_mem.group(1)
                    mem_efficiency = re.search(r'Memory Efficiency \(% of Memory Requested to Maximum Memory Used\):\s+(.+$)', line)
                    if mem_efficiency:
                        self.mem_efficiency = mem_efficiency.group(1)
                    ave_diskr = re.search(r'Average Disk Read:\s+(.+$)', line)
                    if ave_diskr:
                        self.ave_diskr = ave_diskr.group(1)
                    max_diskr = re.search(r'Maximum Disk Read:\s+(.+$)', line)
                    if max_diskr:
                        self.max_diskr = max_diskr.group(1)
                    ave_diskw = re.search(r'Average Disk Write:\s+(.+$)', line)
                    if ave_diskw:
                        self.ave_diskw = ave_diskw.group(1)
                    max_diskw = re.search(r'Maximum Disk Write:\s+(.+$)', line)
                    if max_diskw:
                        self.max_diskw = max_diskw.group(1)

def get_report(job_list_tsv=None):
    """
    Get the job report from the job list
    Args:
        job_list_tsv: path to the job list tsv file
    Returns:
        list of JobStat objects
    """
    job_output_path = os.path.dirname(job_list_tsv)

    try:
        with open(job_list_tsv) as tsvin:
            jobs = csv.reader(tsvin, delimiter='\t')
            report = []
            for job in jobs:
                logger.info(f'Loading {job}')
                report.append(
                    JobStat(
                        output_file=os.path.join(job_output_path, job[3]),
                        job_name=job[1],
                        job_id=job[0],
                        dependency=job[2]
                        )
                    )
    except FileNotFoundError:
        logger.error(f'File {job_list_tsv} not found.')
        report = []
    return report

def parse_time(time_str):
    """
    Parse a time string in the format HH:MM:SS or MM:SS or SS
    Args:
        time_str: time string
    Returns:
        datetime.timedelta object
    """
    if time_str == 'Unknown':
        return None
    parts = time_str.split(':')
    if len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h = '0'
        m, s = parts
    elif len(parts) == 1:
        h = '0'
        m = '0'
        s = parts[0]
    else:
        raise ValueError(f"Unexpected time format: {time_str}")

    # Convert to float to handle decimal seconds
    h = int(h)
    m = int(m)
    s = float(s)

    # Separate seconds and milliseconds
    seconds = int(s)

    return datetime.timedelta(hours=h, minutes=m, seconds=seconds)


def print_report(report, to_stdout=True, to_tsv=None):
    """
    Print the report to the terminal and/or to a tsv file
    Args:
        report: list of JobStat objects
        to_stdout: print to terminal
        to_tsv: path to tsv file
    """
    if not report:
        logger.error('No job report to print.')
        return
    min_start = []
    max_stop = []
    total_machine = datetime.timedelta(0)
    total_machine_core = datetime.timedelta(0)
    status_count = {}

    header = OrderedDict([
        ('id', 'job_id'),
        ('name', 'job_name'),
        ('status', 'status'),
        ('user', 'user'),
        ('node', 'node'),
        ('priority', 'priority'),
        ('submit_time', 'submit_time'),
        ('eligible_time', 'eligible_time'),
        ('start_time', 'start_time'),
        ('queue_time', 'queue_time'),
        ('end_time', 'end_time'),
        ('total_time', 'total_time'),
        ('time_limit', 'time_limit'),
        ('time_efficiency', 'time_efficiency'),
        ('n_cpu', 'n_cpu'),
        ('cpu_time', 'cpu_time'),
        ('cpu_efficiency', 'cpu_efficiency'),
        ('mem', 'mem'),
        ('ave_mem', 'ave_mem'),
        ('max_mem', 'max_mem'),
        ('mem_efficiency', 'mem_efficiency'),
        ('ave_diskr', 'ave_diskr'),
        ('max_diskr', 'max_diskr'),
        ('ave_diskw', 'ave_diskw'),
        ('max_diskw', 'max_diskw'),
        ('output_file_path', 'path')
    ])
    data_table = []
    nb_jobs = len(report)
    for job in report:
        for dep in job.dependency.split(':'):
            if any(d.status not in ['COMPLETED', 'RUNNING', 'PENDING'] for d in report if d.job_id == dep):
                job.status = 'CANCELLED'
        line = []
        for _, value in header.items():
            line.append(getattr(job, value))
        data_table.append(line)
        if job.start_time:
            min_start.append(datetime.datetime.strptime(job.start_time, '%Y-%m-%dT%H:%M:%S'))
        else:
            logger.warning(f'Job Name {job.job_name} with ID {job.job_id} has no Start Time.')
        if job.end_time:
            max_stop.append(datetime.datetime.strptime(job.end_time, '%Y-%m-%dT%H:%M:%S'))
        else:
            logger.warning(f'Job Name {job.job_name} with ID {job.job_id} has no End Time.')
        if job.total_time:
            total_time = parse_time(job.total_time)
            if total_time:
                total_machine += total_time
        if job.cpu_time:
            cpu_time = parse_time(job.cpu_time)
            if cpu_time and job.n_cpu:
                total_machine_core += cpu_time * job.n_cpu
        # Count job status
        if job.status in status_count:
            status_count[job.status] += 1
        else:
            status_count[job.status] = 1


    if to_stdout:
        try:
            min_start = min(min_start)
            max_stop = max(max_stop)
            # total on node
            total_human = max_stop - min_start
        except ValueError:
            logger.error('Missing start or stop time in the report.')
            total_human = None

        status_counts_str = "\n".join([f"    Number of jobs {status}: {count} ({(count / nb_jobs) * 100:.2f}%)" for status, count in status_count.items()])
        print(f"""
{status_counts_str}
    Number of jobs: {nb_jobs}
    Cumulative time spent on compute nodes: {total_machine}
    Cumulative core time: {total_machine_core}
    Human time from beginning of pipeline to its end: {total_human}
    """)

    if to_tsv:
        if not isinstance(to_tsv, str):
            to_tsv = "./log_report.tsv"
        with open(to_tsv, 'w') as output_file:
            writer = csv.writer(output_file, delimiter='\t')
            writer.writerow(header.keys())
            for row in data_table:
                writer.writerow(row)

def main(args=None):
    """
    Main function to run the log_report script
    """
    if args is None:
        parser = argparse.ArgumentParser()
        parser.add_argument('job_list_path', help="Path to a GenPipes job list")
        parser.add_argument('--loglevel', help="Standard Python log level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', "CRITICAL"], default='WARNING')
        parser.add_argument('--tsv', help="output to tsv file")
        parser.add_argument('--quiet', '-q', help="No report printed to terminal", action='store_true', default=False)

        args = parser.parse_args()

    # Configure logging at the module level
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # Update logging level based on the argument
    logger.setLevel(getattr(logging, args.loglevel))

    path = args.job_list_path
    to_tsv = args.tsv
    to_stdout = not args.quiet

    stats = get_report(path)
    print_report(stats, to_tsv=to_tsv, to_stdout=to_stdout)


if __name__ == '__main__':
    main()
