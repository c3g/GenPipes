#!/usr/bin/env python

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
import csv
import logging
import subprocess
import shutil

logger = logging.getLogger(__name__)

def delete_all_jobs(job_list_file):

    jobs = csv.reader(job_list_file, delimiter='\t')

    job_id = [job[0] for job in jobs]
    cmd = []
    if shutil.which('scancel'):
        cmd = ["scancel"] + job_id
    elif shutil.which('qdel'):
        cmd = ["qdel"] + job_id
    else:
        raise OSError('We only support PBS/Torque or SLRUM')

    logger.info(' '.join(cmd))
    try:
        raw_output = subprocess.check_output(cmd).decode("utf-8")
        logger.info('Done')
        logger.info(raw_output)
    except subprocess.CalledProcessError:
        logger.info('Some of the steps are already completed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help='will deleted all jobs submit by Genpipes on a SLURM or PBS like HPC')

    parser.add_argument('job_list_path', help="Path to a GenPipes job list", type=argparse.FileType('r'))

    parser.add_argument('--loglevel', help="Standard Python log level",
                        choices=['ERROR', 'WARNING', 'INFO', "CRITICAL"],
                        default='INFO')

    args = parser.parse_args()

    log_level = args.loglevel
    logging.basicConfig(level=log_level)
    job_list_file = args.job_list_path
    delete_all_jobs(job_list_file=job_list_file)