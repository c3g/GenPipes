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

# Python Standard Modules
import os
from packaging import version

# MUGQIC Modules
from core.config import *
from core.job import *


def guppy_basecalling(fast5_directory,
                      output_directory,
                      transfer_mode=True):
    """
    Use guppy to do basecalling for the output of an ONT sequencer.

    REQUIRES: access to GPU cores
    :return: a job for guppy basecalling
    """

    # Input for Guppy is a directory containing all the raw FAST5 files stored in fast5_directory

    # Output for Guppy is a directory containing the pass/fail FASTQ files (non-demultiplex)
    # As well as a sequencing_summary.txt file that is important for other downstream applications
    read_fastq_dir = output_directory
    sequencing_summary = os.path.join(output_directory, "sequencing_summary.txt")

    if version.parse(global_config_parser.param('guppy_basecall', 'module_guppy').split("gpu-")[1]) < version.parse("4.5.2"):
        qscore_flag = "--qscore_filtering "
    else:
        qscore_flag = ""

    if transfer_mode:
        fast5_directory_tmp = os.path.join(global_config_parser.param('guppy_basecall', 'tmp_dir', required=True), os.path.basename("".join(fast5_directory)))
        transfer_job = "transfer.sh {threads} {fast5_directory} {fast5_directory_tmp} && \\".format(
            threads=global_config_parser.param('guppy_basecall', 'transfer_threads') if global_config_parser.param('guppy_basecall', 'transfer_threads') else "1",
            fast5_directory="".join(fast5_directory),
            fast5_directory_tmp=global_config_parser.param('guppy_basecall', 'tmp_dir', required=True)
            )
    else:
        transfer_job = "\\"
        # fast5_directory = "".join(fast5_directory)

    job = Job(
        fast5_directory,
        [read_fastq_dir, os.path.join(output_directory, "pass"), sequencing_summary],
        [["guppy_basecall", "module_guppy"],
        ["guppy_basecall", "module_mugqic_tools"]],
        command="""\
{transfer_job}
guppy_basecaller --input_path {fast5_directory} \\
    --device auto \\
    --save_path {read_fastq_dir} \\
    --config {basecall_protocol} \\
    {qscore_flag}--min_qscore {min_Q_score} \\
    --records_per_fastq 0 \\
    --compress_fastq \\
    --disable_pings \\
    {other_options}""".format(
            transfer_job=transfer_job,
            fast5_directory=fast5_directory_tmp if transfer_mode else "".join(fast5_directory),
            read_fastq_dir=read_fastq_dir,
            basecall_protocol=global_config_parser.param('guppy_basecall', 'basecall_protocol', required=True),
            qscore_flag=qscore_flag,
            min_Q_score=global_config_parser.param('guppy_basecall', 'min_Q_score', required=True),
            other_options=global_config_parser.param('guppy_basecall', 'other_options', required=False),

        )
    )

    return job

def guppy_demultiplex(fastq_directory,
                      sequencing_summary,
                      output_directory,
                      barcodes_directories,
                      transfer_mode=True):
    """
    Use guppy to do demultiplexing for the output of an ONT sequencer.

    REQUIRES: access to GPU cores
    :return: a job for guppy demultiplexing
    """
    # Input for Guppy Demultiplex is a directory containing all the fastq files stored in fastq_directory

    # Output for Guppy is a directory containing the pass/fail FASTQ files (non-demultiplex)
    # As well as a sequencing_summary.txt file that is important for other downstream applications
    demultiplex_dir = output_directory
    barcoding_summary = os.path.join(output_directory, "barcoding_summary.txt")

    if transfer_mode:
        fastq_directory_tmp = os.path.join(global_config_parser.param('guppy_demultiplex', 'tmp_dir', required=True), os.path.basename(fastq_directory))
        transfer_job = "transfer.sh {threads} {fastq_directory} {fastq_directory_tmp} && \\".format(
            threads=global_config_parser.param('guppy_demultiplex', 'transfer_threads') if global_config_parser.param('guppy_demultiplex', 'transfer_threads') else "1",
            fastq_directory=fastq_directory,
            fastq_directory_tmp=global_config_parser.param('guppy_demultiplex', 'tmp_dir', required=True)
            )
    else:
        transfer_job = "\\"

    job = Job(
        [fastq_directory, sequencing_summary],
        [demultiplex_dir, barcoding_summary] + barcodes_directories,
        [["guppy_demultiplex", "module_guppy"],
        ["guppy_demultiplex", "module_mugqic_tools"]],
        command="""\
{transfer_job}
guppy_barcoder --input_path {fastq_directory} \\
    --device auto \\
    --save_path {demultiplex_dir} \\
    --arrangements_files {arrangements_files} \\
    --records_per_fastq 0 \\
    --compress_fastq \\
    --disable_pings \\
    {other_options}""".format(
            transfer_job=transfer_job,
            fastq_directory=fastq_directory_tmp if transfer_mode else fastq_directory,
            demultiplex_dir=demultiplex_dir,
            arrangements_files=global_config_parser.param('guppy_demultiplex', 'arrangements_files', required=True),
            other_options=global_config_parser.param('guppy_demultiplex', 'other_options', required=False)
        )
    )

    return job
