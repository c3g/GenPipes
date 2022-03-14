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

# MUGQIC Modules
from core.config import *
from core.job import *


def nanopolish_ont(read_fastq_dir, run_name, sample_name, read_fast5_dir, seq_summary, output_dir, ini_section="artic_nanopolish"):
    """
    Run artic nanopolish pipeline for a sample

    :return: a job for an artic nanopolish run
    """
    if not os.path.isabs(read_fast5_dir):
        read_fast5_dir = os.path.join("..", read_fast5_dir)

    output_files = [
        os.path.join(output_dir, sample_name + ".consensus.fasta"),
        os.path.join(output_dir, sample_name + ".pass.vcf.gz"),
        os.path.join(output_dir, sample_name + ".pass.vcf.gz.tbi"),
        os.path.join(output_dir, sample_name + ".trimmed.rg.sorted.bam"),
        os.path.join(output_dir, sample_name + ".trimmed.rg.sorted.bam.bai"),
        os.path.join(output_dir, sample_name + ".primertrimmed.rg.sorted.bam"),
        os.path.join(output_dir, sample_name + ".primertrimmed.rg.sorted.bam.bai"),
        os.path.join(output_dir, run_name + "_" + sample_name + ".fastq")
        ]

    return Job(
        [read_fastq_dir],
        output_files,
        [[ini_section, "module_artic"]],
        command="""\
artic guppyplex \\
 --min-length {min_length} \\
 --max-length {max_length} \\
 --directory {read_fastq_dir} \\
 --prefix {run_name} && \\
artic minion \\
    --normalise {normalise} \\
    --threads {threads} \\
    --scheme-directory {primers_dir} \\
    --read-file {run_name}_{sample_name}.fastq \\
    --fast5-directory {read_fast5_dir} \\
    --sequencing-summary {seq_summary} \\
    {primers_version} \\
    {sample_name}""".format(
            min_length=global_config_parser.param(ini_section, 'min_length', required=True),
            max_length=global_config_parser.param(ini_section, 'max_length', required=True),
            read_fastq_dir=os.path.join("..", "..", read_fastq_dir),
            run_name=run_name,
            normalise=global_config_parser.param(ini_section, 'normalise', required=True),
            threads=global_config_parser.param(ini_section, 'threads', required=True),
            primers_dir=global_config_parser.param(ini_section, 'primers_dir', required=True),
            sample_name=sample_name,
            read_fast5_dir=os.path.join("..", "..", read_fast5_dir),
            seq_summary=os.path.join("..", "..", seq_summary),
            primers_version=global_config_parser.param(ini_section, 'primers_version', required=True)
        )
    )
