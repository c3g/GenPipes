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

# Python Standard Modules
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def create_hicup_conf(name, fastq1, fastq2, sample_output_dir, genome_digest):
    # create HiCUP configuration file:
    configFileContent = """
    Outdir: {sample_output_dir}
    Threads: {threads}
    Quiet:{Quiet}
    Keep:{Keep}
    Zip:1
    Bowtie2: {Bowtie2_path}
    R: {R_path}
    Index: {Genome_Index_hicup}
    Digest: {Genome_Digest}
    Format: {Format}
    Longest: {Longest}
    Shortest: {Shortest}
    {fastq1}
    {fastq2}
    """.format(sample_output_dir = sample_output_dir,
               threads = global_conf.global_get('hicup_align', 'threads'),
               Quiet = global_conf.global_get('hicup_align', 'Quiet'),
               Keep = global_conf.global_get('hicup_align', 'Keep'),
               Bowtie2_path = os.path.expandvars(global_conf.global_get('hicup_align', 'Bowtie2_path')),
               R_path = os.path.expandvars(global_conf.global_get('hicup_align', 'R_path')),
               Genome_Index_hicup = os.path.expandvars(global_conf.global_get('hicup_align', 'Genome_Index_hicup')),
               Genome_Digest = genome_digest,
               Format = global_conf.global_get('hicup_align', 'Format'),
               Longest = global_conf.global_get('hicup_align', 'Longest'),
               Shortest = global_conf.global_get('hicup_align', 'Shortest'),
               fastq1 = fastq1 + ".edited.gz",
               fastq2 = fastq2 + ".edited.gz")

    fileName = "hicup_align." + name + ".conf"

    command_confFile ='echo \"{configFileContent}\" > {fileName}'.format(fileName=fileName, configFileContent=configFileContent)

    return Job(
                input_files=[fastq1 + ".edited.gz", fastq2 + ".edited.gz"],
                name="hicup_align.create_hicup_conf." + name,
                command=command_confFile,
                removable_files=[fileName]
                )


def hicup_run (
    name,
    sample_output_dir,
    fastq1,
    fastq2,
    genome_digest
    ):

    hicup_config = os.path.join(sample_output_dir, "hicup_align." + name + ".conf")

    command_hicup = """\\
rm -rf {sample_output_dir}/* && \\
Bowtie2_path=`which bowtie2`
R_path=`which R`
echo "Outdir: {sample_output_dir}
Threads: {threads}
Quiet:{Quiet}
Keep:{Keep}
Zip:1
Bowtie2: $Bowtie2_path
R: $R_path
Index: {Genome_Index_hicup}
Digest: {Genome_Digest}
Format: {Format}
Longest: {Longest}
Shortest: {Shortest}
{fastq1}
{fastq2}" > {hicup_config} && \\
hicup -c {hicup_config}""".format(
    sample_output_dir=sample_output_dir,
    threads=global_conf.global_get('hicup_align', 'threads'),
    Quiet=global_conf.global_get('hicup_align', 'Quiet'),
    Keep=global_conf.global_get('hicup_align', 'Keep'),
    Genome_Index_hicup=os.path.expandvars(global_conf.global_get('hicup_align', 'Genome_Index_hicup')),
    Genome_Digest=genome_digest,
    Format=global_conf.global_get('hicup_align', 'Format'),
    Longest=global_conf.global_get('hicup_align', 'Longest'),
    Shortest=global_conf.global_get('hicup_align', 'Shortest'),
    fastq1=fastq1 + ".edited.gz",
    fastq2=fastq2 + ".edited.gz",
    hicup_config=hicup_config
    )

    hicup_prefix = ".trim.pair1_2.fastq.gz.edited.hicup.bam"
    hicup_file_output = os.path.join(sample_output_dir, name + hicup_prefix)

    return Job(
                input_files=[fastq1 + ".edited.gz", fastq2 + ".edited.gz"],
                output_files=[hicup_file_output],
                module_entries=[
                    ['hicup_align', 'module_perl'],
                    ['hicup_align', 'module_bowtie2'],
                    ['hicup_align', 'module_samtools'],
                    ['hicup_align', 'module_R'],
                    ['hicup_align', 'module_mugqic_R_packages'],
                    ['hicup_align', 'module_HiCUP']
                    ],
                name="hicup_align." + name,
                command=command_hicup,
                removable_files=[hicup_config]
                )
