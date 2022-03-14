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

# MUGQIC Modules
from core.config import *
from core.job import *

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
               threads = global_config_parser.param('hicup_align', 'threads'),
               Quiet = global_config_parser.param('hicup_align', 'Quiet'),
               Keep = global_config_parser.param('hicup_align', 'Keep'),
               Bowtie2_path = os.path.expandvars(global_config_parser.param('hicup_align', 'Bowtie2_path')),
               R_path = os.path.expandvars(global_config_parser.param('hicup_align', 'R_path')),
               Genome_Index_hicup = os.path.expandvars(global_config_parser.param('hicup_align', 'Genome_Index_hicup')),
               Genome_Digest = genome_digest,
               Format = global_config_parser.param('hicup_align', 'Format'),
               Longest = global_config_parser.param('hicup_align', 'Longest'),
               Shortest = global_config_parser.param('hicup_align', 'Shortest'),
               fastq1 = fastq1 + ".edited.gz",
               fastq2 = fastq2 + ".edited.gz")

    fileName = "hicup_align." + name + ".conf"

    command_confFile ='echo \"{configFileContent}\" > {fileName}'.format(fileName=fileName, configFileContent=configFileContent)

    return Job(input_files = [fastq1 + ".edited.gz", fastq2 + ".edited.gz"],
            name = "hicup_align.create_hicup_conf." + name,
            command = command_confFile,
            removable_files = [fileName]
            )


def hicup_run (name, confFile, sample_output_dir, fastq1, fastq2):

    command_hicup = "rm -rf {sample_output_dir} && \
                    mkdir -p {sample_output_dir} && \
                    hicup -c {fileName} && rm {fileName}".format(sample_output_dir = sample_output_dir, fileName = confFile)

    hicup_prefix = ".trim.pair1_2.fastq.gz.edited.hicup.bam"
    hicup_file_output = os.path.join(sample_output_dir, name + hicup_prefix)

    return Job(input_files = [fastq1 + ".edited.gz", fastq2 + ".edited.gz"],
            output_files = [hicup_file_output],
            module_entries = [['hicup_align', 'module_perl'], ['hicup_align', 'module_bowtie2'], ['hicup_align', 'module_samtools'], ['hicup_align', 'module_R'], ['hicup_align', 'module_mugqic_R_packages'], ['hicup_align', 'module_HiCUP']],
            name = "hicup_align." + name,
            command = command_hicup,
            )




