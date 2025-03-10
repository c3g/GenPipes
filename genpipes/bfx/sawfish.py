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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def discover(
    input_file,
    output_dir,
    ini_section='sawfish'
    ):
    """
    Discover structural variants with sawfish.

    :return: a job for sawfish discover
    """

    genome_fasta = global_conf.global_get(ini_section, 'genome_fasta', required=True)
    outputs = [
        output_dir,
        output_dir + "/contig.alignment.bam",
        output_dir + "/sawfish.log"
    ]

    return Job(
        [input_file],
        outputs,
        [
            [ini_section, "module_sawfish"]
        ],
        command="""\
sawfish discover {other_options} \\
  --ref {genome_fasta} \\
  {expected_cn} \\
  --bam {input_file} \\
  --threads {threads} \\
  --output-dir {output_dir}""".format(
            threads=global_conf.global_get(ini_section, 'threads'),
            other_options=global_conf.global_get(ini_section, 'discover_options', required=False),
            genome_fasta=genome_fasta,
            expected_cn="--expected-cn " + global_conf.global_get(ini_section, 'cn_bed', required=False) if global_conf.global_get(ini_section, 'cn_bed', required=False) else "",
            input_file=input_file,
            output_dir=output_dir
        ),
    )

def joint_call(
    input_dir,
    output_dir,
    ini_section='sawfish'
    ):
    """
    Call structural variants with sawfish.

    :return: a job for sawfish discover
    """

    outputs = [
        output_dir + "/genotyped.sv.vcf.gz"
    ]

    return Job(
        [input_dir],
        outputs,
        [
            [ini_section, "module_sawfish"]
        ],
        command="""\
sawfish joint-call {other_options} \\
  --sample {input_dir} \\
  --threads {threads} \\
  --output-dir {output_dir}""".format(
            threads=global_conf.global_get(ini_section, 'threads'),
            other_options=global_conf.global_get(ini_section, 'call_options', required=False),
            input_dir=input_dir,
            output_dir=output_dir
        ),
    )
