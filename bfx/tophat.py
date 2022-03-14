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

def tophat(
    reads1,
    reads2,
    output_directory,
    rg_id="",
    rg_sample="",
    rg_library="",
    rg_platform_unit="",
    rg_platform="",
    rg_center=""
    ):

    gtf = global_config_parser.param('tophat', 'gtf', required=False, param_type='filepath')
    transcriptome_bowtie_index = global_config_parser.param('tophat', 'transcriptome_bowtie_index', required=False, param_type='prefixpath')

    return Job(
        [reads1, reads2],
        [os.path.join(output_directory, "accepted_hits.bam")],
        [
            ['tophat', 'module_bowtie'],
            ['tophat', 'module_samtools'],
            ['tophat', 'module_tophat']
        ],
        command="""\
mkdir -p {output_directory} && \\
tophat {other_options}{gtf}{transcriptome_index} \\
  --rg-id '{rg_id}' \\
  --rg-sample '{rg_sample}' \\
  --rg-library '{rg_library}' \\
  --rg-platform-unit '{rg_platform_unit}' \\
  --rg-platform '{rg_platform}' \\
  --rg-center '{rg_center}' \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {bowtie_index} \\
  {reads1}{reads2}""".format(
        other_options=global_config_parser.param('tophat', 'other_options', required=False),
        gtf=" \\\n  --GTF " + gtf if gtf else "",
        transcriptome_index=" \\\n  --transcriptome-index " + transcriptome_bowtie_index if transcriptome_bowtie_index else "",
        rg_id=rg_id,
        rg_sample=rg_sample,
        rg_library=rg_library,
        rg_platform_unit=rg_platform_unit,
        rg_platform=rg_platform,
        rg_center=rg_center,
        library_type=global_config_parser.param('DEFAULT', 'strand_info'),
        output_directory=output_directory,
        num_threads=global_config_parser.param('tophat', 'threads'),
        bowtie_index=global_config_parser.param('tophat', 'genome_bowtie_index', param_type='prefixpath'),
        reads1=reads1,
        reads2=" \\\n  " + reads2 if reads2 else ""
        )
    )
