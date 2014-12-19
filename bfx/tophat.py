#!/usr/bin/env python

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

    gtf = config.param('tophat', 'gtf', required=False, type='filepath')
    transcriptome_bowtie_index = config.param('tophat', 'transcriptome_bowtie_index', required=False, type='prefixpath')

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
        other_options=config.param('tophat', 'other_options', required=False),
        gtf=" \\\n  --GTF " + gtf if gtf else "",
        transcriptome_index=" \\\n  --transcriptome-index " + transcriptome_bowtie_index if transcriptome_bowtie_index else "",
        rg_id=rg_id,
        rg_sample=rg_sample,
        rg_library=rg_library,
        rg_platform_unit=rg_platform_unit,
        rg_platform=rg_platform,
        rg_center=rg_center,
        library_type=config.param('DEFAULT', 'strand_info'),
        output_directory=output_directory,
        num_threads=config.param('tophat', 'threads'),
        bowtie_index=config.param('tophat', 'genome_bowtie_index', type='prefixpath'),
        reads1=reads1,
        reads2=" \\\n  " + reads2 if reads2 else ""
        )
    )
