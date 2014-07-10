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

    job = Job(
        [reads1, reads2],
        [os.path.join(output_directory, "accepted_hits.bam")],
        [['tophat', 'moduleVersion.bowtie'], ['tophat', 'moduleVersion.samtools'], ['tophat', 'moduleVersion.tophat']]
    )

    gtf = config.param('tophat', 'referenceGtf', required=False, type='filepath')
    transcriptome_bowtie_index = config.param('tophat', 'transcriptomeBowtieIndex', required=False, type='prefixpath')

    job.command = """\
tophat {other_options}{gtf}{transcriptome_index} \\
  --rg-id '{rg_id}' \\
  --rg-sample '{rg_sample}' \\
  --rg-library '{rg_library}' \\
  --rg-platform_unit '{rg_platform_unit}' \\
  --rg-platform '{rg_platform}' \\
  --rg-center '{rg_center}' \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {bowtie_index} \\
  {reads1}{reads2}""".format(
        other_options=config.param('tophat', 'otherOptions', required=False),
        gtf=" \\\n  --GTF " + gtf if gtf else "",
        transcriptome_index=" \\\n  --transcriptome-index " + transcriptome_bowtie_index if transcriptome_bowtie_index else "",
        rg_id=rg_id,
        rg_sample=rg_sample,
        rg_library=rg_library,
        rg_platform_unit=rg_platform_unit,
        rg_platform=rg_platform,
        rg_center=rg_center,
        library_type=config.param('tophat', 'strandInfo'),
        output_directory=output_directory,
        num_threads=config.param('tophat', 'TBAlnThreads'),
        bowtie_index=config.param('tophat', 'bowtie_index', type='prefixpath'),
        reads1=reads1,
        reads2=" \\\n  " + reads2 if reads2 else ""
    )

    return job
