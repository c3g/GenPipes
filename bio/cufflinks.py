#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def cufflinks(input_bam, output_directory, gtf=None):

    job = Job(
        [input_bam],
        [os.path.join(output_directory, "transcripts.gtf")],
        [['cufflinks', 'moduleVersion.cufflinks']]
    )

    job.command = """\
mkdir -p {output_directory} && \\
cufflinks -q {other_options}{gtf} \\
  --max-bundle-frags {max_bundle_frags} \\
  --library-type {library_type} \\
  --output-dir {output_directory} \\
  --num-threads {num_threads} \\
  {input_bam}""".format(
        other_options=config.param('cufflinks', 'otherOptions', required=False),
        gtf=" \\\n  --GTF " + gtf if gtf else "",
        max_bundle_frags=config.param('cufflinks', 'cufflinksMaxFargs', type='int'),
        library_type=config.param('cufflinks', 'strandInfo'),
        output_directory=output_directory,
        num_threads=config.param('cufflinks', 'cufflinksThreads', type='int'),
        input_bam=input_bam
    )

    return job
