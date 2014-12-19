#!/usr/bin/env python

# Python Standard Modules
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def trimmomatic(
    input1,
    input2,
    paired_output1,
    unpaired_output1,
    paired_output2,
    unpaired_output2,
    single_output,
    quality_offset,
    trim_log
    ):

    if input2:  # Paired end reads
        inputs = [input1, input2]
        outputs = [paired_output1, unpaired_output1, paired_output2, unpaired_output2]
    else:   # Single end reads
        inputs = [input1]
        outputs = [single_output]

    job = Job(inputs, outputs + [trim_log], [['trimmomatic', 'module_java'], ['trimmomatic', 'module_trimmomatic']], removable_files=[paired_output1, unpaired_output1, paired_output2, unpaired_output2, single_output])

    # Retrieve output directories removing duplicates if any
    output_dirs = list(collections.OrderedDict.fromkeys([os.path.dirname(output) for output in outputs]))

    threads = config.param('trimmomatic', 'threads', type='posint')
    adapter_file = config.param('trimmomatic', 'adapter_fasta', type='filepath')
    illumina_clip_settings = config.param('trimmomatic', 'illumina_clip_settings')
    trailing_min_quality = config.param('trimmomatic', 'trailing_min_quality', type='int')
    min_length = config.param('trimmomatic', 'min_length', type='posint')
    headcrop_length = config.param('trimmomatic', 'headcrop_length', required=False, type='posint')

    # CAUTION: Trimmomatic settings order is IMPORTANT! First, illuminaclip settings, then headcrop length, then trailing min quality, then minimum length
    job.command = """\
mkdir -p {output_dirs} && \\
java -XX:ParallelGCThreads=1 -Xmx{ram} -jar $TRIMMOMATIC_JAR {mode} \\
  -threads {threads} \\
  -phred{quality_offset} \\
  {inputs} \\
  {outputs} \\
  ILLUMINACLIP:{adapter_file}{illumina_clip_settings}{headcrop_length} \\
  TRAILING:{trailing_min_quality} \\
  MINLEN:{min_length}""".format(
        output_dirs=" ".join(output_dirs),
        ram=config.param('trimmomatic', 'ram'),
        mode = "PE" if input2 else "SE",
        threads=threads,
        quality_offset=quality_offset if quality_offset == 64 else "33",
        inputs=" \\\n  ".join(inputs),
        outputs=" \\\n  ".join(outputs),
        adapter_file=adapter_file,
        illumina_clip_settings=illumina_clip_settings,
        headcrop_length=" \\\n  HEADCROP:" + str(headcrop_length) if headcrop_length else "",
        trailing_min_quality=trailing_min_quality,
        min_length = min_length
    )

    if quality_offset == 64:
        job.command += " \\\n  TOPHRED33"

    job.command += " \\\n  2> " + trim_log

    return job
