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
    trim_log,
    trim_stats
    ):

    if input2:  # Paired end reads
        inputs = [input1, input2]
        outputs = [paired_output1, unpaired_output1, paired_output2, unpaired_output2]
    else:   # Single end reads
        inputs = [input1]
        outputs = [singleOutput]

    job = Job(inputs, outputs, [["trim", "moduleVersion.java"], ['trim', 'moduleVersion.trimmomatic']])

    # # Retrieve output directories removing duplicates if any
    output_dirs = list(collections.OrderedDict.fromkeys([os.path.dirname(output) for output in outputs]))
    threads = config.param('trim', 'nbThreads', type='int')
    adapter_file = config.param('trim', 'adapterFile', type='filepath')
    clip_settings = config.param('trim', 'clipSettings')
    min_quality = config.param('trim', 'minQuality', type='int')
    min_length = config.param('trim', 'minLength', type='int')
    headcrop = config.param('trim', 'headcrop', required=False, type='int')

    job.command = \
"""mkdir -p {output_dirs} && \\
java -XX:ParallelGCThreads=1 -Xmx2G -jar \$TRIMMOMATIC_JAR {mode} \\
  -threads {threads} \\
  -phred{quality_offset} \\
  {inputs} \\
  {outputs} \\
  ILLUMINACLIP:{adapter_file}{clip_settings} \\
  TRAILING:{min_quality} \\
  MINLEN:{min_length}""".format(
        output_dirs=" ".join(output_dirs),
        mode = "PE" if input2 else "SE",
        threads=threads,
        quality_offset=quality_offset if quality_offset == 64 else "33",
        inputs=" \\\n  ".join(inputs),
        outputs=" \\\n  ".join(outputs),
        adapter_file=adapter_file,
        clip_settings=clip_settings,
        min_quality=min_quality,
        min_length = min_length
    )

    if quality_offset == 64:
        job.command += " \\\n  TOPHRED33"

    if headcrop and headcrop > 0:
        job.command += " \\\n  HEADCROP:" + str(headcrop)

    job.command += " \\\n  2> " + trim_log

    # Compute statistics
    job.command += " && \\\ngrep ^Input " + trim_log + " | perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*$/Raw Fragments,\\1\\nFragment Surviving,\\2\\nSingle Surviving,\\3/' > " + trim_stats

    return job

def normalize(input1, output1):
    job = Job([input1], [output1])
    job.command = "normalize " + input1 + " > " + output1
    return job

def trinity(normalized_readsets):
    job = Job(normalized_readsets, ["Trinity.fasta", "Trinity_stats.csv"])
    job.command = "trinity " + " ".join(normalized_readsets) + " > Trinity.fasta"
    return job

def blastx(input1):
    job = Job([input1], ["blastx_nr.tsv"])
    job.command = "blastx " + input1 + " > blastx_nr.tsv"
    return job

def rsem(reference, fasta):
    job = Job([reference, fasta], ["rsem_" + fasta + ".fpkm"])
    job.command = "rsem -db " + reference + " -input " + fasta
    return job

def trinotate(input1):
    job = Job([input1], ["trinotate.tsv"])
    job.command = "trinotate " + input1 + " > trinotate.tsv"
    return job

def nozzle(input_files, output1):
    job = Job(input_files, [output1])
    job.command = "nozzle " + " ".join(input_files) + " > " + output1
    return job
