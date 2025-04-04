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
import logging

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

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
    adapter_file,
    trim_log
    ):

    if input2:  # Paired end reads
        inputs = [input1, input2]
        outputs = [paired_output1, unpaired_output1, paired_output2, unpaired_output2]
    else:   # Single end reads
        inputs = [input1]
        outputs = [single_output]

    headcrop_length = global_conf.global_get('trimmomatic', 'headcrop_length', required=False, param_type='posint')

    return Job(
        inputs,
        outputs + [trim_log],
        [
            ['trimmomatic', 'module_java'],
            ['trimmomatic', 'module_trimmomatic']
        ],

        # CAUTION: Trimmomatic settings order is IMPORTANT!
        # FIRST Illuminaclip settings, THEN headcrop length, THEN trailing min quality, THEN minimum length
        command="""\
java {java_other_options} -Xmx{ram} -jar $TRIMMOMATIC_JAR {mode} \\
  -threads {threads} \\
  -phred{quality_offset} \\
  {inputs} \\
  {outputs} \\
  ILLUMINACLIP:{adapter_file}{illumina_clip_settings}{headcrop_length} \\
  TRAILING:{trailing_min_quality} \\
  MINLEN:{min_length}{tophred33} \\
  2>&1 | tee -a {trim_log}""".format(
        java_other_options=global_conf.global_get('trimmomatic', 'java_other_options'),
        ram=global_conf.global_get('trimmomatic', 'ram'),
        mode = "PE" if input2 else "SE",
        threads=global_conf.global_get('trimmomatic', 'threads', param_type='posint'),
        quality_offset=quality_offset if quality_offset == 64 else "33",
        inputs=" \\\n  ".join(inputs),
        outputs=" \\\n  ".join(outputs),
        adapter_file=adapter_file,
        illumina_clip_settings=global_conf.global_get('trimmomatic', 'illumina_clip_settings'),
        headcrop_length=" \\\n  HEADCROP:" + str(headcrop_length) if headcrop_length else "",
        trailing_min_quality=global_conf.global_get('trimmomatic', 'trailing_min_quality', param_type='int'),
        min_length=global_conf.global_get('trimmomatic', 'min_length', param_type='posint'),
        tophred33=" \\\n  TOPHRED33" if quality_offset == 64 else "",
        trim_log=trim_log
        ),
        removable_files=[paired_output1, unpaired_output1, paired_output2, unpaired_output2, single_output]
    )

def trimmomatic16S(
    input1,
    input2,
    paired_output1,
    unpaired_output1,
    paired_output2,
    unpaired_output2,
    single_output,
    quality_offset,
    trim_log,
    headcrop
    ):

    if input2:  # Paired end reads
        inputs = [input1, input2]
        outputs = [paired_output1, unpaired_output1, paired_output2, unpaired_output2]
    else:   # Single end reads
        inputs = [input1]
        outputs = [single_output]

    headcrop_length = str(headcrop)

    return Job(
        inputs,
        outputs + [trim_log],
        [
            ['trimmomatic', 'module_java'],
            ['trimmomatic', 'module_trimmomatic']
        ],
        command="""\
java -XX:ParallelGCThreads=1 -Xmx{ram} -jar $TRIMMOMATIC_JAR {mode} \\
  -threads {threads} \\
  -phred{quality_offset} \\
  {inputs} \\
  {outputs} \\
  HEADCROP:{headcrop_length} \\
  2> {trim_log}""".format(
        ram=global_conf.global_get('trimmomatic16S', 'ram'),
        mode = "PE" if input2 else "SE",
        threads=global_conf.global_get('trimmomatic16S', 'threads', param_type='posint'),
        quality_offset=quality_offset if quality_offset == 64 else "33",
        inputs=" \\\n  ".join(inputs),
        outputs=" \\\n  ".join(outputs),
        headcrop_length=str(headcrop_length),
        trim_log=trim_log
        ),
        removable_files=[paired_output1, unpaired_output1, paired_output2, unpaired_output2, single_output]
    )
