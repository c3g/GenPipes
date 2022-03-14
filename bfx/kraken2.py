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
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def kraken2(input1, input2, prefix, other_options=global_config_parser.param('kraken2', 'other_options', required=False), nthread=global_config_parser.param('kraken2', 'threads', required=False), database=global_config_parser.param('kraken2', 'database', required=False)):
    output = prefix + ".kraken2_output"
    report = prefix + ".kraken2_report"

    outputs = [
        output,
        report
        ]

    if input2:  # Paired end reads
        inputs = [input1, input2]
        unclassified_output_1 = prefix + ".unclassified_sequences_1.fastq"
        unclassified_output_2 = prefix + ".unclassified_sequences_2.fastq"
        unclassified_output_paired = prefix + ".unclassified_sequences#.fastq"
        classified_output_1 = prefix + ".classified_sequences_1.fastq"
        classified_output_2 = prefix + ".classified_sequences_2.fastq"
        classified_output_paired = prefix + ".classified_sequences#.fastq"
        outputs.extend((unclassified_output_1, unclassified_output_2, classified_output_1, classified_output_2))
    else:   # Single end reads
        inputs = [input1]
        unclassified_output = prefix + ".unclassified_sequences.fastq"
        classified_output = prefix + ".classified_sequences.fastq"
        outputs.extend((unclassified_output, classified_output))

    return Job(
        inputs,
        outputs,
        [
            ['kraken2', 'module_kraken2']
        ],

        command="""\
kraken2 \\
  {other_options} {paired} \\
  {nthread} \\
  {database} \\
  {unclassified_output} \\
  {classified_output} \\
  {output} \\
  {report} \\
  {inputs}""".format(
      other_options=other_options,
      paired="--paired" if input2 else "",
      nthread="--threads " + nthread,
      database="--db " + database,
      unclassified_output="--unclassified-out " + (unclassified_output_paired if input2 else unclassified_output),
      classified_output="--classified-out " + (classified_output_paired if input2 else classified_output),
      output="--output " + output,
      report="--report " + report,
      inputs=" \\\n  ".join(inputs)
      ),
    )
