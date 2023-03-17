################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

# MUGQIC Modules
from core.config import *
from core.job import *

def trimmer(
        input1,
        input2,
        output
        ):
    return Job(
            [input1, input2],
            [output + "_R1.fastq.gz", output + "_R2.fastq.gz"],
            [
                ['agent_trimmer', 'module_java'],
                ['agent_trimmer', 'module_agent']
            ],
            command="""\
java {java_other_options} -Xmx{ram} -jar $TRIMMER \\
  -fq1 {input1} \\
  -fq2 {input2} \\
  -{sample_type} \\
  -out {output_prefix} {other_options}""".format(
      java_other_options=config.param('agent_trimmer', 'java_other_options'),
      ram=config.param('agent_trimmer', 'ram'),
      input1=input1,
      input2=input2,
      sample_type=config.param('agent_trimmer', 'sample_type'),
      output_prefix=output,
      other_options=config.param('agent_trimmer', 'trimmer_other_options', required=False)
      )
  )


def locatit(input,
            output,
            bed=None,
            mode=None):
    return Job(
            [input],
            [output],
            [
                ['agent_locatit', 'module_java'],
                ['agent_locatit', 'module_agent']
            ],

            command="""\
java {java_other_options} -Xmx{ram} -jar $LOCATIT \\
  {mode} \\
  {other_options} \\
  {bed} \\
  -o {output} \\
  {input}""".format(
      java_other_options=config.param('agent_locatit', 'java_other_options'),
      ram=config.param('agent_locatit', 'ram'),
      mode="-"+mode if mode else "",
      other_options=config.param('agent_locatit', 'locatit_other_options'),
      bed="-l " + bed if bed else "",
      output=output,
      input=input
    ) 
)
