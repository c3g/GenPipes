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

def pileup(input_bam, output):
    
    return Job(
        [input_bam],
        [output],
        [
            ['conpair_concordance_contamination', 'module_java'],
            ['conpair_concordance_contamination', 'module_python'],
            ['conpair_concordance_contamination', 'module_gatk'],
            ['conpair_concordance_contamination', 'module_conpair']
        ],
        command="""\
python3 $CONPAIR_SCRIPTS/run_gatk_pileup_for_sample.py -t {tmp_dir} \\
  -m {ram} \\
  -G $GATK_JAR \\
  -D $CONPAIR_DIR \\
  -R {reference_sequence} \\
  -M {markers} \\
  -B {input} \\
  -O {output} {other_options}""".format(
        tmp_dir=config.param('conpair_concordance_contamination', 'tmp_dir'),
        ram=config.param('conpair_concordance_contamination', 'ram'),
        reference_sequence=config.param('conpair_concordance_contamination', 'genome_fasta', param_type='filepath'),
        markers=config.param('conpair_concordance_contamination', 'markers_bed'),
        input=input_bam,
        other_options=" \\\n  " + config.param('conpair_concordance_contamination', 'other_options', required= False) if config.param('conpair_concordance_contamination', 'other_options', required= False) else "",
        output=output
        )
    )

def concordance(input_normal, input_tumor, output):

    return Job(
        [input_normal, input_tumor],
        [output],
        [
            ['conpair_concordance_contamination', 'module_python'],
            ['conpair_concordance_contamination', 'module_conpair']
        ],
        command="""\
python3 $CONPAIR_SCRIPTS/verify_concordance.py {options} \\
  --markers {markers} \\
  --normal_pileup {input_normal} \\
  --tumor_pileup {input_tumor} \\
  --outfile {output}""".format(
        options=config.param('conpair_concordance_contamination', 'concord_options'),
        markers=config.param('conpair_concordance_contamination', 'markers_txt'),
        input_normal=input_normal,
        input_tumor=input_tumor,
        output= output
        )
    )

def contamination(input_normal, input_tumor, output):

    return Job(
        [input_normal, input_tumor],
        [output],
        [
            ['conpair_concordance_contamination', 'module_python'],
            ['conpair_concordance_contamination', 'module_conpair']
        ],
        command="""\
python3 $CONPAIR_SCRIPTS/estimate_tumor_normal_contamination.py {options} \\
  --markers {markers} \\
  --normal_pileup {input_normal} \\
  --tumor_pileup {input_tumor} \\
  --outfile {output}""".format(
        options=config.param('conpair_concordance_contamination', 'contam_options'),
        markers=config.param('conpair_concordance_contamination', 'markers_txt'),
        input_normal=input_normal,
        input_tumor=input_tumor,
        output=output
        )
    )

