################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

### Start from here ###

def bamcoverage(input_bam, output_file, strand=None):
    return Job(
        input_files=[input_bam],
        output_files=[output_file],
        module_entries=[['default', 'module_deeptools'] 
        ],
        # need to add option for clusters, 
        # each other_option individually & 
        # no need to have forward and reverse commands, all in one, dependent on the ini. 
        # 
        command="""\
bamCoverage --verbose \\
    --outFileFormat bigwig \\
    --numberOfProcessors {cpu} \\
    --binSize {bs} \\
    --normalizeUsing {nu} \\
    --bam {input_bam} \\
    --outFileName {output_file} {strand}""".format(
            output_file=output_file,
            input_bam=input_bam,
            cpu=global_conf.global_get('wiggle', 'cluster_cpu', required=True), 
            bs=global_conf.global_get('wiggle', 'bin_size', required=True),
            nu=global_conf.global_get('wiggle', 'norm_using', required=True),
            strand="--filterRNAstrand " + strand if strand else "", 
        )
  )
