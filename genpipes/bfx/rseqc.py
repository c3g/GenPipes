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
import os
import re

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)

def bam_stat(input, output, ini_section='rseqc'):
    return Job(
        [input],
        [output],
        [
            [ini_section,'module_rseqc']
        ],

    command="""\
    bam_stat.py {options} \\
        -i {input} \\
        > {output}""".format(
            options=global_conf.global_get(ini_section, 'bam_stat_options'),
            input=input,
            output=output,
        ),
    )

def gene_body_coverage(input_bam, output, ini_section='rseqc'):
    return Job(
        [input_bam],
        [output],
        [
            [ini_section, 'module_rseqc']
        ],

        command="""\
        geneBody_coverage.py \\
            -r {geneBody} \\
            -i {input} \\
            -o {output}""".format(
            geneBody=global_conf.global_get(ini_section, 'housekeeping'),
            input=input_bam,
            output=output,
        ),
    )

def infer_experiment(input, output, ini_section='rseqc'):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_rseqc']
        ],

        command="""\
        infer_experiment.py \\
            -r {ref_gene_model} \\
            -i {input} \\
            > {output}""".format(
            ref_gene_model=global_conf.global_get(ini_section, 'ref_gene_model'),
            input=input,
            output=output,
        ),
    )

def inner_distance(input, output, ini_section='rseqc'):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_rseqc']
        ],

        command="""\
        inner_distance.py \\
            -r {ref_gene_model} \\
            -i {input} \\
            -o {output}""".format(
            ref_gene_model=global_conf.global_get(ini_section, 'ref_gene_model'),
            input=input,
            output=output,
        ),
    )

def junction_annotation(input, output, ini_section='rseqc'):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_rseqc']
        ],

        command="""\
        junction_annotation.py \\
            -r {refseq} \\
            -i {input} \\
            -o {output}""".format(
            refseq=global_conf.global_get(ini_section, 'refseq'),
            input=input,
            output=output,
        ),
    )

def junction_saturation(input, output, ini_section='rseqc'):
    output_file=output + ".log"
    
    return Job(
        [input],
        [output_file],
        [
            [ini_section, 'module_rseqc'],
            [ini_section, 'module_R']
        ],

        command="""\
        junction_saturation.py \\
            -r {ref_gene_model} \\
            -i {input} \\
            -o {output} 2> \\
            {output}.log""".format(
            ref_gene_model=global_conf.global_get(ini_section, 'ref_gene_model'),
            input=input,
            output=output,
        ),
    )

def tin(input, output_dir, ini_section='rseqc'):
    file = re.sub(r"\.bam$", ".summary.txt", os.path.basename(input))
    output = os.path.join(output_dir, file)
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_rseqc']
        ],

        command="""\
        tin.py \\
            -r {ref_gene_model} \\
            -i {input}""".format(
            ref_gene_model=global_conf.global_get(ini_section, 'ref_gene_model'),
            input=input
        ),
    )
