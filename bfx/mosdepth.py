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
from core.config import config
from core.job import Job

def mosdepth(input, output_prefix, per_base=False, regions=None):
    outputs = [
            output_prefix + ".mosdepth.global.dist.txt",
            output_prefix + ".mosdepth.summary.txt",
            output_prefix + ".mosdepth.region.dist.txt"
            ]
    return Job(
            [input],
            outputs,
            [
                ['mosdepth', 'module_python'],
                ['mosdepth', 'module_mosdepth']
            ],
            command="""\
mosdepth {other_options} \\
{per_base} \\
{regions} \\
{chrom} \\
{output_prefix} \\
{input}""".format(
    other_options=config.param('mosdepth', 'other_options'),
    per_base="--no-per-base " if not per_base else "",
    regions="--by " + regions if regions else "",
    chrom="--chrom " + config.param('mosdepth', 'chrom') if config.param('mosdepth', 'chrom', required=False) else "",
    output_prefix=output_prefix,
    input=input
        )
    )