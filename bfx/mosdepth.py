################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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

def mosdepth(
        input,
        output_prefix,
        per_base=False,
        regions=None,
        ini_section='mosdepth'
):
    
    outputs = [
            output_prefix + ".mosdepth.global.dist.txt",
            output_prefix + ".mosdepth.summary.txt",
            output_prefix + ".mosdepth.region.dist.txt"
            ]
    return Job(
            [input],
            outputs,
            [
                [ini_section, 'module_python'],
                [ini_section, 'module_mosdepth']
            ],
            command="""\
mosdepth {other_options} \\
{reference} \\
{per_base} \\
{regions} \\
{chrom} \\
{output_prefix} \\
{input}""".format(
    other_options=config.param(ini_section, 'other_options'),
    reference="--fasta " + config.param(ini_section, 'genome_fasta', param_type='filepath'),
    per_base="--no-per-base " if not per_base else "",
    regions="--by " + regions if regions else "",
    chrom="--chrom " + config.param(ini_section, 'chrom') if config.param(ini_section, 'chrom', required=False) else "",
    output_prefix=output_prefix,
    input=input
        )
    )

def run(
        input,
        output_prefix,
        per_base=False,
        regions=None,
        ini_section='mosdepth'
):
    outputs = [
            output_prefix + ".mosdepth.global.dist.txt",
            output_prefix + ".mosdepth.summary.txt",
            output_prefix + ".quantized.bed.gz"
            ]
    
    if regions is not None:
        outputs.append(output_prefix + ".mosdepth.region.dist.txt")
        
    return Job(
            [input],
            outputs,
            [
                [ini_section, 'module_python'],
                [ini_section, 'module_mosdepth']
            ],
            command="""\
export MOSDEPTH_Q0=NO_COVERAGE && \\
export MOSDEPTH_Q1=LOW_COVERAGE && \\
export MOSDEPTH_Q2=CALLABLE && \\
export MOSDEPTH_Q3=HIGH_COVERAGE && \\
mosdepth {other_options} \\
--quantize 0:{Q1}:{Q2}:{Q3} \\
{reference} \\
{regions} \\
{chrom} \\
{output_prefix} \\
{input}""".format(
    Q1=config.param(ini_section, 'Q1'),
    Q2=config.param(ini_section, 'Q2'),
    Q3=config.param(ini_section, 'Q3'),
    other_options=config.param(ini_section, 'other_options'),
    reference="--fasta " + config.param(ini_section, 'genome_fasta', param_type='filepath'),
    per_base="--no-per-base " if not per_base else "",
    regions="--by " + regions if regions else "",
    chrom="--chrom " + config.param(ini_section, 'chrom') if config.param(ini_section, 'chrom', required=False) else "",
    output_prefix=output_prefix,
    input=input
        )
    )
def parse_dedup_coverage_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export dedup_coverage=`awk '$1=="total" {{print $4}}' {input_file}`"""
        )