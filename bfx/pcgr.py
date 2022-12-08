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


def report(
    input_vcf,
    input_cna,
    cpsr_report,
    output_dir,
    tumor_id
    ):
    
    if config.param('report_pcgr', 'module_pcgr').split("/")[2] >= "1":
        call = 'pcgr'
    else:
        call = 'pcgr.py'

    return Job(
        [
            input_vcf,
            input_cna,
        ],
        output_dir,
        [
            ['report_pcgr', 'module_pcgr'],
        ],
        command="""\
{call} {options} \\
    {tumor_type} \\
    {assay} \\
    {tumor_options} \\
    {normal_options} \\
    {mutsig_options} \\
    {tmb_options} \\
    {msi_options} \\
    --input_vcf {input_vcf} \\
    --cpsr_report {cpsr_report} \\
    --input_cna {input_cna} \\
    --pcgr_dir $PCGR_DATA \\
    --output_dir {output_dir} \\
    --genome_assembly {assembly} \\
    --sample_id {tumor_id}""".format(
            call=call,
            options=config.param('report_pcgr', 'options'),
            tumor_type=config.param('report_pcgr', 'tumor_type'),
            assay=config.param('report_pcgr', 'assay'),
            tumor_options=config.param('report_pcgr', 'tumor_options'),
            normal_options=config.param('report_pcgr', 'normal_options'),
            mutsig_options=config.param('report_pcgr', 'mutsig_options'),
            tmb_options=config.param('report_pcgr', 'tmb_options'),
            msi_options=config.param('report_pcgr', 'msi_options'),
            input_vcf=input_vcf,
            cpsr_report=cpsr_report,
            input_cna=input_cna,
            output_dir=output_dir,
            assembly=config.param('report_pcgr', 'assembly'),
            tumor_id=tumor_id
        )
    )

def create_header(output):
            return Job(
                command="""\
`cat > {header} << END
Chromosome\tStart\tEnd\tSegment_Mean
END`""".format(
            header=output,
            )
        )