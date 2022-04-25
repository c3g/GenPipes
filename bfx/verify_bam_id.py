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

def getVersion():
    v2_module_name_pattern = r'\/2(\.\d+)*$'
    v1_module_name_pattern = r'\/1(\.\d+)*$'
    v1_custom_name_pattern = r'\/devMaster_20151216$'
    module_name = config.param('verify_bam_id', 'module_verify_bam_id')
    if (re.search(v1_module_name_pattern, module_name)):
        return 1
    elif (re.search(v1_custom_name_pattern, module_name)):
        return 1
    else:
        if(not re.search(v2_module_name_pattern, module_name)):
            log.warning(f"Cannot determine verifyBamID version from module name '{module_name}'. Assuming v2")
        return 2

def verify(
    input_bam,
    output_prefix,
    var=None,
    ref=None
    ):

    v2_module_name_pattern = r'\/2(\.\d+)*$'
    v1_module_name_pattern = r'\/1(\.\d+)*$'
    v1_custom_name_pattern = r'\/devMaster_20151216$'
    module_name = config.param('verify_bam_id', 'module_verify_bam_id')

    if (getVersion() == 1):
        return verify_v1(input_bam, output_prefix, var)
    elif (getVersion() == 2):
        return verify_v2(input_bam, output_prefix, var, ref)
    else:
        log.error("Could not determine VerifyBamID version from module name in ini")

def verify_v1(
    input_bam,
    output_prefix,
    var=None
    ):

    return Job(
        [input_bam],
        [output_prefix + ".selfSM"],
        [
            ['verify_bam_id', 'module_verify_bam_id']
        ],
        command="""\
verifyBamID \\
--vcf {vcf} \\
--bam {input_bam} \\
--out {output_prefix} \\
{other_options}""".format(
            vcf=var if var else config.param('verify_bam_id', 'vcf', param_type='filepath'),
            input_bam=input_bam,
            output_prefix=output_prefix,
            other_options=config.param('verify_bam_id', 'options')
        )
    )

def verify_v2(
    input_bam,
    output_prefix,
    var=None,
    ref=None
    ):

    return Job(
        [input_bam],
        [output_prefix + ".selfSM"],
        [
            ['verify_bam_id', 'module_verify_bam_id']
        ],
        command="""\
VerifyBamID {other_options} \\
--SVDPrefix {svdprefix} \\
--Reference {reference} \\
--BamFile {input_bam} \\
--Output {output_prefix}""".format(
            svdprefix=var if var else config.param('verify_bam_id', 'svd_dataset'),
            reference=ref if ref else config.param('verify_bam_id', 'genome_fasta', type='filepath'),
            input_bam=input_bam,
            output_prefix=output_prefix,
            other_options=config.param('verify_bam_id', 'options', required=False)
            )
        )
