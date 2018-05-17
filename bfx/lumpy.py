#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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

from core.config import *
from core.job import *

def lumpyexpress_pair(normal_bam, tumor_bam, output_vcf, spl_normal=[], spl_tumor=[], dis_normal=[], dis_tumor=[]):
    
    return Job(
        [normal_bam, tumor_bam, spl_normal, spl_tumor, dis_normal, dis_tumor],
        [output_vcf],
        [
            ['lumpy_paired_sv_calls', 'module_python'],
            ['lumpy_paired_sv_calls', 'module_lumpy'],
            ['lumpy_paired_sv_calls', 'module_samtools'],
        ],
        command="""\
        lumpyexpress {options} \\
        -B {tumor_bam},{normal_bam} \\
        -o {output_vcf} \\
        {splitter_tumor},{splitter_normal} \\
        {discordant_tumor},{discordant_normal} \\
        -K $LUMPY_SCRIPTS/lumpyexpress.config""".format(
            options=config.param('lumpy_paired_sv_calls','options') if config.param('lumpy_paired_sv_calls','options') else "",
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            output_vcf=output_vcf,
            splitter_tumor="-S " + spl_tumor if spl_tumor else "",
            splitter_normal=spl_normal if spl_normal else "",
            discordant_tumor="-D " + dis_tumor if dis_tumor else "",
            discordant_normal=dis_normal if dis_normal else "",
        )
    )
            
