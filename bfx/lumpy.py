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

from core.config import *
from core.job import *

def lumpyexpress_pair(normal_bam, tumor_bam, output_vcf, spl_normal=[], spl_tumor=[], dis_normal=[], dis_tumor=[]):
    if tumor_bam is not None and spl_tumor is not None and dis_tumor is not None:
        inputs = [normal_bam, tumor_bam, spl_normal, spl_tumor, dis_normal, dis_tumor]
    
    else:
        inputs = [normal_bam, spl_normal, dis_normal]
        
    return Job(
        inputs,
        [output_vcf],
        [
            ['lumpy_paired_sv_calls', 'module_python'],
	        #['lumpy_paired_sv_calls', 'module_samblaster'],
            ['lumpy_paired_sv_calls', 'module_samtools'],
            ['lumpy_paired_sv_calls', 'module_sambamba'],
	    ['lumpy_paired_sv_calls', 'module_lumpy'],
        ],
        command="""\
        lumpyexpress {options} \\
        -B {tumor_bam}{normal_bam} \\
        -o {output_vcf} \\
        -S {splitter_tumor}{splitter_normal} \\
        -D {discordant_tumor}{discordant_normal} \\
        -K $LUMPY_SCRIPTS/lumpyexpress.config""".format(
            options=global_config_parser.param('lumpy_paired_sv_calls', 'options') if global_config_parser.param('lumpy_paired_sv_calls', 'options') else "",
            tumor_bam=tumor_bam + "," if tumor_bam else "",
            normal_bam=normal_bam,
            output_vcf=output_vcf,
            splitter_tumor=spl_tumor + "," if spl_tumor else "",
            splitter_normal=spl_normal if spl_normal else "",
            discordant_tumor=dis_tumor + "," if dis_tumor else "",
            discordant_normal=dis_normal if dis_normal else "",
        )
    )
