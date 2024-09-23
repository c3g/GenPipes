################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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


# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def verify(
        input_bam, 
        output_prefix,
        ini_section='verify_bam_id'):
    return Job(
        [input_bam],
        [output_prefix + ".selfSM"],
        [
            [ini_section, 'module_verify_bam_id']
        ],
        command="""\
verifyBamID \\
  --vcf {input_vcf} \\
  --bam {input_bam} \\
  --out {output_prefix} \\
  {other_options}""".format(
            input_vcf=global_conf.global_get(ini_section, 'vcf', param_type='filepath'),
            input_bam=input_bam,
            output_prefix=output_prefix,
            other_options=global_conf.global_get(ini_section, 'options')
        )
    )

def verify2(
        input_bam,
        output_prefix,
        bed_file=None,
        ini_section='verify_bam_id2'
):
    if not isinstance(input_bam, list):
        inputs = [input_bam]
        
    if bed_file is not None:
        inputs.append(bed_file)
        
    return Job(
        inputs,
        [output_prefix + ".selfSM"],
        [
            [ini_section, 'module_verify_bam_id']
        ],
        command="""\
VerifyBamID {other_options} \\
  --SVDPrefix {svdprefix} \\
  --Reference {reference} \\
  --BamFile {input_bam} \\
  --Output {output_prefix}{bed_file}""".format(
            svdprefix=global_conf.global_get(ini_section, 'svd_dataset'),
            reference=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
            input_bam=input_bam,
            output_prefix=output_prefix,
            bed_file=" --BedPath " + bed_file if bed_file else "",
            other_options=global_conf.global_get(ini_section, 'options')
        )
    )

def parse_contamination_freemix_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export contamination_freemix=`grep -v "^#" {input_file} | cut -f7`"""
        )