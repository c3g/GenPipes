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

# Python Standard Modules
import re

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def index(
        input,
        ini_section='samtools_index'
):
    output = [re.sub(r'\b(bam|cram)\b',
                      lambda match: match.group() + (".bai" if match.group() == "bam" else ".crai"),
                      input)]

    return Job(
        [input],
        output,
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools index \\
  {options} \\
  {input}""".format(
            options=global_conf.global_get(ini_section, 'options'),
            input=input
            )
        )

def faidx(
        input,
        filter=None,
        ini_section='samtools_faidx'
):

    return Job(
        [input],
        [input + ".fai"],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools faidx \\
  {input} \\
  {filter}""".format(
            input=input,
            filter=filter if filter else ""
            )
        )

def flagstat(
        input,
        output,
        ini_section='samtools_flagstat'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools flagstat {other_options} \\
  {threads} \\
  {input} \\
  > {output}""".format(
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            threads="--threads " + global_conf.global_get(ini_section, 'threads'),
            input=input,
            output=output
            )
        )

def mpileup(
        inputs,
        output,
        region=None,
        regionFile=None,
        ini_section='rawmpileup'
):
    if not isinstance(inputs, list):
        inputs = [inputs]

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools mpileup {other_options} \\
  {reference_fasta} \\
  {region} \\
  {regionFile} \\
  {input_bams} \\
  {output}""".format(
            other_options = global_conf.global_get(ini_section, 'mpileup_other_options'),
            reference_fasta = "-f " + global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath') if global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath') else "",
            region = "-r " + region if region else "",
            regionFile = "-l " + regionFile if regionFile else "",
            input_bams = " \\\n  ".join([input_bam for input_bam in inputs]),
            output = "> " + output if output else ""
            )
        )

def merge(
        sample_output,
        input_bams,
        ini_section='hicup_align'
):
    """
    merges an array of bams into a single bam
    """
    postfix = sample_output.split('.')[-1].upper()
    
    return Job(
        input_bams,
        [sample_output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools merge -f --write-index \\
  --output-fmt {output_format} \\
  {threads} \\
  -o {sample_output} \\
  {input_bams}""".format(
            output_format=postfix,
            threads="--threads " + global_conf.global_get(ini_section, 'threads'),
            sample_output=sample_output,
            input_bams="".join([" \\\n  " + input_bam for input_bam in input_bams]),
        )
    )

def bcftools_mpileup(
        inputs,
        output,
        options,
        region=None,
        regionFile=None,
        ini_section='rawmpileup'
):

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_samtools']
        ],
    command="""\
bcftools mpileup {options} \\
  -f {reference_fasta}{region}{regionFile}{inputs}{output}""".format(
        options=options if options else "",
        reference_fasta=global_conf.global_get('samtools_mpileup', 'genome_fasta', param_type='filepath'),
        regionFile=" \\\n  -l " + regionFile if regionFile else "",
        region=" \\\n  -r " + region if region else "",
        inputs="".join([" \\\n  " + input for input in inputs]),
        output=" \\\n  > " + output if output else ""
         )
    )

def sort(
        input,
        output,
        sort_by_name=False,
        ini_section='samtools_sort'
):
    
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools sort \\
  {other_options} {sort_by_name} \\
  {reference} \\
  {tmp_dir} \\
  {output_prefix} \\
  {input_bam}""".format(
            other_options=global_conf.global_get(ini_section, 'other_options', required=False),
            sort_by_name="-n " if sort_by_name else " ",
            tmp_dir="-T " + global_conf.global_get(ini_section, 'tmp_dir'),
            reference="--reference " + global_conf.global_get(ini_section, 'genome_fasta'),
            input_bam=input,
            output_prefix="-o " + output
            ),
        )

def view(
        input,
        output=None,
        options="",
        removable=False,
        ini_section='samtools_view'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools view \\
  {options} \\
  {input}  \\
  {output}""".format(
            options=options,
            input=input,
            output="> " + output if output else ""
            ),
        removable_files=[output if removable else ""]
    )

def fixmate(
        input,
        output=None,
        options="",
        ini_section='samtools_fixmate'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools fixmate \\
  {options} \\
  {input} \\
  {output}""".format(
            options=options,
            input=input,
            output="> " + output if output else "-"
            )
        )

def bcftools_cat(
        inputs,
        output,
        ini_section='samtools_cat'
):

    if not isinstance(inputs, list):
        inputs = [inputs]
        
    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_bcftools']
        ],
        command="""\
bcftools concat \\
  {inputs} \\
  {output}""".format(
            inputs=" \\\n  ".join(inputs),
            output="> " + output if output else ""
            )
        )

def bcftools_view(
        input,
        output,
        options="",
        pair_calling=False,
        ini_section='samtools_view'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_bcftools']
        ],
        command="""\
bcftools view \\
  {pair_calling} {options} \\
  {input} \\
  {output}""".format(
            options=options,
            pair_calling="-T pair" if pair_calling else "",
            input=input,
            output="> " + output if output else ""
            )
        )

def bcftools_call(
        input,
        output,
        options="",
        pair_calling=False,
        ini_section='samtools_call'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_bcftools']
        ],
        command="""\
bcftools call \\
  {pair_calling} {options} \\
  {input} \\
  {output}""".format(
            options=options,
            pair_calling="-T pair" if pair_calling else "",
            input=input,
            output="> " + output if output else ""
            )
        )

def bcftools_call_pair(
        input,
        output,
        options="",
        pair_calling=False,
        ini_section='samtools_paired'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
$BCFTOOLS_BIN/bcftools view \\
  {pair_calling} {options} \\
  {input} \\
  {output}""".format(
            options=options,
            pair_calling="-T pair" if pair_calling else "",
            input=input,
            output="> " + output if output else ""
            )
        )

def bcftools_cat_pair(
        inputs,
        output,
        ini_section='samtools_paired'
):

    if not isinstance(inputs, list):
        inputs = [inputs]

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
$BCFTOOLS_BIN/bcftools cat \\
  {inputs} \\
  {output}""".format(
            inputs=" \\\n  ".join(inputs),
            output="> " + output if output else ""
            )
        )

def bcftools_view_pair(
        input,
        output,
        options="",
        pair_calling=False,
        ini_section='samtools_paired'
):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
$BCFTOOLS_BIN/bcftools view \\
  {pair_calling} {options} \\
  {input} \\
  {output}""".format(
            options=options,
            pair_calling="-T pair" if pair_calling else "",
            input=input,
            output="> " + output if output else ""
            )
        )


def mapped_count(
        bam,
        output=None,
        bed=None,
        options=None,
        ini_section='samtools_count'
):

    return Job(
        [bam],
        [output],
        [
            [ini_section, 'module_samtools']
        ],
        command="""\
samtools view -F4 {options} -c \\
  {input} \\
  {target_option} \\
  {output}""".format(
            options=options if options else " ",
            input=bam,
            output="> " + output if output else " ",
            target_option="-L " + bed if bed else " "
            ),
        )

def bam2fq(
        input_bam,
        output_pair1,
        output_pair2,
        output_other,
        output_single,
        ini_section='samtools_bam2fq'
):
    if output_pair2:  # Paired end reads
        outputs = [output_pair1, output_pair2, output_other, output_single]
    else:   # Single end reads
        outputs = [output_other]

    return Job(
        [input_bam],
        outputs,
        [
            ['samtools', 'module_samtools'],
        ],
        command="""\
samtools bam2fq {other_options} \\
  {nthread} \\
  {output_pair1}{output_pair2}{output_single}{output_other} \\
  {input_bam}""".format(
          other_options=global_conf.global_get(ini_section, 'samtools_bam2fq_other_options', required=False),
          nthread="-@ " + global_conf.global_get(ini_section, 'samtools_bam2fq_threads', required=False),
          output_pair1="-1 " + output_pair1 if output_pair2 else "",
          output_pair2=" -2 " + output_pair2 if output_pair2 else "",
          output_single= " -s " + output_single if output_pair2 else "",
          output_other = " -0 " + output_other,
          input_bam=input_bam
          )
        )

def quickcheck(
        input,
        output=None,
        options=None,
        ini_section='samtools'
):

    return Job(
            [input],
            [output],
            [
                [ini_section, 'module_samtools'],
            ],
            command="""\
samtools quickcheck {options} \\
  {input} {output}""".format(
      options=options if options else "",
      input=input,
      output="> " + output if output else ""
      )
  )

