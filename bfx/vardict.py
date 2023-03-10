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

# MUGQIC Modules
from core.config import *
from core.job import *

def paired(input_normal, input_tumor, tumor_name, output=None, region=[]):
    return Job(
        [input_normal, input_tumor, region],
        [output],
        [
        ['vardict_paired', 'module_vardict'],
        ['vardict_paired', 'module_samtools'],
        ['vardict_paired', 'module_perl'],
        ['vardict_paired', 'module_R']
        ],
        command="""\
vardict \\
  -G {reference_fasta} \\
  -N {tumor_name} \\
  -b "{paired_samples}" \\
  {vardict_options}{region}{output}""".format(
        reference_fasta=config.param('vardict_paired', 'genome_fasta', param_type='filepath'),
        tumor_name=tumor_name,
        paired_samples=input_tumor + "|" + input_normal,
        vardict_options=config.param('vardict_paired', 'vardict_options'),
        region=" \\\n  " + region if region else "",
        output=" \\\n  > " + output if output else ""
        )
    )

def paired_java(input_normal, input_tumor, tumor_name, output=None, region=[]):
    return Job(
        [input_normal, input_tumor, region],
        [output],
        [
        ['vardict_paired', 'module_java'],
        ['vardict_paired', 'module_vardict_java'],
        ['vardict_paired', 'module_samtools'],
        ['vardict_paired', 'module_perl'],
        ['vardict_paired', 'module_R']
        ],
        command="""\
java {java_other_options} -Djava.io.tmpdir={tmp_dir} -Xms768m -Xmx{ram} -classpath {classpath} \\
  -G {reference_fasta} \\
  -N {tumor_name} \\
  -b "{paired_samples}" \\
  {vardict_options}{region}{output}""".format(
        tmp_dir=config.param('vardict_paired', 'tmp_dir'),
        reference_fasta=config.param('vardict_paired', 'genome_fasta', param_type='filepath'),
        tumor_name=tumor_name,
        paired_samples=input_tumor + "|" + input_normal,
        java_other_options=config.param('vardict_paired', 'java_other_options'),
        ram=config.param('vardict_paired', 'ram'),
        classpath=config.param('vardict_paired', 'classpath'),
        vardict_options=config.param('vardict_paired', 'vardict_options'),
        region=" \\\n  " + region if region else "",
        output=" \\\n  > " + output if output else ""
        )
    )

def single_java(input_bam, sample_name, output=None, nosv=False, freq=None, region=[]):
    return Job(
            [input_bam],
            [output],
            [
            ['vardict_single', 'module_java'],
            ['vardict_single', 'module_vardict_java'],
            ['vardict_single', 'module_perl'],
            ['vardict_single', 'module_R']
            ],
            command="""\
java {java_other_options} -Djava.io.tmpdir={tmp_dir} -Xms768m -Xmx{ram} -classpath {classpath} \\
  -N {sample_name} \\
  -G {reference_fasta} \\
  -b {input_bam} \\
  {nosv} \\
  -f {freq} \\
  {region} \\
  {vardict_options} \\
  -r {min_reads} \\
  -q {min_phred}{output}""".format(
      java_other_options=config.param('vardict_single', 'java_other_options'),
      tmp_dir=config.param('vardict_single', 'tmp_dir'),
      ram=config.param('vardict_single', 'ram'),
      classpath=config.param('vardict_single', 'classpath'),
      sample_name=sample_name,
      reference_fasta=config.param('vardict_single', 'genome_fasta', param_type='filepath'),
      input_bam=input_bam,
      nosv=" \\\n --nosv " if nosv else "",
      freq=freq if freq else "0.01",
      region=" \\\n " + region if region else "",
      vardict_options=config.param('vardict_single', 'vardict_options'),
      min_reads=config.param('vardict_single', 'min_reads'),
      min_phred=config.param('vardict_single', 'min_phred'),
      output=" \\\n > " + output if output else ""
      )
    )


def testsomatic(input=None, output=None):
    return Job(
        [input],
        [output],
        [
        ['vardict_paired', 'module_vardict_java'],
        ['vardict_paired', 'module_R']
        ],
        command="""\
$VARDICT_BIN/testsomatic.R {input} {output}""".format(
        input=" \\\n " + input if input else "",
        output=" \\\n  > " + output if output else ""
        )
    )


def teststrandbias(input=None, output=None):
    return Job(
            [input],
            [output],
            [
           # ['vardict_single', 'module_vardict_java'],
            ['vardict_single', 'module_R']
            ],
            command="""\
$VARDICT_BIN/teststrandbias.R {input} {output}""".format(
            input=" \\\n " + input if input else "",
            output=" \\\n " + output if output else ""
            )
    )


def var2vcf(output, normal_name, tumor_name, input=None):
    return Job(
        [input],
        [output],
        [
        ['vardict_paired', 'module_vardict_java'],
        ['vardict_paired', 'module_perl']
        ],
        command="""\
perl $VARDICT_BIN/var2vcf_paired.pl \\
    -N "{pairNames}" \\
    {var2vcf_options}{input}{output}""".format(
        pairNames=tumor_name + "|" + normal_name,
        var2vcf_options=config.param('vardict_paired', 'var2vcf_options'),
        input=" \\\n " + input if input else "",
        output=" \\\n  > " + output if output else ""
        )
    )


def var2vcf_valid(output, sample_name, freq, input=None):
    return Job(
            [input],
            [output],
            [
           # ['vardict_single', 'module_vardict_java'],
            ['vardict_single', 'module_perl']
            ],
            command="""\
perl $VARDICT_BIN/var2vcf_valid.pl \\
    -N {sample_name} \\
    -f {freq} \\
    {other_options}{input}{output}""".format(
        sample_name=sample_name,
        freq=freq if freq else 0.01,
        other_options=config.param('vardict_single', 'var2vcf_valid_options'),
        input=" \\\n " + input if input else "",
        output=" \\\n > " + output if output else ""
        )
    )


def dict2beds(dictionary,beds):
    return Job(
        [dictionary],
        beds,
        [
            ['vardict_paired', 'module_mugqic_tools'],
            ['vardict_paired', 'module_python']
        ],
        command="""\
dict2BEDs.py \\
  --dict {dictionary} \\
  --beds {beds} {dict2bed_options}""".format(
        dictionary=dictionary if dictionary else config.param('DEFAULT', 'genome_dictionary', param_type='filepath'),
        beds=' '.join(beds),
        dict2bed_options=config.param('vardict_paired', 'dict2bed_options')
        )
    )

