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
        reference_fasta=global_config_parser.param('vardict_paired', 'genome_fasta', param_type='filepath'),
        tumor_name=tumor_name,
        paired_samples=input_tumor + "|" + input_normal,
        vardict_options=global_config_parser.param('vardict_paired', 'vardict_options'),
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
        tmp_dir=global_config_parser.param('vardict_paired', 'tmp_dir'),
        reference_fasta=global_config_parser.param('vardict_paired', 'genome_fasta', param_type='filepath'),
        tumor_name=tumor_name,
        paired_samples=input_tumor + "|" + input_normal,
        java_other_options=global_config_parser.param('vardict_paired', 'java_other_options'),
        ram=global_config_parser.param('vardict_paired', 'ram'),
        classpath=global_config_parser.param('vardict_paired', 'classpath'),
        vardict_options=global_config_parser.param('vardict_paired', 'vardict_options'),
        region=" \\\n  " + region if region else "",
        output=" \\\n  > " + output if output else ""
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
        var2vcf_options=global_config_parser.param('vardict_paired', 'var2vcf_options'),
        input=" \\\n " + input if input else "",
        output=" \\\n  > " + output if output else ""
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
        dictionary=dictionary if dictionary else global_config_parser.param('DEFAULT', 'genome_dictionary', param_type='filepath'),
        beds=' '.join(beds),
        dict2bed_options=global_config_parser.param('vardict_paired', 'dict2bed_options')
        )
    )

